! Anderson Da Silva, Louis Marelle, 2026
!
! Purpose:
!  Diagnose ice nucleating particle (INP) number concentrations from the
!  model's transported aerosol populations, as a function of a prescribed
!  vector of activation temperatures, for comparison with in-situ immersion
!  freezing observations (e.g. MOSAiC, MOCCHA, ARTofMelt, Tara Polar Station).
!
!  This is a DIAGNOSTIC module: it does not advect or remove INP, and it does
!  not feed back on the microphysics (offline goal #1). The parent aerosols
!  (dust, sea spray, biological, soot) are already transported and removed by
!  the host aerosol scheme. n_INP(T) is recomputed every time step from the
!  current aerosol and a species-dependent ice-active surface-site density
!  n_s(T) (dust, marine, soot) or a per-particle frozen fraction (biological).
!
!  Architecture (two layers, so GOCART now / MOSAIC later reuse the physics):
!    - Layer A (scheme adapter): convert the raw chem array into standardized
!      per-class aerosol number (# kg-1) and surface area (m2 kg-1). Only the
!      GOCART adapter is implemented here; a MOSAIC-4bin adapter is the next step.
!    - Layer B (scheme-agnostic physics): map {N, SA, T} -> n_INP_species(T).
!
! Outputs (registry 4D arrays, dimensioned i,k,j,num_inp_act_temps; # kg-1):
!    inp_dust, inp_marine, inp_bio, inp_soot, inp_total  (cumulative, immersion)
!  Convert to per-litre for obs comparison externally: n[#/L] = n[#/kg]*rho*1e-3.
!
! Activation temperatures (deg C), warmest first:
!    Tc(k) = inp_act_temp_min - (k-1)*inp_act_temp_step,  k = 1..num_inp_act_temps
!    defaults: min = -5 C, step = 2.5 C, n = 11  ->  -5 .. -30 C
!
! ============================================================================
! !!! COEFFICIENT VERIFICATION REQUIRED BEFORE PRODUCTION RUNS !!!
!   VERIFIED here (published coefficients):
!     - Niemand et al. (2012)  desert-dust immersion n_s
!     - DeMott et al. (2015)   dust, number-based (D15)
!   PLACEHOLDERS — functional form coded, COEFFICIENTS MUST BE SET FROM THE
!   PAPERS by the user (search for "TODO VERIFY"):
!     - Tobo et al. (2019)      Icelandic/high-latitude dust immersion n_s
!     - McCluskey et al. (2018) sea-spray-aerosol immersion n_s
!     - Ullrich et al. (2017)   soot immersion n_s
!     - Hummel et al. (2018)    bacteria / fungal-spore frozen fractions
!   The namelist default dust parameterization is therefore Niemand-2012
!   (inp_dust_param=2); switch to Tobo-2019 (=1) once its coefficients are set.
! ============================================================================
!
! References:
!  (1) Da Silva, A., Marelle, L., Lapere, R., and Raut, J.-C. Modelling Ice
!      Nucleating Particles from High-Latitude Sources to Reproduce Arctic In
!      Situ Concentrations. In preparation, 2026
!  (2) Niemand et al. (2012)  - https://doi.org/10.1175/JAS-D-11-0249.1
!  (3) DeMott et al. (2015)   - https://doi.org/10.5194/acp-15-393-2015
!  (4) Tobo et al. (2019)     - https://doi.org/10.1038/s41561-019-0314-x
!  (5) McCluskey et al. (2018)- https://doi.org/10.1175/JAS-D-17-0155.1
!  (6) Ullrich et al. (2017)  - https://doi.org/10.1175/JAS-D-16-0074.1
!  (7) Hummel et al. (2018)   - https://doi.org/10.5194/acp-18-15437-2018

MODULE module_inp_diag

  USE module_model_constants, ONLY: PICONST

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: inp_diag_driver

  !---- Dust parameterization selector (namelist inp_dust_param)
  INTEGER, PARAMETER :: INP_DUST_TOBO2019    = 1
  INTEGER, PARAMETER :: INP_DUST_NIEMAND2012 = 2
  INTEGER, PARAMETER :: INP_DUST_DEMOTT2015  = 3
  INTEGER, PARAMETER :: INP_DUST_ULLRICH2017 = 4

  !---- Activation temperatures (deg C), warmest first. The cumulative immersion
  !     spectrum is evaluated at each of these. This is an analysis choice, not
  !     a runtime option: edit the vector to match the observation report
  !     temperatures. NUM_INP_ACT_TEMPS MUST equal the 'inpt' dimension constant
  !     in Registry/registry.dimspec.
  INTEGER, PARAMETER :: NUM_INP_ACT_TEMPS = 11
  REAL, PARAMETER :: INP_ACT_TEMP_C(NUM_INP_ACT_TEMPS) =                  &
       (/ -5.0, -7.5, -10.0, -12.5, -15.0, -17.5, -20.0,                 &
          -22.5, -25.0, -27.5, -30.0 /)

  !---- GOCART size bins: dry log-centre radii (m), copied from the NWFA/NIFA
  !     construction in phys/module_microphysics_driver.F so that the INP
  !     diagnostic uses exactly the same assumed sizes as the model NIFA.
  REAL, PARAMETER :: rcen_dust(5) = &
       (/ 0.3162E-6, 1.342E-6, 2.324E-6, 4.243E-6, 7.746E-6 /)
  REAL, PARAMETER :: rcen_seas(4) = &
       (/ 0.2236E-6, 0.8660E-6, 2.739E-6, 7.071E-6 /)

  !---- Soot (GOCART bc1+bc2 are bulk mass): assumed monodisperse size/density.
  !     TODO VERIFY: pick a value consistent with the BC size assumed elsewhere.
  REAL, PARAMETER :: rcen_soot   = 0.05E-6   ! m  (0.1 um diameter, placeholder)
  REAL, PARAMETER :: dens_soot   = 1800.0    ! kg m-3

CONTAINS

!------------------------------------------------------------------------------
  SUBROUTINE inp_diag_driver( chem, t_phy, rho_phy, chem_opt,            &
                  inp_dust_param,                                        &
                  inp_dust, inp_marine, inp_bio, inp_soot, inp_total,    &
                  ids,ide, jds,jde, kds,kde,                            &
                  ims,ime, jms,jme, kms,kme,                            &
                  its,ite, jts,jte, kts,kte                             )

    ! Purpose:
    !  Fill the INP diagnostic 4D arrays from the current aerosol field.
    !
    ! Arguments:
    !  chem        = advected chem species (aerosols ug/kg, numbers #/kg)
    !  t_phy       = air temperature (K), unused for n_s (evaluated at the
    !                prescribed activation temperatures, not the local T), kept
    !                for future use (e.g. capping above 0 C)
    !  rho_phy     = air density (kg m-3)
    !  inp_dust_param      = dust n_s scheme selector (see INP_DUST_* params)
    !  inp_dust/marine/bio/soot/total = OUTPUT n_INP (# kg-1), cumulative,
    !                immersion mode, dimensioned (i,k,j,NUM_INP_ACT_TEMPS).
    !                Activation temperatures are the module constant INP_ACT_TEMP_C.

    USE module_state_description, ONLY: num_chem, param_first_scalar,      &
         p_dust_1, p_dust_2, p_dust_3, p_dust_4, p_dust_5,                 &
         p_seas_1, p_seas_2, p_seas_3, p_seas_4,                           &
         p_bc1, p_bc2, p_bact1, p_spore1,                                  &
         MOZCART_KPP, T1_MOZCART_KPP, GOCART_SIMPLE, GOCARTRADM2,          &
         GOCARTRACM_KPP

    IMPLICIT NONE

    !---- Arguments
    INTEGER, INTENT(IN) :: chem_opt, inp_dust_param
    INTEGER, INTENT(IN) :: ids,ide, jds,jde, kds,kde,                      &
                           ims,ime, jms,jme, kms,kme,                      &
                           its,ite, jts,jte, kts,kte
    REAL, DIMENSION(ims:ime, kms:kme, jms:jme, num_chem), INTENT(IN) :: chem
    REAL, DIMENSION(ims:ime, kms:kme, jms:jme), INTENT(IN) :: t_phy, rho_phy
    REAL, DIMENSION(ims:ime, kms:kme, jms:jme, NUM_INP_ACT_TEMPS),         &
         INTENT(INOUT) :: inp_dust, inp_marine, inp_bio, inp_soot, inp_total

    !---- Local
    INTEGER :: i, j, k, it
    REAL    :: tc                          ! activation temperature (deg C)
    REAL    :: ndust(5), sadust(5)         ! # kg-1, m2 kg-1
    REAL    :: nseas(4),  saseas(4)        ! # kg-1, m2 kg-1
    REAL    :: nbact, nspore               ! # kg-1
    REAL    :: nsoot, sasoot               ! # kg-1, m2 kg-1
    REAL    :: sadust_tot, saseas_tot, ndust_tot
    REAL    :: ns_d, ns_m, ns_s            ! surface site densities (m-2)
    REAL    :: n05_cm3                     ! dust number > 0.5 um (cm-3), D15
    LOGICAL :: is_gocart

    !---- Initialize outputs
    inp_dust   = 0.0
    inp_marine = 0.0
    inp_bio    = 0.0
    inp_soot   = 0.0
    inp_total  = 0.0

    ! GOCART-only for now (MOSAIC adapter is the next implementation step)
    is_gocart = chem_opt == MOZCART_KPP    .OR. &
                chem_opt == T1_MOZCART_KPP .OR. &
                chem_opt == GOCART_SIMPLE  .OR. &
                chem_opt == GOCARTRADM2    .OR. &
                chem_opt == GOCARTRACM_KPP
    IF (.NOT. is_gocart) RETURN

    DO j = jts, jte
    DO k = kts, kte
    DO i = its, ite

      !-------- Layer A: GOCART adapter -> number and surface area per class
      ! Dust (5 bins): N = mass / (mass per particle), SA = N * 4 pi r^2
      ndust  = 0.0 ; sadust = 0.0
      CALL gocart_mass_to_num_sa(chem(i,k,j,p_dust_1), rcen_dust(1), 2650.0, ndust(1), sadust(1))
      CALL gocart_mass_to_num_sa(chem(i,k,j,p_dust_2), rcen_dust(2), 2650.0, ndust(2), sadust(2))
      CALL gocart_mass_to_num_sa(chem(i,k,j,p_dust_3), rcen_dust(3), 2650.0, ndust(3), sadust(3))
      CALL gocart_mass_to_num_sa(chem(i,k,j,p_dust_4), rcen_dust(4), 2650.0, ndust(4), sadust(4))
      CALL gocart_mass_to_num_sa(chem(i,k,j,p_dust_5), rcen_dust(5), 2650.0, ndust(5), sadust(5))
      ! NOTE: 2650 kg m-3 is the GOCART dust density; TODO VERIFY against
      ! module_data_gocart_dust den_dust(:) (per-bin values) and replace.

      ! Sea spray (4 bins), den_seas = 2200 kg m-3 (module_data_gocart_seas)
      nseas = 0.0 ; saseas = 0.0
      CALL gocart_mass_to_num_sa(chem(i,k,j,p_seas_1), rcen_seas(1), 2200.0, nseas(1), saseas(1))
      CALL gocart_mass_to_num_sa(chem(i,k,j,p_seas_2), rcen_seas(2), 2200.0, nseas(2), saseas(2))
      CALL gocart_mass_to_num_sa(chem(i,k,j,p_seas_3), rcen_seas(3), 2200.0, nseas(3), saseas(3))
      CALL gocart_mass_to_num_sa(chem(i,k,j,p_seas_4), rcen_seas(4), 2200.0, nseas(4), saseas(4))

      ! Biological: number tracers (# kg-1), zero if not in this chem_opt
      nbact  = 0.0 ; nspore = 0.0
      IF (p_bact1  >= param_first_scalar) nbact  = MAX(chem(i,k,j,p_bact1),  0.0)
      IF (p_spore1 >= param_first_scalar) nspore = MAX(chem(i,k,j,p_spore1), 0.0)

      ! Soot: GOCART bc1+bc2 bulk mass -> assumed monodisperse number/SA
      nsoot = 0.0 ; sasoot = 0.0
      IF (p_bc1 >= param_first_scalar .AND. p_bc2 >= param_first_scalar) THEN
        CALL gocart_mass_to_num_sa( chem(i,k,j,p_bc1)+chem(i,k,j,p_bc2),    &
                                    rcen_soot, dens_soot, nsoot, sasoot )
      END IF

      sadust_tot = SUM(sadust)
      saseas_tot = SUM(saseas)
      ndust_tot  = SUM(ndust)
      ! Dust number > 0.5 um diameter for DeMott-2015: all GOCART dust bins
      ! have centre diameter > 0.5 um, so n05 = total dust number, in cm-3
      n05_cm3 = ndust_tot * rho_phy(i,k,j) * 1.0E-6

      !-------- Layer B: physics, looped over activation temperatures
      DO it = 1, NUM_INP_ACT_TEMPS
        tc = INP_ACT_TEMP_C(it)

        ! Dust
        IF (inp_dust_param == INP_DUST_DEMOTT2015) THEN
          inp_dust(i,k,j,it) = ninp_demott2015_dust(n05_cm3, tc+273.15)    &
                               / MAX(rho_phy(i,k,j)*1.0E-3, 1.0E-12)
          ! D15 gives # per standard litre; convert L-1 -> kg-1 (approx)
        ELSE
          ns_d = ns_dust(tc, inp_dust_param)
          inp_dust(i,k,j,it) = sadust_tot * ns_d
        END IF

        ! Marine (sea spray)
        ns_m = ns_marine(tc)
        inp_marine(i,k,j,it) = saseas_tot * ns_m

        ! Soot
        ns_s = ns_soot(tc)
        inp_soot(i,k,j,it) = sasoot * ns_s

        ! Biological (per-particle frozen fraction)
        inp_bio(i,k,j,it) = nbact*ff_bact(tc) + nspore*ff_spore(tc)

        inp_total(i,k,j,it) = inp_dust(i,k,j,it) + inp_marine(i,k,j,it)     &
                            + inp_bio(i,k,j,it)  + inp_soot(i,k,j,it)
      END DO

    END DO
    END DO
    END DO

  END SUBROUTINE inp_diag_driver

!------------------------------------------------------------------------------
  SUBROUTINE gocart_mass_to_num_sa(mass_ugkg, rcen, dens, num, sa)
    ! Convert a bulk/bin aerosol mass mixing ratio (ug kg-1) to number
    ! (# kg-1) and surface area (m2 kg-1), assuming all mass sits at the bin
    ! centre radius rcen (m) with density dens (kg m-3). This mirrors the
    ! GOCART branch of the NWFA/NIFA construction in module_microphysics_driver.
    IMPLICIT NONE
    REAL, INTENT(IN)  :: mass_ugkg, rcen, dens
    REAL, INTENT(OUT) :: num, sa
    REAL :: mass_per_particle_ug

    ! mass per particle (ug) = 1e9 * (4/3) pi r^3 rho
    mass_per_particle_ug = 1.0E9 * (4.0/3.0)*PICONST * dens * rcen**3
    num = MAX(mass_ugkg, 0.0) / MAX(mass_per_particle_ug, 1.0E-30)
    sa  = num * 4.0*PICONST * rcen**2
  END SUBROUTINE gocart_mass_to_num_sa

!==============================================================================
! Layer B: ice-active surface-site densities n_s(T) [m-2] and bio frozen
! fractions. Tc is the activation temperature in deg C.
!==============================================================================

  REAL FUNCTION ns_dust(tc, param)
    IMPLICIT NONE
    REAL,    INTENT(IN) :: tc
    INTEGER, INTENT(IN) :: param
    REAL :: tcl

    SELECT CASE (param)

    CASE (INP_DUST_NIEMAND2012)
      ! VERIFIED: Niemand et al. (2012), desert-dust immersion.
      ! n_s = exp(-0.517*Tc + 8.934) m-2, valid -36 <= Tc <= -12 C.
      IF (tc > -12.0) THEN
        ns_dust = 0.0
      ELSE
        tcl = MAX(tc, -36.0)
        ns_dust = EXP(-0.517*tcl + 8.934)
      END IF

    CASE (INP_DUST_TOBO2019)
      ! TODO VERIFY: Tobo et al. (2019) Icelandic/high-latitude dust immersion.
      ! Functional form exp(a*Tc + b); SET a,b FROM THE PAPER. Placeholder below
      ! reuses the Niemand slope with a +1 decade offset (NOT the paper values).
      IF (tc > -7.0) THEN
        ns_dust = 0.0
      ELSE
        tcl = MAX(tc, -35.0)
        ns_dust = EXP(-0.517*tcl + 8.934 + 2.3026)   !!! TODO VERIFY coefficients
      END IF

    CASE (INP_DUST_ULLRICH2017)
      ! TODO VERIFY: Ullrich et al. (2017) dust immersion n_s polynomial.
      ns_dust = 0.0                                   !!! TODO VERIFY: implement

    CASE DEFAULT
      ns_dust = 0.0
    END SELECT
  END FUNCTION ns_dust

  REAL FUNCTION ns_marine(tc)
    ! TODO VERIFY: McCluskey et al. (2018) sea-spray-aerosol immersion n_s.
    ! Form n_s = exp(a*Tc + b) m-2; SET a,b FROM THE PAPER.
    IMPLICIT NONE
    REAL, INTENT(IN) :: tc
    REAL :: tcl
    IF (tc > -5.0) THEN
      ns_marine = 0.0
    ELSE
      tcl = MAX(tc, -30.0)
      ns_marine = EXP(-0.545*tcl - 1.2)               !!! TODO VERIFY coefficients
    END IF
  END FUNCTION ns_marine

  REAL FUNCTION ns_soot(tc)
    ! TODO VERIFY: Ullrich et al. (2017) soot immersion n_s. Minor INP source at
    ! mixed-phase temperatures; left at zero until coefficients are set.
    IMPLICIT NONE
    REAL, INTENT(IN) :: tc
    ns_soot = 0.0                                     !!! TODO VERIFY: implement
  END FUNCTION ns_soot

  REAL FUNCTION ff_bact(tc)
    ! TODO VERIFY: Hummel et al. (2018) bacterial frozen fraction f(T) in [0,1].
    ! Bacteria carry the warm tail (active up to ~ -5 C for ice-active strains).
    IMPLICIT NONE
    REAL, INTENT(IN) :: tc
    ff_bact = 0.0                                     !!! TODO VERIFY: implement
    ! Placeholder smooth ramp (NOT paper values), so the field is non-trivial:
    IF (tc <= -2.0) ff_bact = MIN(MAX((-2.0 - tc)/18.0, 0.0), 1.0)
  END FUNCTION ff_bact

  REAL FUNCTION ff_spore(tc)
    ! TODO VERIFY: Hummel et al. (2018) fungal-spore frozen fraction f(T) in [0,1].
    ! Spores activate colder than bacteria (~ -15 to -25 C).
    IMPLICIT NONE
    REAL, INTENT(IN) :: tc
    ff_spore = 0.0                                    !!! TODO VERIFY: implement
    IF (tc <= -12.0) ff_spore = MIN(MAX((-12.0 - tc)/15.0, 0.0), 1.0)
  END FUNCTION ff_spore

  REAL FUNCTION ninp_demott2015_dust(n05_cm3, tk)
    ! VERIFIED form: DeMott et al. (2015) D15, dust.
    !   n_INP[std L-1] = cf * (n05)^(a*(273.16-Tk)+b) * exp(c*(273.16-Tk)+d)
    !   dust: cf=3, a=0, b=1.25, c=0.46, d=-11.6 ; n05 in cm-3
    IMPLICIT NONE
    REAL, INTENT(IN) :: n05_cm3, tk
    REAL, PARAMETER :: cf=3.0, a=0.0, b=1.25, c=0.46, d=-11.6
    REAL :: dt
    IF (n05_cm3 <= 0.0 .OR. tk >= 273.16) THEN
      ninp_demott2015_dust = 0.0
      RETURN
    END IF
    dt = 273.16 - tk
    ninp_demott2015_dust = cf * (n05_cm3**(a*dt + b)) * EXP(c*dt + d)
  END FUNCTION ninp_demott2015_dust

END MODULE module_inp_diag
