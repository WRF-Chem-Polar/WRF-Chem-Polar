! 2026/03 Louis Marelle
!
! Purpose:
!  Calculate the sedimentation velocities for aerosols and sediment the
!  aerosols through the tmospheric column
!
! Remarks:
!  - Only compatible with GOCART for now
!  - Sedimentation velocity calculations are the same as in aer_drydep_simple
!  - Sedimentation to the surface is already done in aer_drydep_simple
!  - Hygroscopic growth is only taken into account for seasalt and sulfate
!
! References:
! ------------------------------------------------------------------------------
!  (1) [Gerber1985] Gerber (1985) - Relative humidity parameterization of the
!      Navy aerosol model (NAM), NRL Rep. 8956. National Research Laboratory,
!      Washington DC. https://apps.dtic.mil/sti/tr/pdf/ADA163209.pdf
!  (2) [DL99] Després, B., & Lagoutière, F. (1999). Un schéma non linéaire
!      anti-dissipatif pour l'équation d'advection linéaire. Comptes Rendus de
!      l'Académie des Sciences-Series I-Mathematics, 328(10), 939-943.

MODULE module_aer_settling_simple

  USE module_model_constants, ONLY: g,piconst

  IMPLICIT NONE

  ! Number of wet-scavenged aerosol species for chem_opt
  INTEGER :: naerspecies
  ! List of wet-scavenged aerosol species for chem_opt
  INTEGER, DIMENSION(:), ALLOCATABLE :: aer_pointers
  ! Aerosol dry diameter for aerosol species in aer_pointers (m)
  REAL, DIMENSION(:), ALLOCATABLE :: aer_dry_diam
  ! Aerosol densities for aerosol species in aer_pointers (kg m-3)
  REAL, DIMENSION(:), ALLOCATABLE :: aer_dry_dens
  ! Deposition fluxes for aerosol species in aer_pointers (kg m-2)
  REAL, DIMENSION(:), ALLOCATABLE :: deposition_flux

  ! chem_opt aerosol model type
  LOGICAL :: aer_is_gocart

CONTAINS

!----------------------------------------------------------------------------

  SUBROUTINE aer_settling_simple(chem, t_phy, rho_phy, dz8w, rh, &
                               dtstep, aer_settling_opt,            &
                               dd_dust, dd_seasalt,                 &
                               ids,ide, jds,jde, kds,kde,           &
                               ims,ime, jms,jme, kms,kme,           &
                               its,ite, jts,jte, kts,kte            )
   ! Purpose:
   !  Calculates sedimentation for aerosols
   !
   ! Arguments:
   !  Input:
   !    t_phy = air temperature (degree celsius)
   !    rho_phy = air density (kg/m3)
   !    dz8w = model vertical level depth (m)
   !    rh = relative humidity (0-1)
   !    dtstep = model time step (s)
   !    aer_settling_opt = aerosol settling option config flag (category)
   !    dd_dust = accumulated dust deposition at the surface (kg m-2)
   !    dd_seasalt = accumulated seasalt deposition at the surface (kg m-2)
   !    ids,ims,its etc. = domain, memory, and tile start and end indices
   !  Input/output:
   !    chem = atmospheric component mixing ratio (µg kg-1 for aerosols)

   USE module_state_description, ONLY: p_sulf, p_seas_1, p_seas_2, p_seas_3, &
    p_seas_4, p_dust_1, p_dust_2, p_dust_3, p_dust_4, p_dust_5, num_chem, &
    param_first_scalar

   !---- Arguments
   REAL, DIMENSION(ims:ime, kms:kme, jms:jme), INTENT(IN) :: t_phy, rho_phy, &
     dz8w, rh
   REAL, DIMENSION(ims:ime, jms:jme), INTENT(INOUT) :: dd_dust, dd_seasalt
   REAL, DIMENSION(ims:ime, kms:kme, jms:jme, num_chem),          &
       INTENT(INOUT) :: chem
   REAL, INTENT(IN) :: dtstep
   INTEGER, INTENT(IN) :: aer_settling_opt
   INTEGER, INTENT(IN) :: ids,ide, jds,jde, kds,kde,    &
                          ims,ime, jms,jme, kms,kme,    &
                          its,ite, jts,jte, kts,kte

   !---- Local variables
   ! Settling velocity (m s-1)
   REAL :: vsettl
   ! Settling velocity in column (m s-1)
   REAL, DIMENSION(kts:kte) :: vsettl_1d
   ! Vertical levels depths in column (m)
   REAL, DIMENSION(kts:kte) :: dz_1d
   ! Aerosol load in column, vertically resolved (kg m-2)
   REAL, DIMENSION(kts:kte) :: aerload_1d
   ! Air density in column (kg m-3)
   REAL, DIMENSION(kts:kte) :: rho_1d
   ! Surface aerosol concentration before deposition (µg kg-1)
   REAL :: aerconc_before
   ! Air temperature (K), air density (g cm-3), relative humidity (0-1)
   REAL :: tk, rho_gcm3, relhum
   ! Aerosol dry and wet densities (g cm-3)
   REAL :: aer_dry_dens_gcm3, aer_wet_dens_gcm3
   ! Aerosol dry and wet diameters (cm)
   REAL :: aer_dry_d_cm, aer_wet_d_cm
   ! Sub-timestep (s)
   REAL :: sub_dtstep
   ! Courant number at surface
   REAL :: cfl_sfc
   ! Conversion factor from µg kg-1 to kg m-2
   REAL :: conv
   ! Floor of time step (s), max timestep satisfying CFL (s)
   INTEGER :: dtstep_floor, dtstep_max
   ! Number of subtimesteps sub_dtstep in dtstep, counter for subtimestep loop
   INTEGER :: n_sub_dtsteps, i_sub_dtstep
   ! x,y,z and chemical species index, chemical species pointer
   INTEGER :: i, j, k, iaer, p_chem
   ! Switches for sedimentation option, set depending on aer_settling_opt
   LOGICAL :: sediment_above_surface
   LOGICAL :: sediment_at_surface

   !---- Parameters
   REAL, PARAMETER :: PI = piconst
   ! Parameters for the wet diameter calculation of Gerber1985
   REAL, PARAMETER :: C1_WET_SEASALT = 0.7674
   REAL, PARAMETER :: C2_WET_SEASALT = 3.079
   REAL, PARAMETER :: C3_WET_SEASALT = 2.572E-11
   REAL, PARAMETER :: C4_WET_SEASALT = -1.424
   REAL, PARAMETER :: C1_WET_SULFATE = 0.4809
   REAL, PARAMETER :: C2_WET_SULFATE = 3.082
   REAL, PARAMETER :: C3_WET_SULFATE = 3.110E-11
   REAL, PARAMETER :: C4_WET_SULFATE = -1.428
   ! Only perform sedimentation for particles larger than this
   REAL, PARAMETER :: AER_DRY_DIAM_CRIT = 0.5E-6 ! m

   !---- Set scheme options depending on aer_settling_opt in namelist
   IF(aer_settling_opt == 3) THEN
     sediment_above_surface = .TRUE.
     sediment_at_surface = .TRUE.
   ELSEIF(aer_settling_opt == 2) THEN
     sediment_above_surface = .TRUE.
     sediment_at_surface = .FALSE.
   ELSE
     CALL wrf_error_fatal('aer_settling_simple: Incorrect value for option aer_settling_opt')
   ENDIF

   !-------- Calculate sedimentation for aerosol species -------
   DO j=jts,jte
   DO i=its,ite
     deposition_flux(:) = 0.0
     ! Conversion factor from µg kg-1 to kg m-2 at surface
     conv = rho_phy(i,kts,j)*dz8w(i,kts,j)*1.E-9
     DO iaer = 1,naerspecies
       p_chem = aer_pointers(iaer)

       IF(p_chem >= PARAM_FIRST_SCALAR .AND. aer_dry_diam(iaer) > AER_DRY_DIAM_CRIT) THEN
         vsettl_1d(:) = 0.0
         dz_1d(:) = dz8w(i,:,j)
         rho_1d(:) = rho_phy(i,:,j)

         !---- Calculate sedimentation velocities vsettl in the column
         DO k=kte,kts,-1
           !-- Initialize physical parameters and convert units
           tk = t_phy(i,k,j) + 273.15
           rho_gcm3 = rho_phy(i,k,j) * 1E-3
           relhum = rh(i,k,j)
           aer_dry_d_cm = aer_dry_diam(iaer) * 1.0E2
           !-- Calculate the aerosol wet diameter (cm)
           ! Calculate wet diameter with Gerber1985 eq. 15 and Table 1
           ! This should only be done for hydrophilic aerosols, sea salt, sulfate
           ! This could be enabled with oc2, bc2, using 'rural' parameters from
           ! Gerber1985
           IF(p_chem == p_seas_1 .OR. p_chem == p_seas_2 .OR. p_chem == p_seas_3 &
              .OR. p_chem == p_seas_4) THEN
             aer_wet_d_cm = ( C1_WET_SEASALT*aer_dry_d_cm**C2_WET_SEASALT &
                           / (C3_WET_SEASALT*aer_dry_d_cm**C4_WET_SEASALT - log(relhum)) &
                           + aer_dry_d_cm**3.0 ) ** (1./3.)
           ELSEIF(p_chem == p_sulf) THEN
             aer_wet_d_cm = ( C1_WET_SULFATE*aer_dry_d_cm**C2_WET_SULFATE &
                           / (C3_WET_SULFATE*aer_dry_d_cm**C4_WET_SULFATE - log(relhum)) &
                           + aer_dry_d_cm**3.0 ) ** (1./3.)
           ELSE
             aer_wet_d_cm = aer_dry_d_cm
           ENDIF
           !-- Calculate the aerosol wet density (g/cm3)
           ! Using the dry density for now
           aer_wet_dens_gcm3 = aer_dry_dens(iaer) * 1.E-3
           !---- Calculate the sedimentation velocity
           CALL calc_vsettl(vsettl, tk, rho_gcm3, aer_wet_d_cm, aer_wet_dens_gcm3)
           vsettl_1d(k) = vsettl
         END DO ! k

         !---- Advect the aerosols from upper to lower levels by vsettl
         !-- Determine the maximum time-step satisying the CFL condition in
         !   the column
         dtstep_max = MAXVAL(dz_1d(:)/vsettl_1d(:))
         ! Calculate the subtimestep needed to stay below dtstep_max
         dtstep_floor = INT(dtstep)
         n_sub_dtsteps = MAX(1, INT(dtstep_floor/dtstep_max)+1)
         sub_dtstep = REAL(dtstep) / REAL(n_sub_dtsteps)
         !-- Perform sedimentation over n_sub_dtsteps subtimesteps
         DO i_sub_dtstep = 1,n_sub_dtsteps
           !- Sediment the aerosol aloft
           ! Cannot advect µg kg-1 directly, because densities and vertical
           ! level depths are not uniform. We convert first to µg m-2 which
           ! can be advected safely.
           ! sulf is in ppmv but can also be transformed this way.
           aerload_1d(:) = chem(i,:,j,p_chem) * dz_1d(:) * rho_1d(:)
           IF(sediment_above_surface) THEN
             CALL advect_aer(aerload_1d=aerload_1d, vsettl_1d=vsettl_1d, &
                      dz_1d=dz_1d, timestep=sub_dtstep, kts=kts, kte=kte)
           ENDIF
           ! Convert back to µg kg-1 and ppmv
           chem(i,:,j,p_chem) = aerload_1d(:) / dz_1d(:) / rho_1d(:)
           !- Sediment the aerosols for deposition to the surface, if needed
           IF(sediment_at_surface) THEN
             cfl_sfc = vsettl_1d(kts)*sub_dtstep/dz_1d(kts)
             aerconc_before = chem(i,kts,j,p_chem)
             chem(i,kts,j,p_chem) = MAX(chem(i,kts,j,p_chem) &
                                    - cfl_sfc*chem(i,kts,j,p_chem), 0.)
             deposition_flux(iaer) = deposition_flux(iaer) &
                                   + (aerconc_before-chem(i,kts,j,p_chem))*conv
           ENDIF
         ENDDO ! i_sub_dtstep

       END IF !aer_dry_diam(iaer) > AER_DRY_DIAM_CRIT

       !-- Output deposition fluxes from sedimentation at the surface
       IF(sediment_at_surface) THEN
         IF(p_chem == p_dust_1 .OR. p_chem == p_dust_2 .OR. p_chem == p_dust_3 &
            .OR. p_chem == p_dust_4 .OR. p_chem == p_dust_5) THEN
           dd_dust(i,j) = dd_dust(i,j) + deposition_flux(iaer)
         ELSEIF(p_chem == p_seas_1 .OR. p_chem == p_seas_2 &
            .OR. p_chem == p_seas_3 .OR. p_chem == p_seas_4) THEN
           dd_seasalt(i,j) = dd_seasalt(i,j) + deposition_flux(iaer)
         ENDIF
       ENDIF

     END DO !iaer

   END DO ! i
   END DO ! j

  END SUBROUTINE aer_settling_simple

!----------------------------------------------------------------------------

  SUBROUTINE calc_vsettl(vsettl, tk, rho_gcm3, aer_wet_d_cm, aer_wet_dens_gcm3)
    ! Calculate the aerosol sedimentation velocity
    !
    ! Arguments:
    !  Input:
    !    tk = air temperature (K)
    !    rho_gcm3 = air density (g cm-3)
    !    aer_wet_d_cm = aerosol wet diameter (cm)
    !    aer_wet_dens_gcm3 = aerosol wet density (g cm-3)
    !  Output:
    !    vsettl = aerosol sedimentation velocity (m s-1)
    !
    IMPLICIT NONE
    !---- Arguments
    REAL, INTENT(OUT) :: vsettl
    REAL, INTENT(IN) :: tk, rho_gcm3, aer_wet_d_cm, aer_wet_dens_gcm3
    !---- Local variables
    ! Air kinematic viscosity (cm2 s-1)
    REAL :: airkinvisc
    ! Mean free path of air molecules (cm)
    REAL :: freepath
    ! Cunningham correction factor for the stokes velocity (dimensionless)
    REAL :: cunningham_fact
    !---- Parameters
    REAL, PARAMETER :: GRAVITY = g*100.0 ! cm s-2
    !
    airkinvisc = (1.83E-4 * (416.16/(tk+120.0)) * ((tk/296.16)**1.5)) &
                 / rho_gcm3
    freepath = 7.39758E-4 * airkinvisc / SQRT(tk)
    vsettl = (GRAVITY*aer_wet_dens_gcm3*aer_wet_d_cm*aer_wet_d_cm) &
             / (18.*airkinvisc*rho_gcm3)
    cunningham_fact = 1.0 + 2.0 * freepath/aer_wet_d_cm &
                      * (1.257 + 0.400*EXP(-0.55*aer_wet_d_cm/freepath) )
    vsettl = vsettl * cunningham_fact * 0.01 ! m s-1
  END SUBROUTINE calc_vsettl

!----------------------------------------------------------------------------

  SUBROUTINE advect_aer(aerload_1d, vsettl_1d, dz_1d, timestep, kts, kte)
    ! Advect aerosols from sedimentation from levels k+1 to k
    !
    ! Arguments:
    !   Input:
    !     vsettl_1d = aerosol sedimentation velocity at each level of 1D column
    !                 (m s-1)
    !     dz_1d = vertical level depth in 1D column (m)
    !     kts, kte = start and end indices for vertical level k
    !     timestep = model time step (s)
    !   Input/Output:
    !     aerload_1d = aerosol load in each level of 1D column (µg m-2)
    !
    !---- Arguments
    REAL, DIMENSION(kts:kte), INTENT(INOUT) ::  aerload_1d
    REAL, DIMENSION(kts:kte), INTENT(IN) :: vsettl_1d, dz_1d
    REAL, INTENT(IN) :: timestep
    INTEGER, INTENT(IN) :: kts, kte
    !---- Local variables
    ! Aerosol load in column at previous time step (µg m-2)
    REAL, DIMENSION(kts:kte+1) :: aerload_prev_1d
    ! Courant number (dimensionless) at level k+1
    REAL :: cfl_above
    ! Aerosol load used for advection (µg m-2)
    REAL :: aerload_prime
    ! k-loop index
    INTEGER :: k
    ! Switch for advection option 2 for DL99 scheme, 1 for upwind scheme.
    ! DL99 is less diffusive and extremely computationally simple, so it
    ! should always be used except for tests.
    INTEGER :: settling_adv_opt = 1
    !
    ! Advect aerosols by vsettl from levels k+1 to k
    aerload_prev_1d(kts:kte) = aerload_1d(kts:kte)
    aerload_prev_1d(kte+1) = 0.0
    DO k=kte,kts,-1
      IF(k < kte) THEN
        cfl_above = vsettl_1d(k+1) * timestep / dz_1d(k+1)
        IF(settling_adv_opt == 1) THEN
          ! Advection scheme of Despres and Lagoutiere (1999)
          aerload_prime = cprime_dl99(cupup=aerload_prev_1d(k+2),cup=aerload_prev_1d(k+1),cdo=aerload_prev_1d(k),cfl=cfl_above)
        ELSEIF(settling_adv_opt == 2) THEN
          ! Traditional upwind scheme
          aerload_prime = aerload_prev_1d(k+1)
        ELSE
          CALL wrf_error_fatal('aer_settling_simple: invalid settling_adv_opt value')
        ENDIF
        aerload_1d(k+1) = MAX(aerload_1d(k+1) - cfl_above*aerload_prime, 0.)
        aerload_1d(k) = MAX(aerload_1d(k) + cfl_above*aerload_prime, 0.)
      ENDIF
    END DO
  END SUBROUTINE advect_aer

!----------------------------------------------------------------------------

  FUNCTION cprime_dl99(cupup,cup,cdo,cfl)
    ! Purpose:
    !  This functionc calculates the upwind concentration for advection,
    !  following the Despres and Lagoutiere (1999) formulation.
    !  Adapted from the WRF-CHIMERE implementation by Sylvain Mailler.
    !
    ! References:
    !  [DL99] Després, B., & Lagoutière, F. (1999). Un schéma non linéaire
    !    anti-dissipatif pour l'équation d'advection linéaire. Comptes Rendus
    !    de l'Académie des Sciences-Series I-Mathematics, 328(10), 939-943.
    !  [Mailler2021] Mailler, S., Pennel, R., Menut, L., and Lachâtre, M.:
    !    Using the Després and Lagoutière (1999) antidiffusive transport
    !    scheme: a promising and novel method against excessive vertical
    !    diffusion in chemistry-transport models, Geosci. Model Dev., 14,
    !    2221–2233, https://doi.org/10.5194/gmd-14-2221-2021, 2021.
    !
    !---- Arguments
    ! Aerosol loads in the i-2, i-1, and i cells relative to the advection
    ! direction (advection direction positive) (here in µg m-2)
    REAL, INTENT(IN) :: cupup, cup, cdo
    ! Courant number (dimensionless)
    REAL, INTENT(IN) :: cfl
    !---- Local variables
    ! Upwind concentration/load used in the advection equation (here in µg m-2)
    REAL :: cprime_dl99
    !
    IF((cfl.LE.0.) .OR. (cfl.GE.1.) .OR. ((cup-cupup)*(cdo-cup).LE.0.)) THEN
      ! Fallback to standard upwind in the limit-cases of cfl
      ! or if the considered point is a concentration extremum
      cprime_dl99 = cup
    ELSE
      ! Standard case
      ! Concentration estimate at the output facet.
      cprime_dl99 = MAX(cup+(1.-cfl)/2.*(MAX(0.,MIN(2.*(cup-cupup)/(cdo-cup)/cfl,2./(1-cfl)))*(cdo-cup)),0.)
    ENDIF
  END FUNCTION cprime_dl99

!----------------------------------------------------------------------------

  SUBROUTINE aer_settling_simple_init(chem_opt)
    ! Purpose
    !  Initialize aer_settling_simple parameters, which depend on chem_opt
    !
    ! Arguments
    !  chem_opt = chemistry-aerosol option from namelist
    !
    ! References
    !   GOCART parameters:
    !   - Chin et al. (2002): Tropospheric Aerosol Optical Thickness from
    !     the GOCART Model and Comparisons with Satellite and Sun Photometer
    !     Measurements, J. Atmos. Sci., 59, 461–483,
    !     https://doi.org/10.1175/1520-0469(2002)059<0461:TAOTFT>2.0.CO;2

    USE module_state_description, ONLY:p_sulf, p_seas_1, p_seas_2, p_seas_3, &
      p_seas_4, p_bc1, p_oc1, p_bc2, p_oc2, p_dust_1, p_dust_2, p_dust_3, &
      p_dust_4, p_dust_5, MOZCART_KPP, GOCART_SIMPLE, T1_MOZCART_KPP, &
      GOCARTRADM2, GOCARTRACM_KPP

    IMPLICIT NONE

    !---- Arguments
    !  Input
    INTEGER, INTENT(IN) :: chem_opt

    !---- Local variables
    INTEGER :: status_id

    !---- Initialize chem mechanism type
    aer_is_gocart = chem_opt == MOZCART_KPP .OR. &
                    chem_opt == T1_MOZCART_KPP .OR. &
                    chem_opt == GOCART_SIMPLE .OR. &
                    chem_opt == GOCARTRADM2 .OR. &
                    chem_opt == GOCARTRACM_KPP

    !-------- Initialize aerosol tracer number for chem_opt --------
    ! This could be recovered from the aerosol scheme data directly
    IF(aer_is_gocart) THEN
      ! MOZCART
      naerspecies = 14
    ELSE
      CALL wrf_error_fatal('aer_settling_simple_init: chem_opt not supported')
    ENDIF

    !---- Allocate aerosol mappings and parameters for dry deposition
    IF (.NOT. ALLOCATED(aer_pointers)) THEN
      ALLOCATE(aer_pointers(naerspecies), stat=status_id)
      IF (status_id /= 0) CALL wrf_error_fatal('aer_settling_simple_init: Could not allocate aer_pointers')
    ENDIF
    IF (.NOT. ALLOCATED(aer_dry_diam)) THEN
      ALLOCATE(aer_dry_diam(naerspecies), stat=status_id)
      IF (status_id /= 0) CALL wrf_error_fatal('aer_settling_simple_init: Could not allocate aer_dry_diam')
    ENDIF
    IF (.NOT. ALLOCATED(aer_dry_dens)) THEN
      ALLOCATE(aer_dry_dens(naerspecies), stat=status_id)
      IF (status_id /= 0) CALL wrf_error_fatal('aer_settling_simple_init: Could not allocate aer_dry_dens')
    ENDIF
    IF (.NOT. ALLOCATED(deposition_flux)) THEN
      ALLOCATE(deposition_flux(naerspecies), stat=status_id)
      IF (status_id /= 0) CALL wrf_error_fatal('aer_settling_simple_init: Could not allocate deposition_flux')
    ENDIF

    !---- Initialize aerosol mappings and parameters for dry deposition
    IF(aer_is_gocart) THEN
      !-- GOCART
      ! Species that experience dry deposition
      aer_pointers(:) = (/ p_sulf, p_seas_1, p_seas_2, p_seas_3, p_seas_4, &
                           p_bc1, p_oc1, p_bc2, p_oc2, p_dust_1, &
                           p_dust_2, p_dust_3, p_dust_4, p_dust_5 /)
      ! Effective aerosol dry diameter, from Chin et al. (2002) except dust
      ! and sea salt from WRF-Chem GOCART parameters
      aer_dry_diam(:) = (/ 0.312E-6, 0.6E-6, 2.0E-6, 6.5E-6, 15.0E-6, &
                           0.078E-6, 0.174E-6, 0.078E-6, 0.174E-6, 1.46E-6, &
                           2.8E-6, 4.8E-6, 9.0E-6, 16.0E-6 /)
      ! Aerosol dry densities, from Chin et al. (2002). Value of 1000 kg/m3
      ! for BC might not be realistic but was left like this for consistency
      ! with the rest of GOCART
      aer_dry_dens(:) = (/ 1700., 2200., 2200., 2200., 2200., &
                           1000., 1800., 1000., 1800., 2600., &
                           2600., 2600., 2600., 2600. /)
    ELSE
      CALL wrf_error_fatal ('aer_settling_simple_init: chem_opt not supported')
    ENDIF

  END SUBROUTINE aer_settling_simple_init

!----------------------------------------------------------------------------

END MODULE module_aer_settling_simple

