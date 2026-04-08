! 2022/12/09 Louis Marelle for LATMOS and CICERO
!
! Purpose:
! Perform cloud-phase aqueous chemistry
!
! Remarks:
!  - Only enabled for MOZCART for now
!  - For GOCART, cloud chemistry can also be performed in gocart_chem but the
!    approach there is too simple and strongly overestimates the speed of so2
!    oxidation
!
! References:
! ==============================================================================
! (1) Chin, M., Ginoux, P., Kinne, S., Torres, O., Holben, B. N., Duncan, B.
!     N., Martin, R. V., Logan, J. A., Higurashi, A., & Nakajima, T. (2002).
!     Tropospheric Aerosol Optical Thickness from the GOCART Model and
!     Comparisons with Satellite and Sun Photometer Measurements, Journal of
!     the Atmospheric Sciences, 59(3), 461-483
! (2) Andre K. Pattantyus, Steven Businger, Steven G. Howell, Review of sulfur
!     dioxide to sulfate aerosol chemistry at Kīlauea Volcano, Hawaii,
!     Atmospheric Environment, Volume 185, 2018, Pages 262-271, ISSN
!     1352-2310, https://doi.org/10.1016/j.atmosenv.2018.04.055.
! (3) Luo, G., Yu, F., and Moch, J. M.: Further improvement of wet process
!     treatments in GEOS-Chem v12.6.0: impact on global distributions of
!     aerosols and aerosol precursors, Geosci. Model Dev., 13, 2879–2903,
!     https://doi.org/10.5194/gmd-13-2879-2020, 2020.
!
!TODO use double precision everywhere?

MODULE module_aer_aqchem_simple

   ! USE module_wrf_error, ONLY : wrf_debug

   IMPLICIT NONE

   PRIVATE

   PUBLIC :: aer_aqchem_simple, aer_aqchem_simple_init

   !---- Global module-wide variables
   !-- Allocatable size-resolved aerosol arrays
   ! Dim is size bin number and depends on aerosol mechanism.
   ! These are initialized in aer_aqchem_simple_init
   REAL, DIMENSION(:), ALLOCATABLE :: bin_diameters_cen
   ! Number of sulfate aerosol bins/modes for chem_opt in use
   INTEGER :: naersizes
   ! chem_opt aerosol model family
   LOGICAL :: aer_is_gocart
   LOGICAL :: aer_is_mosaic4bin
   LOGICAL :: aer_is_mosaic8bin
   !-- Parameters
   ! H+ concentration in aqueous phase (mol L-1)
   !TODO should solve for pH to calculate cloudborne_h, here I use pH =
   ! 4.5 which probably overestimates reaction rates
   REAL, PARAMETER :: cloudborne_h = 3.16E-5



CONTAINS

!----------------------------------------------------------------------------

   SUBROUTINE aer_aqchem_simple(                     &
              id, dtstep, config_flags,              &
              p_phy, t_phy, rho_phy,                 &
              moist, chem, cldfra,                   &
              ids,ide, jds,jde, kds,kde,             &
              ims,ime, jms,jme, kms,kme,             &
              its,ite, jts,jte, kts,kte )

   ! Purpose:
   ! Calculates and updates cloud chemistry tendencies for aerosols
   !
   ! Arguments:
   ! Input:
   !   id = domain ID
   !   dtstep = model time step (s)
   !   config_flags = namelist configuration flags
   !   p_phy = 3D pressure (Pa)
   !   t_phy = 3D physical temperature (K)
   !   rho_phy = 3D air density (kg/m3)
   !   cldfra = 3D cloud fraction (unitless, 0-1)
   !   moist = advected water tracers (kg/kg dry air)
   !   ids,ims,its etc. = domain, memory, and tile start and end indices
   ! Input/output:
   !   chem = advected chemical species (gases in ppm, aerosols in µg/kg)
   !
   ! Remarks:
   !  - Only done for MOZCART for now
   !  - I assume here that all aqueous species are released to the gas phase
   !    at the end of the call, but in reality it only happens by desorption
   !    or evaporation of cloud droplets

   ! USE associations:
   !TODO Use ONLY:
   USE module_configure
   USE module_state_description
   ! molar weight of dry air in g/mol, water density in kg/m3 = g/L
   USE module_model_constants, ONLY: mwdry,rhowater

   !---- Arguments
   INTEGER, INTENT(IN) :: ids,ide, jds,jde, kds,kde,    &
                          ims,ime, jms,jme, kms,kme,    &
                          its,ite, jts,jte, kts,kte,    &
                          id
   REAL, INTENT(IN) :: dtstep
   TYPE(grid_config_rec_type), INTENT(IN) :: config_flags
   REAL, DIMENSION( ims:ime, kms:kme, jms:jme, num_moist ),         &
         INTENT(IN ) ::                                   moist
   REAL, DIMENSION(ims:ime, kms:kme, jms:jme)         ,        &
       INTENT(IN) :: t_phy, p_phy, rho_phy
   REAL, DIMENSION(ims:ime, kms:kme, jms:jme)         ,        &
       INTENT(IN) :: cldfra
   REAL, DIMENSION(ims:ime, kms:kme, jms:jme, num_chem),          &
       INTENT(INOUT) :: chem

   !---- Local variables and parameters
   INTEGER :: i, j, k
   ! Cloud fraction, temperature (K), cloud+rainwater mmr (kg/kg dry air), dry
   ! air density (kg/m3), fraction of condensed water in liquid phase (0-1)
   REAL :: fc, tk, condensed_water, rhoair, liquid_fraction
   ! Concentrations of so2 & h2o2 in cloud water (mol/L), production of
   ! sulfate by cloud chemistry  (mol/L)
   REAL :: cloudborne_so2, cloudborne_h2o2, sulf_aq_prod
   ! Conversion factor from ppm_air to mol/L (aqueous phase)
   REAL :: conv
   ! in-cloud liquid water content (m3 water/m3 dry air), fraction of dissolved so2 &
   ! h2o2 (unitless)
   REAL :: lwc, fw_so2, fw_h2o2
   ! Reaction rate for h2o2+so2 (mol-2 L^2 s-1)
   ! https://doi.org/10.1016/j.atmosenv.2018.04.055
   ! dsulf_aq/dt = k_so2_h2o2[H+][H2O2][so2]_aq
   REAL, PARAMETER :: k_so2_h2o2 = 9.1E+7

   !---- Cloud chemistry code
   DO j=jts,jte
     DO k=kts,kte
       DO i=its,ite
         ! Initialize variables
         fc = cldfra(i,k,j)
         tk = t_phy(i,k,j)
         ! Calculate condensed water for aqueous chemistry, assuming it only
         ! occurs in qc or qr
         condensed_water = moist(i, k, j, p_qc) + moist(i, k, j, p_qr)

         ! Only perform cloud chemistry in cloudy grid cells and above 237K,
         IF (fc > 0.001 .AND. condensed_water > 1.0E-12 .AND. tk > 237.0) THEN
           rhoair = rho_phy(i, k, j)
           IF(tk > 258.0) THEN
             liquid_fraction = 1.0
           ELSE
             liquid_fraction = condensed_water / (condensed_water &
                               + moist(i, k, j, p_qi) + moist(i, k, j, p_qs))
           ENDIF

           ! Calculate conversion factor from ppm_air to in-cloud (or rain) mol/L_aq
           !TODO The true in-cloud mixing ratio is condensed_water/fc
           conv = 1.0E-6 / mwdry / condensed_water / fc * rhowater

           !---- Calculate so2 and h2o2 concencentrations in cloud+rain water in mol/L
           ! ! Assuming instantaneous total dilution - disabled
           ! cloudborne_so2 = chem(i, k, j, p_so2) * conv * fc
           ! cloudborne_h2o2 = chem(i, k, j, p_h2o2) * conv * fc
           ! Calculate the fraction of dissolved gases
           lwc = condensed_water / fc / rhowater * rhoair
           CALL dissolved_gas_frac(tk, lwc, fw_so2, fw_h2o2)
           ! Convert species to mol/L in aqueous phase
           cloudborne_so2 =  liquid_fraction * fw_so2 * chem(i, k, j, p_so2) * conv
           cloudborne_h2o2 =  liquid_fraction * fw_h2o2 * chem(i, k, j, p_h2o2) * conv

           !---- Simple cloud oxidation of dissolved so2 by h2o2
           IF ( cloudborne_so2 > 1.0E-30 .AND. cloudborne_h2o2 > 1.0E-30) THEN
             ! dsulf_aq/dt = k[H+][H2O2][so2]_aq
             sulf_aq_prod = k_so2_h2o2 * cloudborne_h * cloudborne_h2o2 * cloudborne_so2 * dtstep
             ! Do not deplete more than 100% of reactants
             sulf_aq_prod = MIN(sulf_aq_prod, cloudborne_h2o2, cloudborne_so2)
             ! Release species from cloud phase to insterstitial phase, and
             ! convert units from in-cloud mol/L to ppm in air
             chem(i,k,j,p_sulf) = chem(i,k,j,p_sulf) + sulf_aq_prod / conv
             chem(i,k,j,p_h2o2) = chem(i,k,j,p_h2o2) - sulf_aq_prod / conv
             chem(i,k,j,p_so2) = chem(i,k,j,p_so2) - sulf_aq_prod / conv
             ! Prevent very low values
             chem(i,k,j,p_sulf) = MAX(1.0E-32, chem(i,k,j,p_sulf))
             chem(i,k,j,p_so2) = MAX(1.0E-32, chem(i,k,j,p_so2))
             chem(i,k,j,p_h2o2) = MAX(1.0E-32, chem(i,k,j,p_h2o2))

           END IF ! cloudborne_so2 > threshold
        END IF ! fc > fc_threshold .AND. condensed_water > q_threshold

       END DO ! i
     END DO ! k
   END DO ! j

   END SUBROUTINE aer_aqchem_simple

!----------------------------------------------------------------------------

   SUBROUTINE dissolved_gas_frac(tk, lwc, fw_so2, fw_h2o2)

     ! Purpose:
     ! Calculates fw, the fraction of dissolved so2 and h2o2 gas in aqueous
     ! phase (unitless), following Luo et al. (2020)
     !
     ! Arguments:
     ! Input:
     !   tk = air temperature (K)
     !   lwc = in-cloud liquid water content (m3 water / m3 air)
     ! Output:
     !   fw_so2 = fraction of dissolved so2 (unitless)
     !   fw_so2 = fraction of dissolved h2o2 (unitless)

     IMPLICIT NONE

     ! Arguments
     REAL, INTENT(IN) :: tk
     REAL, INTENT(IN) :: lwc
     REAL, INTENT(OUT) :: fw_so2
     REAL, INTENT(OUT) :: fw_h2o2

     ! Local variables and parameters
     ! Henry's law constants H and effective Henry's law constants H*
     ! (mol L-1 atm-1) for so2 and h2o2
     REAL :: h_so2,h_star_so2,h_h2o2,h_star_h2o2
     ! Rate coefficients for H* calculation (mol L-1)
     REAL :: k1,k2,k3
     REAL :: tempfactor
     ! gas constant in L atm K-1 mol -1
     REAL, PARAMETER :: r_h = 0.08205

     !---- Calculate effective henry's law constants and dissolved fractions
     !for so2 and h2o2
     tempfactor = 298.15/tk - 1.0

     ! dissolved so2
     h_so2 = 1.22 * EXP(10.55 * tempfactor)
     k1 = 1.3E-2 * EXP(6.75 * tempfactor)
     k2 = 6.31E-8 * EXP(5.05 * tempfactor)
     h_star_so2 = h_so2 * (1.0 + k1/cloudborne_h + k1*k2/(cloudborne_h*cloudborne_h))
     fw_so2 = 1.0 - 1.0/(1.0 + h_star_so2 * r_h * tk * lwc)
     fw_so2 = MIN(1.0, MAX(0.0, fw_so2))

     ! dissolved h2o2
     h_h2o2 = 8.3E+4 * (24.82 * tempfactor)
     k3 = 2.2E-12 * EXP(12.52 * tempfactor)
     h_star_h2o2 = h_h2o2 * (1.0 + k3/cloudborne_h)
     fw_h2o2 = 1.0 - 1.0/(1.0 + h_star_h2o2 * r_h * tk * lwc)
     fw_h2o2 = MIN(1.0, MAX(0.0, fw_h2o2))

   END SUBROUTINE dissolved_gas_frac

!----------------------------------------------------------------------------

   SUBROUTINE aer_aqchem_simple_init(chem_opt)

     ! Purpose
     !  Initialize aer_wetscav_simple parameters, which depend on chem_opt
     !
     ! Arguments
     !  chem_opt TODO
     !
     ! Remarks
     !  TODO
     !  TODO check that species exist, define naerbins (1 for GOCART)

     USE module_state_description ! For chem_opt option values

     IMPLICIT NONE

     ! Arguments
     !  Input
     INTEGER, INTENT(IN) :: chem_opt

    !---- Initialize chem mechanism families IDs
    aer_is_gocart = chem_opt == MOZCART_KPP .OR. &
                    chem_opt == T1_MOZCART_KPP .OR. &
                    chem_opt == GOCART_SIMPLE .OR. &
                    chem_opt == GOCARTRADM2 .OR. &
                    chem_opt == GOCARTRACM_KPP

    aer_is_mosaic4bin = chem_opt == CBMZ_MOSAIC_4BIN .OR. &
                        chem_opt == CBMZ_MOSAIC_DMS_4BIN .OR. &
                        chem_opt == CBMZ_MOSAIC_4BIN_AQ .OR. &
                        chem_opt == CBMZ_MOSAIC_DMS_4BIN_AQ .OR. &
                        chem_opt == MOZART_MOSAIC_4BIN_KPP .OR. &
                        chem_opt == MOZART_MOSAIC_4BIN_AQ_KPP .OR. &
                        chem_opt == MOZART_MOSAIC_4BIN_DMS_AQ_KPP .OR. &
                        chem_opt == CRI_MOSAIC_4BIN_AQ_KPP .OR. &
                        chem_opt == SAPRC99_MOSAIC_4BIN_VBS2_KPP

    aer_is_mosaic8bin = chem_opt == CBMZ_MOSAIC_8BIN .OR. &
                        chem_opt == CBMZ_MOSAIC_DMS_8BIN .OR. &
                        chem_opt == CBMZ_MOSAIC_8BIN_AQ .OR. &
                        chem_opt == CBMZ_MOSAIC_DMS_8BIN_AQ .OR. &
                        chem_opt == CRI_MOSAIC_8BIN_AQ_KPP .OR. &
                        chem_opt == SAPRC99_MOSAIC_8BIN_VBS2_AQ_KPP .OR. &
                        chem_opt == SAPRC99_MOSAIC_8BIN_VBS2_KPP


    !-------- Initialize aerosol bin/mode number for chem_opt --------
    ! This could be recovered from the aerosol scheme data directly
    IF(aer_is_gocart) THEN
      ! GOCART cases with 1 bulk sulfate value
      naersizes = 1
    ELSEIF(aer_is_mosaic8bin) THEN
      ! MOSAIC-8bin cases
      naersizes = 8
    ELSEIF(aer_is_mosaic4bin) THEN
      ! MOSAIC-4bin cases
      naersizes = 4
    ELSE
      CALL wrf_error_fatal ('aer_aqchem_simple_init: chem_opt not supported')
    ENDIF

   END SUBROUTINE aer_aqchem_simple_init

!----------------------------------------------------------------------------

END MODULE module_aer_aqchem_simple

