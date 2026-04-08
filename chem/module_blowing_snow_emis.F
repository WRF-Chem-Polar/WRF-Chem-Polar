! Louis Marelle, Jennie Thomas, Shaddy Ahmed, 2022/09/13
!
! Purpose:
!  Compute size-resolved blowing snow aerosol emissions over sea ice (mass and number) for
!  the sectional 4-bin or 8-bin MOSAIC aerosols, or for the GOCART 4 sea salt
!  bins.
!
! Remarks:
!  - TODO this should be coupled to the dynamics, emiting QSNOW and QVAPOR and changing
!    snowh, snowc, snowdens. Blowing snow should also depend on snow density
!    and blowing snow events should increase density.
!  - TODO meltedsnow and snowh_melt should (1) be calculated in the physics
!    code when chem is used or at least in an external subroutine in chem (2)
!    drift with sea ice (3) be reset in the fall, e.g. when sea-ice reforms in
!    october and later or with new snow fall
!  - TODO This should ideally read in observed sea ice state (snow state, snow
!    salinity etc.), or include some more detailed snow salinity profiles in
!    snow

! References:
! =============================================================================
! (1) Marelle, Thomas, Ahmed et al., JAMES, 2021:
!     https://doi.org/10.1029/2020MS002391
! (2) Yang et al., 2008: doi:10.1029/2008GL034536.
! (3) Déry & Yau, 2001: https://doi.org/10.1023/A:1018965008049
! (4) Huang and Jaeglé, 2017: https://doi.org/10.5194/acp-17-3699-2017
!

MODULE module_blowing_snow_emis

  IMPLICIT NONE

  !---- Global module-wide variables
  ! These are initialized in blowing_snow_emis_init

  ! Allocatable size-resolved aerosol arrays (dim is size bin number and depends
  ! on aerosol mechanism, emiss_massfrac_nacl is unitless (0-1),
  ! bin_diameters_cen in um)
  !
  ! Mass proportion of total NaCl emissions going into each aerosol size bin,
  ! calculated offline from the size+salinity integral in Yang et al. (2008),
  ! eq. 8
  ! emiss_massfrac_nacl depends on the snow salinity and aerosol chem mechanism
  ! (bin limits). If the snow_salinity parameter is changed this needs to be
  ! recomputed offline (see Yang et al., 2008). Values are based on a fixed snow
  ! snow salinity (see snow_salinity) N=5 aerosol fractioning, Huang and Jaeglé
  ! (2017)
  REAL, DIMENSION(:), ALLOCATABLE :: emiss_massfrac_nacl
  !
  ! Aerosol bin center diameters in micro m
  REAL, DIMENSION(:), ALLOCATABLE :: bin_diameters_cen
  !
  ! Average snow salinity (psu). Do not update snow_salinity without updating
  ! emiss_massfrac_nacl in init routine
  ! Huang and Jaeglé, 2017, Arctic, rang 0.01 to 1.7 psu
  REAL :: snow_salinity = 0.1
  !
  ! Number of sea salt aerosol bins for chem_opt in use
  INTEGER :: naersizes

  ! chem_opt aerosol model family
  LOGICAL :: aer_is_gocart
  LOGICAL :: aer_is_mosaic4bin
  LOGICAL :: aer_is_mosaic8bin

  CONTAINS

!------------------------------------------------------------------------------

  SUBROUTINE blowing_snow_emis( id, dtstep, dz8w, moist, t_phy, p_phy, &
         rho_phy, chem, u10, v10, ust, T2,            &
         snowh, tsk, xice, lakemask, meltedsnow, snowh_melt, &
         e_nacl_blow,       &
         ids,ide, jds,jde, kds,kde,                                 &
         ims,ime, jms,jme, kms,kme,                                 &
         its,ite, jts,jte, kts,kte                                  )

    ! Purpose:
    !  Compute size-resolved blowing snow aerosol number and mass emissions over
    !  sea ice

    ! Arguments:
    ! Input:
    !  id = domain ID number
    !  dtstep = model time step (s)
    !  dz8w = 3D vertical level depth (m)
    !  moist = advected water species (kg/kg)
    !  t_phy = 3D air temperature (K)
    !  p_phy = 3D air pressure (Pa)
    !  rho_phy = 3D air density (kg/m3)
    !  u10 = windspeed at 10m in x (positive west-east) direction (m/s)
    !  v10 = windspeed at 10m in y (positive south-north) direction (m/s)
    !  ust = friction velocity (m/s)
    !  T2 = 2-m temperature (K)
    !  snowh = snow height (m)
    !  tsk = surface skin temperature (K)
    !  xice = sea ice or lake ice fraction (0-1)
    !  lakemask = lake mask (1 for lakes, 0 for rest)
    !  meltedsnow = past snow melt mask (1 for past melt, 0 for no melt)
    !  snowh_melt = snow height at the end of the last melt period (m)
    !  ids,ims,its etc. = domain, memory, and tile start and end indices
    ! Input/output
    !  chem = advected chemical species (gases in ppm, aerosols in microg/kg)
    !  e_nacl_blow = accumulated NaCl emission flux from blowing snow (microg/m2) !TODO should be mg/m2

    !TODO include all explicit ONLY: associations for module_state_description
    USE module_configure
    USE module_state_description, ONLY: param_first_scalar
    USE module_data_mosaic_asect, ONLY: numptr_aer, lptr_cl_aer, lptr_na_aer
    USE module_data_mosaic_constants, ONLY: pi
    USE module_model_constants, ONLY: mwdry, ep_2, r_v ! physical constants - molar weight of dry air etc.
    USE module_mp_morr_two_moment, ONLY: POLYSVP

    IMPLICIT NONE

    !---- Subroutine arguments
    ! Input
    INTEGER, INTENT(IN) :: id,                                      &
                           ids,ide, jds,jde, kds,kde,               &
                           ims,ime, jms,jme, kms,kme,               &
                           its,ite, jts,jte, kts,kte
    REAL, INTENT(IN) :: dtstep
    REAL, DIMENSION(ims:ime, kms:kme, jms:jme),                     &
          INTENT(IN) :: dz8w, t_phy, p_phy, rho_phy
    REAL, DIMENSION( ims:ime, kms:kme, jms:jme, num_moist ),        &
          INTENT(IN ) :: moist
    REAL, DIMENSION(ims:ime, jms:jme)               ,               &
          INTENT(IN) :: u10, v10, ust, t2, snowh, tsk, xice, lakemask
    ! Input/output
    REAL, DIMENSION(ims:ime, jms:jme)               ,               &
          INTENT(INOUT) :: meltedsnow, snowh_melt
    REAL, DIMENSION(ims:ime, jms:jme)               ,               &
          INTENT(INOUT) :: e_nacl_blow
    REAL, DIMENSION(ims:ime, kms:kme, jms:jme, num_chem),           &
          INTENT(INOUT) :: chem

    !---- Parameters
    ! Temperature threshold for emissions (K), sea ice cover
    ! threshold for blowing snow emissions (area fraction)
    REAL, PARAMETER :: tskcrit = 273.15, seaice_crit = 0.75
    ! Conversion factor from kg/m2/day to kg/m2/s
    REAL, PARAMETER :: kconv_day_s = 1.1574E-5
    ! Snow age (hours), value from Huang and Jaeglé (2017) - this value might be
    ! better than using actual snow age
    REAL, PARAMETER :: snow_age = 72.0
    ! Reference 10m-wind speed (m/s)
    REAL, PARAMETER :: w10_0 = 6.975
    ! Parameters for qsnow_subl_prime parameterization, from Dery and Yau (2001)
    REAL, PARAMETER :: a0 = 3.78407E-1, a1 = -8.64089E-2, a2 = -1.60570E-2, &
            a3 = 7.25516E-4, a4 = -1.25650E-1, a5 = 2.48430E-2, &
            a6 = -9.56871E-4, a7 = 1.24600E-2, a8 = 1.56862E-3, &
            a9 = -2.93002E-4
    ! Conversion factor for the ksi parameter (unitless)
    REAL, PARAMETER :: convert_ksi = -1.0E12
    ! Density of ice (kg/m3) and sea salt (kg/m3)
    REAL, PARAMETER :: rho_ice = 914.0, rho_nacl = 2160.
    REAL, PARAMETER :: NACL_AERDEN = 2.17 * 1.0E6 * 1.0E-12 ! g/cm3 -> ug/cm3 -> ug/um3


    !---- Local variables
    CHARACTER (LEN=1000) :: msg
    ! Loop indices i j k
    INTEGER :: i, j, k
    ! Conversion factor from emissions to µg/kg, Conversion factor from number
    ! emissions to #/kg
    REAL :: conv2, conv3
    ! Conversion factor from 1/m2/s to 1/kg_dry_air
    REAL :: conv
    ! aerosol bin center diameter (um)
    REAL :: dcen
    ! NaCl emissions (kg/m2/s), Snow sublimation flux (kg/m2/s), Normalized
    ! column-integrated snow sublimation flux (kg/m2/day), Snow mixing ratio in
    ! saltation layer (kg/kg), Snow mixing ratio at reference wind speed (kg/kg)
    REAL :: qsnow_subl, qsnow_subl_prime, qsnow_saltation, qsnow_0
    ! 10m-windspeed (m/s), critical 10m-windspeed (m/s), temperature-moisture
    ! parameter from Dery and Yau 2001 (-1.0E-12 m2/s-1), ice saturation mixing
    ! ratio (kg/kg), conductivity term in the sublimation process (m s kg−1),
    ! diffusion term in the sublimation process (m s kg−1), saturation water vapor
    ! pressure wrt ice (Pa)
    REAL :: w10, w10_crit, ksi, qvapor_sat_ice, fconduction, fdiffusion, eis
    ! Molecular diffusion coefficient of water vapor (m2 s-1), latent heat of
    ! sublimation (J kg-1), coefficient of thermal conductivity of dry air (J m-1
    ! s-1 K-1),
    REAL :: d_v, l_s, k_d
    ! WRF-Chem aerosol aerosol type, size, phase indexes, and
    ! "pointers" to number, Na and Cl in chemical composition array chem
    INTEGER :: itype=1, iaerbin, iphase=1, l_num, l_na, l_cl, l_seasalt
    ! Adjustment factor for the emission flux (unitless)
    REAL :: a_prime
    ! Total sea salt emission flux (ug/m2/s)
    REAL :: nacl_emiss_tot
    ! Sea salt emission fluxes um3/m2/s, ug/m2/s, #/m2/s in bin
    REAL :: nacl_vol_emiss, nacl_mass_emiss, nacl_num_emiss
    ! Na and Cl molar weights in g/mol
    REAL, PARAMETER :: MW_NA=22.990, MW_CL=35.453


    !-------- Calculate emissions from blowing snow --------
    DO i = its,ite
      DO j = jts,jte
        ! Conversion factor from 1/m2/s to 1/kg_dry_air
        conv = 1.0/(rho_phy(i,kts,j) * dz8w(i,kts,j)) * dtstep

        !---- Blowing snow: calculate the total NaCl emission flux from blowing snow
        nacl_emiss_tot = 0.0
        ! Only valid over sea-ice, excluding lake ice
        IF (xice(i, j) .GT. seaice_crit .AND. lakemask(i,j) .LT. 0.5) THEN
          !SAhmed - if snow melts, turn off surface emissions there until
          ! temperatures are back below freezing AND fresh snow has been
          ! deposited (at least 1 cm).
          IF(tsk(i, j) .GT. tskcrit) THEN
            snowh_melt(i,j) = snowh(i,j)
            meltedsnow(i,j) = 1.
          ELSEIF (snowh(i,j) .GT. snowh_melt(i,j)+0.01) THEN
            snowh_melt(i,j) = 0.0
            meltedsnow(i,j) = 0.0
          END IF
          ! Only valid below freezing temperatures, and over snow that has not
          ! experienced past melt
          IF (tsk(i, j) .LE. tskcrit .AND. meltedsnow(i,j) .LT. 0.5) THEN
            !-- Calculate the critical 10m wind speed for blowing snow, from Yang
            ! et al. (2008)
            w10_crit = w10_0 + 0.0033 * (t2(i,j) - 273.15 + 27.27) ** 2.0
            ! Wind speed in m s-1
            w10 = sqrt( u10(i,j) ** 2.0 + v10(i,j) ** 2.0 )
            ! Only valid above critical wind speed
            IF (w10 .gt. w10_crit) THEN
              !-- Calculate the snow sublimation flux, qsnow_subl, based on Yang et
              ! al. (2008) and Dery and Yau (2001), with help from Markus Frey
              ! for the fconduction and fdiffusion terms
              l_s = 2.83458E6 - (t2(i,j) - 273.15) * (340. + 10.46 * (t2(i,j) - 273.15))
              k_d = 0.023807 + 7.1128E-5 * (t2(i,j) - 273.15)
              d_v = 2.11E-5 * (t2(i,j) / 273.15) ** 1.94 * (1013.25 / (p_phy(i, kts, j) * 1.0E-2))
              ! qvapor_sat_ice calculated from the morrison microphysics module
              ! POLYSVP function
              eis = min(0.99 * p_phy(i,kts,j), POLYSVP(t_phy(i,kts,j), 1))
              qvapor_sat_ice = ep_2 * eis / (p_phy(i,kts,j) - eis)
              fconduction = (l_s / (r_v * t2(i,j)) - 1.0) * l_s / (k_d * t2(i,j))
              fdiffusion = (r_v * t2(i,j)) / (d_v * eis)
              ! From Dery and Yau (2001) eq. 11
              ksi = convert_ksi * (moist(i, kts, j, p_qv) / qvapor_sat_ice - 1.0) &
                    / (2.0 * rho_ice * (fconduction + fdiffusion))
              ! From Dery and Yau (2001) eq. 14
              qsnow_subl_prime = a0 + a1 * ksi + a2 * ksi ** 2.0 + a3 * ksi ** 3.0 &
                + a4 * w10 + a5 * ksi * w10 + a6 * ksi ** 2.0 * w10 + a7 * w10 ** 2.0 &
                + a8 * ksi * w10 ** 2 + a9 * w10 ** 3.0
              ! From Yang et al., (2008)
              a_prime = 1.0 / (1.038 + 0.03758 * snow_age - 0.00014349 * snow_age ** 2.0 &
                        + 1.911315E-7 * snow_age ** 3.0)
              qsnow_saltation = 0.385 * (1. - w10_crit / w10) ** 2.59 / ust(i,j)
              qsnow_0 = 0.385 * (1. - w10_0 / w10) ** 2.59 / ust(i,j)
              qsnow_subl = kconv_day_s * a_prime * qsnow_subl_prime * qsnow_saltation / qsnow_0
              !-- Calculate the total sea salt emission flux (ug/m2/s), from Yang et al.,
              ! (2008), scaled by the ice fraction.
              nacl_emiss_tot = 1.0E9 * qsnow_subl * snow_salinity / 1000. * xice(i, j)
              IF (nacl_emiss_tot < 0.) THEN
                ! TODO understand why this happens sometimes
                WRITE(msg,*) 'module_blowing_snow_emis: nacl_emiss_tot<0, ',nacl_emiss_tot,' at i,j=',i,j
                CALL wrf_debug(100, TRIM(msg))
              ENDIF
              nacl_emiss_tot = MAX(nacl_emiss_tot, 0.)

            END IF ! w10 .gt. w10_crit
          END IF ! tsk <= tskcrit
        END IF ! xice(i, j) .GT. seaice_crit .AND. meltedsnow


        !---- Blowing snow: distribute sea salt emissions to aer. size bins
        IF(nacl_emiss_tot > 0.) THEN
          DO iaerbin = 1, naersizes
            dcen = bin_diameters_cen(iaerbin)
            nacl_mass_emiss = nacl_emiss_tot * emiss_massfrac_nacl(iaerbin)
            nacl_vol_emiss = nacl_mass_emiss / NACL_AERDEN
            nacl_num_emiss = nacl_vol_emiss / ((PI/6.0)*dcen**3.0)
            !---- Accumulate nacl emissions from blowing snow, output diagnosis
            e_nacl_blow(i, j) = e_nacl_blow(i, j) + nacl_mass_emiss * dtstep

            !---- Emit blowing snow emissions (add them to chem) --------
            !---- Emit MOSAIC aerosols to chem
            IF(aer_is_mosaic4bin .OR. aer_is_mosaic8bin) THEN
              ! MOSAIC sodium aerosol
              l_na = lptr_na_aer(iaerbin,itype,iphase)
              IF (l_na >= param_first_scalar) THEN
                chem(i,kts,j,l_na) = chem(i,kts,j,l_na) + nacl_mass_emiss &
                                     * MW_NA/(MW_NA+MW_CL) * conv
              END IF
              ! MOSAIC chlorine aerosol
              l_cl = lptr_cl_aer(iaerbin,itype,iphase)
              IF (l_cl >= param_first_scalar) THEN
                chem(i,kts,j,l_cl) = chem(i,kts,j,l_cl) + nacl_mass_emiss &
                                     * MW_CL/(MW_NA+MW_CL) * conv
              END IF
              ! MOSAIC aerosol number
              l_num = numptr_aer(iaerbin,itype,iphase)
              IF (l_num >= param_first_scalar) THEN
                chem(i,kts,j,l_num) = chem(i,kts,j,l_num) + nacl_num_emiss * conv
              END IF
            !---- Emit GOCART aerosols to chem
            ELSEIF(aer_is_gocart) THEN
              ! GOCART seasalt aerosol
              IF(iaerbin == 1 ) THEN
                l_seasalt = p_seas_1
              ELSEIF(iaerbin == 2 ) THEN
                l_seasalt = p_seas_2
              ELSEIF(iaerbin == 3 ) THEN
                l_seasalt = p_seas_3
              ELSEIF(iaerbin == 4 ) THEN
                l_seasalt = p_seas_4
              ENDIF
              IF (l_seasalt >= param_first_scalar) THEN
                chem(i,kts,j,l_seasalt) = chem(i,kts,j,l_seasalt) + nacl_mass_emiss * conv
              ENDIF
            !---- Other chem_opt cases are not supported
            ELSE
              CALL wrf_error_fatal ('blowing_snow_emis: chem_opt not supported')
            ENDIF ! aer_is_
          ENDDO ! iaerbin
        ENDIF ! nacl_emiss_tot > 0.

      ENDDO ! j
    ENDDO ! i

  END SUBROUTINE blowing_snow_emis

!------------------------------------------------------------------------------

  SUBROUTINE blowing_snow_emis_init(chem_opt)

    ! Purpose:
    !  Initialize the parameters for the blowing snow aerosol emission scheme,
    !  which depend on chem_opt.
    !
    ! Remarks:
    ! - emiss_massfrac_nacl was precalculated by integrating the source functions
    !   size distributions from the original publications over the WRF-Chem
    !   aerosol size ranges. For this reason if the aerosol bin limits are changed
    !   emiss_massfrac_nacl needs to be recalculated.
    ! - TODO Center diameters and naersizes should be compared to the ones used at
    !   runtime (for MOSAIC, dcen_sect(iaerbin,itype) and nsize_aer(itype) from
    !   module_data_mosaic_asect), and model should crash if there is a mismatch
    !   - TODO If comparing to MOSAIC nsizes and dcen_sect, careful the MOSAIC top
    !     level initializations should have been done before.
    ! - TODO Allocatables here are never deallocated
    ! - TODO some of this could be in a separate data module

    !TODO only use with only:, specify all variables here
    USE module_state_description ! For chem_opt option values

    IMPLICIT NONE

    ! Dummy arguments:
    ! Input
    INTEGER, INTENT(IN) :: chem_opt

    !---- Initialize chem mechanism families IDs
    !TODO Could be automated by checking the presence of mechanism-specific pointers
    !TODO this needs to be updated if a new mechanism is defined
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

    !-------- Initialize sea salt aerosol bin number for chem_opt --------
    ! This could be recovered from the aerosol scheme data directly, but
    ! emiss_massfrac_nacl is already chem_opt-dependent and depends on naersizes so
    ! there is no harm in redefining it here just for this module
    IF(aer_is_gocart) THEN
      ! GOCART with 4 sea salt bins cases
      naersizes = 4
    ELSEIF(aer_is_mosaic8bin) THEN
      ! MOSAIC-8bin cases
      naersizes = 8
    ELSEIF(aer_is_mosaic4bin) THEN
      ! MOSAIC-4bin cases
      naersizes = 4
    ELSE
      CALL wrf_error_fatal ('blowing_snow_emis_init: chem_opt not supported')
    ENDIF

    !-------- Initialize size-resolved aerosol arrays for SSA emissions --------
    ! Initialize emiss_massfrac_nacl (0-1 fraction) and bin_diameters_cen
    ! (micrometer) for each family of aerosol options selected by chem_opt

    ! Allocate the size-resolved aerosol arrays
    IF (.NOT. ALLOCATED(emiss_massfrac_nacl)) THEN
      ALLOCATE(emiss_massfrac_nacl(naersizes))
    ENDIF
    IF (.NOT. ALLOCATED(bin_diameters_cen)) THEN
      ALLOCATE(bin_diameters_cen(naersizes))
    ENDIF

    !-- Initialize NaCl size bin massfrac for GOCART aerosols
    IF(aer_is_gocart) THEN
      ! GOCART sea salt bin center diameters (micrometer)
      ! These should be log-averaged but I follow the GOCART definition instead
      bin_diameters_cen(:) = (/ 0.60, 2.00, 6.5, 15.0 /)
      ! Size-resolved NaCl mass repartition in GOCART size bins
      emiss_massfrac_nacl(:) = (/ 0.8737, 0.1261, 1.491E-04, 0. /)

    !-- Initialize NaCl size bin massfrac for MOSAIC-8bin aerosols
    ELSEIF(aer_is_mosaic8bin) THEN
      ! MOSAIC-8bin bin center diameters (micrometer)
      bin_diameters_cen(:) = (/ 0.0552, 0.1105, 0.2210, 0.4419, &
                             0.8839, 1.7678, 3.5355, 7.0711 /)
      ! Size-resolved NaCl mass repartition in MOSAIC-8bin size bins
      emiss_massfrac_nacl(:) = (/ 2.711E-02, 8.600E-02, 2.181E-01, 0.3580, &
                               2.629E-01, 4.721E-02, 7.154E-04, 0. /)

    !-- Initialize NaCl size bin massfrac for MOSAIC-4bin aerosols
    ELSEIF(aer_is_mosaic4bin) THEN
      ! MOSAIC-4bin bin center diameters (micrometer)
      bin_diameters_cen(:) = (/ 0.0781, 0.3125, 1.2500, 5.0000 /)
      ! Size-resolved NaCl mass repartition in MOSAIC-4bin size bins
      emiss_massfrac_nacl(:) = (/ 0.1131, 0.5761, 0.3101, 0.7155E-03 /)

    ELSE ! aer_is_
      CALL wrf_error_fatal ('blowing_snow_emis_init: chem_opt not supported')

    ENDIF ! aer_is_

  END SUBROUTINE blowing_snow_emis_init

END MODULE module_blowing_snow_emis

