! 2025/06 Louis Marelle
!
! Purpose:
!  Calculate dry deposition and sedimentation velocities for aerosols at the
!  surface, following Emerson et al. (2020), a modified version of the Zhang
!  et al. (2001) parameterization with updated parameters.
!
! Remarks:
!  - Only compatible with GOCART for now
!  - Sedimentation above the first level is already done in gocart_settling or
!    mosaic_settling
!  - Hygroscopic growth is only taken into account for seasalt and sulfate
!
! References:
! ------------------------------------------------------------------------------
!  (1) Emerson et al. (2020) -  Revisiting particle dry deposition and its
!      role in radiative effect estimates, Proc. Natl. Acad. Sci. U.S.A.  117
!      (42) 26076-26082, https://doi.org/10.1073/pnas.2014761117
!  (2) Zhang et al., (2001) - A size-segregated particle dry deposition scheme
!      for an atmospheric aerosol module. Atmospheric environment, 35(3),
!      549-560, https://doi.org/10.1016/S1352-2310(00)00326-5
!  (3) Gerber (1985) - Relative humidity parameterization of the Navy aerosol
!      model (NAM), NRL Rep. 8956. National Research Laboratory, Washington
!      DC. https://apps.dtic.mil/sti/tr/pdf/ADA163209.pdf

MODULE module_aer_drydep_simple

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: aer_drydep_simple, aer_drydep_simple_init

  INTEGER, PARAMETER :: nland_cats_zhang = 15
  INTEGER, PARAMETER :: nseasons = 5

  ! Number of wet-scavenged aerosol species for chem_opt
  INTEGER :: naerspecies
  ! List of wet-scavenged aerosol species for chem_opt
  INTEGER, DIMENSION(:), ALLOCATABLE :: aer_pointers
  ! Aerosol dry diameter for aerosol species in aer_pointers (m)
  REAL, DIMENSION(:), ALLOCATABLE :: aer_dry_diam
  ! Aerosol densities for aerosol species in aer_pointers (kg/m3)
  REAL, DIMENSION(:), ALLOCATABLE :: aer_dry_dens
  ! Land cover type name 
  CHARACTER(LEN=256) :: mminlu_loc
  ! Radius (m) of surface collectors for interception
  ! Interception efficiency land-dependent coefficient, from Zhang2001
  REAL, SAVE :: RAD_COLLECTOR_ZHANG(5,15)
  DATA RAD_COLLECTOR_ZHANG/ &
         0.002, 0.002, 0.002, 0.002, 0.002, &  ! luc 1
         0.005, 0.005, 0.005, 0.005, 0.005, &  ! luc 2
         0.002, 0.002, 0.005, 0.005, 0.002, &  ! luc 3
         0.005, 0.005, 0.010, 0.010, 0.005, &  ! luc 4
         0.005, 0.005, 0.005, 0.005, 0.005, &  ! luc 5
         0.002, 0.002, 0.005, 0.005, 0.002, &  ! luc 6
         0.002, 0.002, 0.005, 0.005, 0.002, &  ! luc 7
         -9.99, -9.99, -9.99, -9.99, -9.99, &  ! luc 8
         -9.99, -9.99, -9.99, -9.99, -9.99, &  ! luc 9
         0.010, 0.010, 0.010, 0.010, 0.010, &  ! luc 10
         0.010, 0.010, 0.010, 0.010, 0.010, &  ! luc 11
         -9.99, -9.99, -9.99, -9.99, -9.99, &  ! luc 12
         -9.99, -9.99, -9.99, -9.99, -9.99, &  ! luc 13
         -9.99, -9.99, -9.99, -9.99, -9.99, &  ! luc 14
         0.010, 0.010, 0.010, 0.010, 0.010  /  ! luc 15

  ! Maps MODIS surface categories to zhang categories
  INTEGER, SAVE  :: ZHANG_SURF_FROM_MODIS(21)
  DATA ZHANG_SURF_FROM_MODIS / &
    !   ! zhang category        <---- modis category
    1, &! evergreen needle-leaf    1 'Evergreen Needleleaf Forest'
    2, &! evergreen broadleaf      2 'Evergreen Broadleaf Forest'
    3, &! deciduous needle-leaf    3 'Deciduous Needleleaf Forest'
    4, &! deciduous broadleaf      4 'Deciduous Broadleaf Forest'
    5, &! mixed broad/needle-leaf  5 'Mixed Forests'
   10, &! shrubs & woodland        6 'Closed Shrublands'
   10, &! shrubs & woodland        7 'Open Shrublands'
   10, &! shrubs & woodland        8 'Woody Savannas'
   10, &! shrubs & woodland        9 'Savannas'
    6, &! grass                   10 'Grasslands'
   11, &! wetland w plants        11 'Permanent wetlands'
    7, &! crops, mixed farming    12 'Croplands'
   15, &! urban                   13 'Urban and Built-Up'
    7, &! crops, mixed farming    14 'cropland/natural vegetation mosaic'
   12, &! ice cap & glacier       15 'Snow and Ice'
    8, &! desert                  16 'Barren or Sparsely Vegetated'
   14, &! ocean                   17 'Water'
   10, &! shrubs & woodland       18 'Wooded Tundra'
    9, &! tundra                  19 'Mixed Tundra'
    8, &! desert                  20 'Barren Tundra'
    13/ ! inland water            21 'Lakes'

  ! Maps USGS surface categories to zhang categories
  INTEGER, SAVE :: ZHANGSURF_FROM_USGS(25)
  DATA ZHANGSURF_FROM_USGS / &
    !     ! zhang category        <---  usgs category
   15, &  ! urban                    1: urban and built-up land
    7, &  ! crops, mixed farming     2: dryland cropland & pasture
    7, &  !                          3: irrigated crop & pasture
    7, &  !                          4: mix dry/irrig crop & pasture
    7, &  !                          5: cropland/grassland mosaic
   10, &  ! shrubs & woodland        6: cropland/woodland mosaic
    6, &  ! grass                    7: grassland
    6, &  !                          8: shrubland
    6, &  !                          9: mixed shrubland/grassland
    6, &  !                         10: savanna
    4, &  ! deciduous broadleaf     11: deciduous broadleaf forest
    3, &  ! deciduous needle-lf     12: deciduous needleleaf forest
    2, &  ! evergreen broadleaf     13: evergreen broadleaf forest
    1, &  ! evergreen needle-lf     14: evergreen needleleaf forest
    5, &  ! mixed broad/needle-lf   15: mixed forest
   14, &  ! ocean                   16: water bodies
   11, &  ! wetland w plants        17: herbaceous wetland
   10, &  !                         18: wooded wetland
    8, &  ! desert                  19: barren or sparsely vegetated
    9, &  ! tundra                  20: herbaceous tundra
   10, &  !                         21: wooded tundra
    9, &  !                         22: mixed tundra
    8, &  !                         23: bare ground tundra
   12, &  ! ice cap & glacier       24: snow or ice
    8/    !                         25: no data

  ! chem_opt aerosol model type
  LOGICAL :: aer_is_gocart

CONTAINS

!----------------------------------------------------------------------------

  SUBROUTINE aer_drydep_simple(ddvel, t_phy, rho_phy, ust, rh, &
                               xlat, snowh, xice, ivgtyp,      &
                               res_adyn, julday, dtstep,       &
                               ids,ide, jds,jde, kds,kde,      &
                               ims,ime, jms,jme, kms,kme,      &
                               its,ite, jts,jte, kts,kte       )
   ! Purpose:
   !  Calculates dry deposition velocities for aerosols
   !  
   ! Arguments:
   !  Input:
   !    t_phy = air temperature (degree celsius)
   !    rho_phy = air density (kg/m3)
   !    ust = friction velocity (m/s)
   !    rh = relative humidity (0-1)
   !    xlat = latitude (-90 to 90)
   !    snowh = snow depth (m)
   !    xice = sea ice concentration (fraction)
   !    ivgtyp = land use category (categorical)
   !    res_adyn = aerodynamical deposition resistance from wesely (s/m)
   !    julday = julian day (1-366)
   !    dtstep = model time step (s)
   !    ids,ims,its etc. = domain, memory, and tile start and end indices
   !  Input/output:
   !    ddvel = dry deposition velocities (m/s)

   ! USE associations:
   USE module_model_constants, ONLY: g,piconst
   USE module_state_description, ONLY: p_sulf, p_seas_1, p_seas_2, p_seas_3, &
     p_seas_4, p_dust_1, p_dust_2, p_dust_3, p_dust_4, p_dust_5, num_chem, &
     param_first_scalar

   IMPLICIT NONE

   !---- Arguments
   REAL, DIMENSION(ims:ime, kms:kme, jms:jme), INTENT(IN) :: t_phy, rho_phy, &
                                                             rh
   REAL, DIMENSION(ims:ime, jms:jme), INTENT(IN) :: ust, xlat, snowh, xice
   INTEGER, DIMENSION(ims:ime, jms:jme), INTENT(IN) :: ivgtyp
   REAL, DIMENSION(its:ite, jts:jte), INTENT(IN) :: res_adyn
   REAL, DIMENSION(its:ite, jts:jte, num_chem),          &
       INTENT(INOUT) :: ddvel
   REAL, INTENT(IN) :: dtstep
   INTEGER, INTENT(IN) :: julday
   INTEGER, INTENT(IN) :: ids,ide, jds,jde, kds,kde,    &
                          ims,ime, jms,jme, kms,kme,    &
                          its,ite, jts,jte, kts,kte

   !---- Local variables
   ! Dry deposition velocity (m/s)
   REAL :: vdep
   ! Settling velocity in lowest level (m/s)
   REAL :: vsettl
   ! Surface deposition resistance (s/m)
   REAL :: res_sfc
   ! Air temperature (K), air density (g/cm3), relative humidity (0-1)
   REAL :: tk, rho_gcm3, relhum
   ! Friction velocity (cm/s)
   REAL :: ust_cms
   ! Mean free path of air molecules (cm)
   REAL :: freepath
   ! Aerodynamical deposition resistance (s/cm)
   REAL :: res_adyn_scm
   ! Aerosol dry and wet densities (g/cm3)
   REAL :: aer_dry_dens_gcm3, aer_wet_dens_gcm3
   ! Stokes and schmidt numbers (dimensionless)
   REAL :: stokes, schmidt
   ! Particle Brownian diffusivity (cm2/s)
   REAL :: aerodiffus
   ! Aerosol dry and wet diameters (cm)
   REAL :: aer_dry_d_cm, aer_wet_d_cm
   ! Air kinematic viscosity (cm2/s)
   REAL :: airkinvisc
   ! Cunningham correction factor
   REAL :: cunningham_fact
   ! Interception collection radius parameter (cm)
   REAL :: aintercept
   ! Deposition efficiencies from brownian motion, impaction, and interception
   REAL :: ebrown, eimpact, eintercept
   ! Particle bounce correction factor
   REAL :: r1fac
   ! x,y, and chemical species index, chemical species pointer
   INTEGER :: i, j, iaer, p_chem
   ! Zhang2001 surface category
   INTEGER :: zhangsurfcat
   ! WRF (USGS or MODIS) surface category
   INTEGER :: wrfsurfcat
   ! Zhang2001 Seasonal category
   INTEGER :: iseason
   ! Flag for vegetated surfaces (true if vegetation is present)
   LOGICAL :: sfc_is_vegetated
   ! Flag for frozen surfaces (true if ice or snow is present)
   LOGICAL :: sfc_is_frozen

   !---- Parameters
   ! Surface resistance parameterization empirical coefficient
   REAL, PARAMETER :: EPSZ = 3.
   ! Parameters for the deposition efficiencies, from Emerson2020
   REAL, PARAMETER :: CBROWN = 0.2
   REAL, PARAMETER :: CIMPACT = 0.4, ALPHAIMPACT = 0.8, BETAIMPACT = 1.7
   REAL, PARAMETER :: CINTERCEPT = 2.5, NUINTERCEPT = 0.8
   REAL, PARAMETER :: BOLTZMANN = 1.3807E-16 ! erg/K = g cm2 s-2 K-1
   REAL, PARAMETER :: PI = piconst
   REAL, PARAMETER :: GRAVITY = g*100.0 ! cm s-2
   ! Parameters for the wet diameter calculation of Gerber (1985)
   REAL, PARAMETER :: C1_WET_SEASALT = 0.7674
   REAL, PARAMETER :: C2_WET_SEASALT = 3.079
   REAL, PARAMETER :: C3_WET_SEASALT = 2.572E-11
   REAL, PARAMETER :: C4_WET_SEASALT = -1.424
   REAL, PARAMETER :: C1_WET_SULFATE = 0.4809
   REAL, PARAMETER :: C2_WET_SULFATE = 3.082
   REAL, PARAMETER :: C3_WET_SULFATE = 3.110E-11
   REAL, PARAMETER :: C4_WET_SULFATE = -1.428

   !-------- Calculate dry deposition velocities for aerosol species -------
   DO j=jts,jte
     DO i=its,ite
       !---- Initialize physical parameters and convert units
       tk = t_phy(i,kts,j) + 273.15 ! K
       rho_gcm3 = rho_phy(i,kts,j) * 1E-3 ! g/cm3
       ust_cms = ust(i,j) * 100.0 ! cm/s
       relhum = rh(i,kts,j) ! 0-1
       ! Aerodynamic resistance
       res_adyn_scm = res_adyn(i,j) * 0.01 ! s/cm
       ! Air kinematic viscosity (cm2/s)
       airkinvisc = (1.83E-4 * (416.16/(tk+120.0)) * ((tk/296.16)**1.5)) &
                    / rho_gcm3
       ! Air molecular freepath (cm)
       freepath = 7.39758E-4 * airkinvisc / SQRT(tk)
       ! Land cover category from WRF
       wrfsurfcat = ivgtyp(i,j)
       ! is the surface frozen? (snow or sea ice)
       sfc_is_frozen = xice(i,j) > 0.15 .OR. snowh(i,j) > 0.1

       !---- Initialize the seasonal category parameter
       ! This is only used to set up the interception parameters if the grid
       ! cell is vegetated, so no need to worry about "lush vegetation" at the
       ! poles or in other deserts
       !
       ! If there is snow on the ground or sea ice at the surface, use the
       ! winter/frozen ground seasonal category
       IF(sfc_is_frozen) THEN
          iseason = 4 ! winter, snow on ground
       ELSEIF(xlat(i,j) > 25.0) THEN
         !-- Northern midlatitudes and the Arctic
         IF(julday >= 59 .and. julday <= 149) THEN
           iseason = 5 ! Transitional spring with partially green short annuals.
         ELSEIF(julday > 149 .and. julday < 244) THEN
           iseason = 1 ! Midsummer with lush vegetation.
         ELSEIF(julday >= 244 .and. julday <= 304) THEN
           iseason = 2 ! Autumn with cropland that has not been harvested.
         ELSEIF((julday >=0 .and. julday < 59) &
                .or. (julday > 304 .and. julday <= 365)) THEN
           iseason = 3 ! Late autumn after frost, no snow.
         ELSE
           CALL wrf_error_fatal('aer_drydep_simple: Invalid julday')
         ENDIF
       ELSEIF(xlat(i,j) > -25.0) THEN
         !-- Tropical regions
         iseason = 1 ! Midsummer with lush vegetation.
       ELSE ! xlat(i,j) <= -25.0
         !-- Southern midlatitudes and Antarctica
         iseason = 2 ! Default: autumn with cropland that has not been harvested.
         IF(julday >= 239 .and. julday <= 329) THEN
           iseason = 5 ! Transitional spring with partially green short annuals.
         ELSEIF((julday >= 0 .and.  julday < 64) &
                .or. (julday > 329 .and. julday <= 365)) THEN
           iseason = 1 ! Midsummer with lush vegetation.
         ELSEIF(julday >= 64 .and. julday <= 124) THEN
           iseason = 2 ! Autumn with cropland that has not been harvested.
         ELSEIF (julday > 124 .and. julday < 239) THEN
           iseason = 3 ! Late autumn after frost, no snow.
         ELSE
           CALL wrf_error_fatal('aer_drydep_simple: Invalid julday')
         ENDIF
       ENDIF

       !---- Initialize the surface category parameters
       ! If there is snow on the ground or sea ice at the surface, use the
       ! snow/ice surface category instead
       IF(sfc_is_frozen) THEN
          zhangsurfcat = 12 ! ice
       ELSE IF ( trim(mminlu_loc) == 'USGS' ) THEN
         ! Map from usgs to zhang surface category
         IF ((wrfsurfcat >= 1) .and. (wrfsurfcat <= 25)) THEN
            zhangsurfcat = ZHANGSURF_FROM_USGS(wrfsurfcat)
         ELSE
           CALL wrf_error_fatal('aer_drydep_simple: Invalid wrfsurfcat')
         ENDIF
       ELSEIF( trim(mminlu_loc) == 'MODIFIED_IGBP_MODIS_NOAH' ) THEN
         ! Map from MODIS to zhang surface category
         IF ((wrfsurfcat >= 1) .and. (wrfsurfcat <= 20)) THEN
            zhangsurfcat = ZHANG_SURF_FROM_MODIS(wrfsurfcat)
         ELSE
           CALL wrf_error_fatal('aer_drydep_simple: Invalid wrfsurfcat')
         ENDIF
       ELSEIF( trim(mminlu_loc) == 'MODIS' ) THEN
         ! Map from MODIS+lakes to zhang surface category
         IF ((wrfsurfcat >= 1) .and. (wrfsurfcat <= 21)) THEN
            zhangsurfcat = ZHANG_SURF_FROM_MODIS(wrfsurfcat)
         ELSE
           CALL wrf_error_fatal('aer_drydep_simple: Invalid wrfsurfcat')
         ENDIF
       ELSE
           CALL wrf_error_fatal('aer_drydep_simple: Invalid mminlu_loc')
       ENDIF
       sfc_is_vegetated = zhangsurfcat /= 8 &
                     .AND. zhangsurfcat /= 9 &
                     .AND. zhangsurfcat /= 12 &
                     .AND. zhangsurfcat /= 13 &
                     .AND. zhangsurfcat /= 14

       !---- Land-surface dependent interception parameter converted to cm
       aintercept = RAD_COLLECTOR_ZHANG(iseason, zhangsurfcat) * 1.E2

       DO iaer = 1,naerspecies
         !---- Initialize aerosol parameters and convert units
         vdep = 0.0
         p_chem = aer_pointers(iaer)

         !-- Aerosol wet diameter (cm)
         aer_dry_d_cm = aer_dry_diam(iaer) * 1.0E2
         ! Calculate wet diameter with Gerber1985 eq. 15 and Table 1
         ! This should only be done for hydrophilic aerosols, sea salt, sulfate
         ! This could be enabled with oc2, bc2, using 'rural' parameters from
         ! Gerber1985
         IF(p_chem == p_seas_1 .OR. p_chem == p_seas_2 .OR. p_chem == p_seas_3 &
            .OR. p_chem == p_seas_4) THEN
           aer_wet_d_cm = ( C1_WET_SEASALT*aer_dry_d_cm**C2_WET_SEASALT &
                         / (C3_WET_SEASALT*aer_dry_d_cm**C4_WET_SEASALT - log(relhum)) &
                         + aer_dry_d_cm**3.0 ) ** (1./3.)
           ! WRITE(*,*) 'i,j,iaer,rh,aer_diam_dry,aer_diam_wet: ', i,j,iaer,relhum,aer_dry_diam(iaer),aer_wet_d_cm
         ELSEIF(p_chem == p_sulf) THEN
           aer_wet_d_cm = ( C1_WET_SULFATE*aer_dry_d_cm**C2_WET_SULFATE &
                         / (C3_WET_SULFATE*aer_dry_d_cm**C4_WET_SULFATE - log(relhum)) &
                         + aer_dry_d_cm**3.0 ) ** (1./3.)
           ! WRITE(*,*) 'i,j,iaer,rh,aer_diam_dry,aer_diam_wet: ', i,j,iaer,relhum,aer_dry_diam(iaer),aer_wet_d_cm
         ELSE
           aer_wet_d_cm = aer_dry_d_cm
         ENDIF

         !-- Aerosol wet density (g/cm3)
         ! Using the dry density for now
         aer_dry_dens_gcm3 = aer_dry_dens(iaer) * 1.E-3

         !---- Calculate the sedimentation velocity at the surface
         ! Settling (Stokes) velocity (cm/s)
         vsettl = (GRAVITY*aer_dry_dens_gcm3*aer_wet_d_cm*aer_wet_d_cm) &
                  / (18.*airkinvisc*rho_gcm3)
         ! Cunningham correction factor for the stokes velocity
         cunningham_fact = 1.0 + 2.0 * freepath/aer_wet_d_cm &
                           * (1.257 + 0.400*EXP(-0.55*aer_wet_d_cm/freepath) )
         vsettl = vsettl*cunningham_fact

         !---- Calculate the surface deposition resistance
         aerodiffus = BOLTZMANN * tk * cunningham_fact &
                      / (3.0*PI*airkinvisc*rho_gcm3*aer_wet_d_cm)
         schmidt = airkinvisc/aerodiffus
         IF(sfc_is_vegetated) THEN
           IF(aintercept < 0.) CALL wrf_error_fatal('aer_wetscav_simple: Bad interception coefficient')
           stokes = ust_cms*vsettl/(GRAVITY*aintercept)
           eintercept = CINTERCEPT * (aer_wet_d_cm/aintercept)**NUINTERCEPT
         ELSE
           stokes = ust_cms*ust_cms*vsettl/(GRAVITY*airkinvisc)
           eintercept = 0.0
         ENDIF
         eimpact = CIMPACT * (stokes/(ALPHAIMPACT+stokes))**BETAIMPACT
         ebrown = CBROWN * schmidt**(-2.0/3.0)
         r1fac = EXP(-(stokes**0.5)) ! Zhang2001
         res_sfc = 1. / (EPSZ * ust(i,j) * (ebrown+eimpact+eintercept) * r1fac) ! s/m

         !---- Calculate the total deposition velocity in m/s
         vsettl = vsettl * 0.01 ! m/s
         vdep = vsettl + 1./(res_adyn(i,j) + res_sfc)
         ! WRITE(*,*) 'i,j,iaer,aerd,aerdens,vdep,vsettl,res_sfc,res_adyn',i,j,iaer,aer_wet_d_cm,aer_dry_dens_gcm3,vdep,vsettl,res_sfc,res_adyn(i,j)
         IF(p_chem > param_first_scalar) THEN
           ddvel(i,j,p_chem) = vdep
         ENDIF

       END DO !iaer
       ! WRITE(*,*) 'i,j,iaer,vdep,vsettl,res_sfc,res_adyn',i,j,iaer,aer_wet_d_cm,aer_dry_dens_gcm3,vdep,vsettl,res_sfc,res_adyn(i,j)
     END DO ! i
   END DO ! j

  END SUBROUTINE aer_drydep_simple
  
!----------------------------------------------------------------------------

  SUBROUTINE aer_drydep_simple_init(chem_opt)
    ! Purpose
    !  Initialize aer_drydep_simple parameters, which depend on chem_opt
    !
    ! Arguments
    !  chem_opt = chemistry-aerosol option from namelist
    !
    ! References
    !  (1) For GOCART parameters:
    !      Chin et al. (2002) - Tropospheric Aerosol Optical Thickness from
    !      the GOCART Model and Comparisons with Satellite and Sun Photometer
    !      Measurements, J. Atmos. Sci., 59, 461–483,
    !      https://doi.org/10.1175/1520-0469(2002)059<0461:TAOTFT>2.0.CO;2

    USE module_state_description, ONLY:p_sulf, p_seas_1, p_seas_2, p_seas_3, &
      p_seas_4, p_bc1, p_oc1, p_bc2, p_oc2, p_dust_1, p_dust_2, p_dust_3, &
      p_dust_4, p_dust_5, MOZCART_KPP, GOCART_SIMPLE, T1_MOZCART_KPP, &
      GOCARTRADM2, GOCARTRACM_KPP

    IMPLICIT NONE

    ! Arguments:
    !  Input
    INTEGER, INTENT(IN) :: chem_opt

    ! Local variables
    INTEGER :: status_id


    !---- Initialize chem mechanism type
    aer_is_gocart = chem_opt == MOZCART_KPP .OR. &
                    chem_opt == T1_MOZCART_KPP .OR. &
                    chem_opt == GOCART_SIMPLE .OR. &
                    chem_opt == GOCARTRADM2 .OR. &
                    chem_opt == GOCARTRACM_KPP

    !-------- Initialize land cover type --------
    CALL nl_get_mminlu( 1, mminlu_loc )
    CALL wrf_debug(15, 'aer_drydep_simple_init: Using the following land cover type:')
    CALL wrf_debug(15, TRIM(mminlu_loc))

    !-------- Initialize aerosol tracer number for chem_opt --------
    ! This could be recovered from the aerosol scheme data directly
    IF(aer_is_gocart) THEN
      ! MOZCART
      naerspecies = 14
    ELSE
      CALL wrf_error_fatal ('aer_drydep_simple_init: chem_opt not supported')
    ENDIF

    !---- Allocate aerosol mappings and parameters for dry deposition
    IF (.NOT. ALLOCATED(aer_pointers)) THEN
      ALLOCATE(aer_pointers(naerspecies), stat=status_id)
      IF (status_id /= 0) CALL wrf_error_fatal('aer_drydep_simple_init: Could not allocate aer_pointers')
    ENDIF
    IF (.NOT. ALLOCATED(aer_dry_diam)) THEN
      ALLOCATE(aer_dry_diam(naerspecies), stat=status_id)
      IF (status_id /= 0) CALL wrf_error_fatal('aer_drydep_simple_init: Could not allocate aer_dry_diam')
    ENDIF
    IF (.NOT. ALLOCATED(aer_dry_dens)) THEN
      ALLOCATE(aer_dry_dens(naerspecies), stat=status_id)
      IF (status_id /= 0) CALL wrf_error_fatal('aer_drydep_simple_init: Could not allocate aer_dry_dens')
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
      CALL wrf_error_fatal ('aer_drydep_simple_init: chem_opt not supported')
    ENDIF

  END SUBROUTINE aer_drydep_simple_init

!----------------------------------------------------------------------------

END MODULE module_aer_drydep_simple



