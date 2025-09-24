! 2025/06 Louis Marelle
!
! Purpose:
!  Perform dry deposition and sedimentation for aerosols
!
! Remarks:
!  - Only compatible with GOCART for now
!  - Sedimentation above the first level is done in gocart_settling or mosaic_settling
!  - Hygroscopic growth is only taken into account for seasalt and sulfate
!
! References:
! ==============================================================================
!  (1) Emerson et al., 2020 - https://doi.org/10.1073/pnas.2014761117
!  (2) Zhang et al., 2001 - https://doi.org/10.1016/S1352-2310(00)00326-5
!  (3) Gerber, 1985, https://apps.dtic.mil/sti/tr/pdf/ADA163209.pdf


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
  CHARACTER*256 :: mminlu_loc
  ! Radius (m) of surface collectors for interception

  !TODO following parameters could be moved to a data module
  ! Interception efficiency land-dependent coefficient, from Zhang2001
  REAL, SAVE :: rad_collector_zhang(5,15)
  DATA rad_collector_zhang/ &
         0.002, 0.002, 0.002, 0.002, 0.002, &  ! luc 1
         0.005, 0.005, 0.005, 0.005, 0.005, &  ! luc 2
         0.002, 0.002, 0.005, 0.005, 0.002, &  ! luc 3
         0.005, 0.005, 0.010, 0.010, 0.005, &  ! luc 4
         0.005, 0.005, 0.005, 0.005, 0.005, &  ! luc 5
         0.002, 0.002, 0.005, 0.005, 0.002, &  ! luc 6
         0.002, 0.002, 0.005, 0.005, 0.002, &  ! luc 7
         9.999, 9.999, 9.999, 9.999, 9.999, &  ! luc 8
         9.999, 9.999, 9.999, 9.999, 9.999, &  ! luc 9
         0.010, 0.010, 0.010, 0.010, 0.010, &  ! luc 10
         0.010, 0.010, 0.010, 0.010, 0.010, &  ! luc 11
         9.999, 9.999, 9.999, 9.999, 9.999, &  ! luc 12
         9.999, 9.999, 9.999, 9.999, 9.999, &  ! luc 13
         9.999, 9.999, 9.999, 9.999, 9.999, &  ! luc 14
         0.010, 0.010, 0.010, 0.010, 0.010  /  ! luc 15

  ! Maps MODIS surface categories to zhang categories
  INTEGER, SAVE  :: zhangsurf_from_modis(21)
  DATA zhangsurf_from_modis / &
    !   ! zhang category        <---- modis category
    1, &! evergreen needle-leaf    1 'Evergreen Needleleaf Forest'
    2, &! evergreen broadleaf      2 'Evergreen Broadleaf Forest'
    3, &! deciduous needle-leaf    3 'Deciduous Needleleaf Forest'
    4, &! deciduous broadleaf      4 'Deciduous Broadleaf Forest'
    5, &! mixed broad/needle-leaf  5 'Mixed Forests'
    6, &! grass                    6 'Closed Shrublands'
    6, &! grass                    7 'Open Shrublands'
    6, &! grass                    8 'Woody Savannas'
    6, &! grass                    9 'Savannas'
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
  INTEGER, SAVE  :: zhangsurf_from_usgs(25)
  DATA zhangsurf_from_usgs / &
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

CONTAINS

!----------------------------------------------------------------------------

  SUBROUTINE aer_drydep_simple(chem, t_phy, rho_phy, ust, dz8w, rh, &
                               snowh, xice, ivgtyp, &
                               res_adyn, julday, dtstep, &
                               ids,ide, jds,jde, kds,kde,                   &
                               ims,ime, jms,jme, kms,kme,                   &
                               its,ite, jts,jte, kts,kte                    )
   ! Purpose:
   !  Calculates and applies dry deposition tendencies for aerosols
   !  
   ! Arguments:
   !  Input:
   !    t_phy = air temperature (degree celsius)
   !    rho_phy = air density (kg/m3)
   !    ust = friction velocity (m/s)
   !    dz8w = model vertical level depth (m)
   !    rh = relative humidity (0-1)
   !    snowh = snow depth (m)
   !    xice = sea ice concentration (fraction)
   !    ivgtyp = land use category (categorical)
   !    res_adyn = aerodynamical deposition resistance from wesely (s/m)
   !    julday = julian day (1-366)
   !    dtstep = model time step (s)
   !    ids,ims,its etc. = domain, memory, and tile start and end indices
   !  Input/output:
   !    chem = advected chemical species (gases in ppm, aerosols in ug/kg)

   ! USE associations:
   !TODO Use ONLY:
   USE module_configure
   USE module_state_description

   IMPLICIT NONE

   !---- Arguments
   REAL, DIMENSION(ims:ime, kms:kme, jms:jme), INTENT(IN) :: t_phy
   REAL, DIMENSION(ims:ime, kms:kme, jms:jme), INTENT(IN) :: rho_phy
   REAL, DIMENSION(ims:ime, jms:jme), INTENT(IN) :: ust
   REAL, DIMENSION(ims:ime, kms:kme, jms:jme), INTENT(IN) :: dz8w
   REAL, DIMENSION(ims:ime, kms:kme, jms:jme), INTENT(IN) :: rh
   REAL, DIMENSION(ims:ime, jms:jme), INTENT(IN) :: snowh
   REAL, DIMENSION(ims:ime, jms:jme), INTENT(IN) :: xice
   INTEGER, DIMENSION(ims:ime, jms:jme), INTENT(IN) :: ivgtyp
   REAL, DIMENSION(its:ite, jts:jte), INTENT(IN) :: res_adyn
   REAL, DIMENSION(ims:ime, kms:kme, jms:jme, num_chem),          &
       INTENT(INOUT) :: chem
   REAL, INTENT(IN) :: dtstep
   INTEGER, INTENT(IN) :: julday
   INTEGER, INTENT(IN) :: ids,ide, jds,jde, kds,kde,    &
                          ims,ime, jms,jme, kms,kme,    &
                          its,ite, jts,jte, kts,kte

   !---- Local variables
   ! Aerosol mass mixing ratio (microg/kg)
   REAL :: aer_mmr, aer_mmr_old
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
   ! Aerosol wet density (g/cm3)
   REAL :: aerodens
   ! Stokes and schmidt numbers (dimensionless)
   REAL :: stokes, schmidt
   ! Particle Brownian diffusivity (cm2/s)
   REAL :: aerodiffus
   ! Aerosol wet diameter (cm)
   REAL :: aer_diam_cm
   ! Air kinematic viscosity (cm2/s)
   REAL :: airkinvisc
   ! Cunningham correction factor
   REAL :: cunningham_fact
   ! Dry deposition rate (s-1)
   REAL :: drydeprate
   ! Interception collection radius parameter (cm)
   REAL :: aintercept
   ! Deposition efficiencies from brownian motion, impaction, and interception
   REAL :: ebrown, eimpact, eintercept
   ! Particle bounce correction factor
   REAL :: r1fac
   !TODO describe
   REAL :: c1_wet,c2_wet,c3_wet,c4_wet
   ! x,y,z, and chemical species index, chemical species pointer
   INTEGER :: i, j, k, nv, p_chem
   !TODO
   INTEGER :: zhangsurfcat, wrfsurfcat, iseason, iseason_in
   LOGICAL :: is_vegetated

   !---- Parameters
   ! Surface resistance parameterization empirical coefficient
   !TODO This is varying by process in Zhang2001
   REAL, PARAMETER :: EPSZ = 3.
   ! Parameters for the deposition efficiencies, from Emerson2020
   REAL, PARAMETER :: cbrown = 0.2
   REAL, PARAMETER :: cimpact = 0.4, alphaimpact = 0.8, betaimpact = 1.7
   REAL, PARAMETER :: cintercept=2.5, nuintercept=0.8
   ! REAL, PARAMETER :: boltzmann = 1.3807E-23 ! J/K = kg m-2 s-2 K-1
   REAL, PARAMETER :: boltzmann = 1.3807E-16 ! erg/K = g cm2 s-2 K-1
   REAL, PARAMETER :: pi = 3.1416
   REAL, PARAMETER :: gravity = 980.616 ! cm s-2

   !---- Initialize seasonal category parameter
   !TODO this is hardcoded for the northern hemisphere midlatitudes
   !TODO read in XLAT
   ! This should not cause too many issues in the high NH latitudes but will
   ! cause problems in the southern hemisphere and the tropics
   iseason_in = 2 ! Default: autumn with cropland that has not been harvested.
   IF(julday >= 60 .and. julday <= 150) THEN
     iseason_in = 5 ! Transitional spring with partially green short annuals.
   ELSEIF(julday > 151 .and. julday < 245) THEN
     iseason_in = 1 ! Midsummer with lush vegetation.
   ELSEIF(julday >= 245 .and. julday < 305) THEN
     iseason_in = 2 ! Autumn with cropland that has not been harvested.
   ELSEIF (julday < 60 .or. julday > 300) THEN
     iseason_in = 3 ! Late autumn after frost, no snow.
   ENDIF
   ! category 4 (Winter, snow on ground and sub-freezing.) can be set below
   ! depending on snow/ice cover

   !-------- Perform dry deposition for aerosol species -------
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

       !---- Get the surface and seasonal category
       ! Get the surface category
       zhangsurfcat = 12 ! default to ice
       IF ( trim(mminlu_loc) == 'USGS' ) THEN
         ! Map from usgs to zhang surface category
         IF ((wrfsurfcat >= 1) .and. (wrfsurfcat <= 25)) &
            zhangsurfcat = zhangsurf_from_usgs(wrfsurfcat)
       ENDIF
       IF( trim(mminlu_loc) == 'MODIFIED_IGBP_MODIS_NOAH' ) THEN
         ! Map from MODIS to zhang surface category
         IF ((wrfsurfcat >= 1) .and. (wrfsurfcat <= 20)) &
            zhangsurfcat = zhangsurf_from_modis(wrfsurfcat)
       ENDIF
       IF( trim(mminlu_loc) == 'MODIS' ) THEN
         ! Map from MODIS+lakes to zhang surface category
         IF ((wrfsurfcat >= 1) .and. (wrfsurfcat <= 21)) &
            zhangsurfcat = zhangsurf_from_modis(wrfsurfcat)
       ENDIF
       is_vegetated = zhangsurfcat == 8 &
                      .OR. zhangsurfcat == 9 &
                      .OR. zhangsurfcat == 12 &
                      .OR. zhangsurfcat == 13 &
                      .OR. zhangsurfcat == 14
       ! Get the seasonal category
       iseason = iseason_in
       ! Force snow/ice surface and winter if there is snow on ground or sea
       ! ice
       IF(xice(i,j) > 0.15 .OR. snowh(i,j) > 0.1) THEN
          zhangsurfcat = 12 ! ice
          iseason = 4 ! winter, snow on ground
       ENDIF

       !---- Land-surface dependent interception parameter converted to cm
       aintercept = rad_collector_zhang(iseason, zhangsurfcat) * 1.E2

       DO nv = 1,naerspecies
         !---- Initialize aerosol parameters and convert units
         vdep = 0.0
         p_chem = aer_pointers(nv)

         !-- Aerosol wet diameter (cm)
         aer_diam_cm = aer_dry_diam(nv) * 1.0E2
         ! Calculate wet diameter with eq. 10 in Zhang2001, also see Gerber1985
         !TODO only for hygrophilic aerosols, sea salt, sulfate, and oc2, bc2
         IF(p_chem == p_seas_1 .OR. p_chem == p_seas_2 .OR. p_chem == p_seas_3 .OR. p_chem == p_seas_4) THEN
           !TODO move to data module?
           c1_wet = 0.7674
           c2_wet = 3.079
           c3_wet = 2.573E-11
           c4_wet = -1.424
           aer_diam_cm = ( c1_wet*aer_diam_cm**c2_wet &
                         / (c3_wet*aer_diam_cm**c4_wet - log(relhum)) &
                         + aer_diam_cm**3.0 ) ** (1./3.)
           ! WRITE(*,*) 'i,j,nv,rh,aer_diam_dry,aer_diam_wet: ', i,j,nv,relhum,aer_dry_diam(nv),aer_diam_cm
         ELSEIF(p_chem == p_sulf) THEN
           !TODO move to data module?
           c1_wet = 0.4809
           c2_wet = 3.082
           c3_wet = 3.110E-11
           c4_wet = -1.428
           aer_diam_cm = ( c1_wet*aer_diam_cm**c2_wet &
                         / (c3_wet*aer_diam_cm**c4_wet - log(relhum)) &
                         + aer_diam_cm**3.0 ) ** (1./3.)
           ! WRITE(*,*) 'i,j,nv,rh,aer_diam_dry,aer_diam_wet: ', i,j,nv,relhum,aer_dry_diam(nv),aer_diam_cm
         ENDIF

         !-- Aerosol wet density (g/cm3)
         !TODO using dry density for now
         aerodens = aer_dry_dens(nv) * 1.E-3

         !---- Calculate the sedimentation velocity at the surface
         ! Settling (Stokes) velocity (cm/s)
         vsettl = (gravity*aerodens*aer_diam_cm*aer_diam_cm) &
                  / (18.*airkinvisc*rho_gcm3)
         ! Cunningham correction factor for the stokes velocity
         cunningham_fact = 1.0 + 2.0 * freepath/aer_diam_cm &
                           * (1.257 + 0.400*EXP(-0.55*aer_diam_cm/freepath) )
         vsettl = vsettl*cunningham_fact

         !---- Calculate the surface deposition resistance
         !TODO Check zhang for modal factors etc.
         ! g cm2 s-2 K-1 / *** -> g cm2 s-2 / ***
         ! -> g cm2 s-2 / (cm2 s-1 g cm-3 * cm) ->  s-1 / (cm-2) -> cm2/s
         aerodiffus = boltzmann * tk * cunningham_fact &
                      / (3.0*pi*airkinvisc*rho_gcm3*aer_diam_cm)
         schmidt = airkinvisc/aerodiffus
         IF(is_vegetated) THEN
           stokes = ust_cms*vsettl/(gravity*aintercept)
           eintercept = cintercept * (aer_diam_cm/aintercept)**nuintercept
         ELSE
           stokes = ust_cms*ust_cms*vsettl/(gravity*airkinvisc)
           eintercept = 0.0
         ENDIF
         eimpact = cimpact * (stokes/(alphaimpact+stokes))**betaimpact
         ebrown = cbrown * schmidt**(-0.66666666)
         r1fac = EXP(-(stokes**0.5)) ! Zhang2001
         res_sfc = 1. / (EPSZ * ust(i,j) * (ebrown+eimpact+eintercept) * r1fac) ! s/m

         !---- Calculate the total deposition velocity in m/s
         vsettl = vsettl * 0.01 ! m/s
         vdep = vsettl + 1./(res_adyn(i,j) + res_sfc)

         !---- Perform dry deposition
         ! WRITE(*,*) 'i,j,nv,aerd,aerdens,vdep,vsettl,res_sfc,res_adyn',i,j,nv,aer_diam_cm,aerodens,vdep,vsettl,res_sfc,res_adyn(i,j)
         drydeprate = vdep/dz8w(i,kts,j) ! s -1
         aer_mmr_old = chem(i,kts,j,p_chem)
         aer_mmr = aer_mmr_old * EXP(-dtstep*drydeprate)
         chem(i,kts,j,p_chem) = MAX(0.0, aer_mmr)

       END DO !nv
       ! WRITE(*,*) 'i,j,nv,vdep,vsettl,res_sfc,res_adyn',i,j,nv,aer_diam_cm,aerodens,vdep,vsettl,res_sfc,res_adyn(i,j)
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
    !  (1) GOCART parameters: Chin et al. (2002) -
    !      https://doi.org/10.1175/1520-0469(2002)059%3C0461:TAOTFT%3E2.0.CO;2

    USE module_state_description ! For chem_opt and p_ option values

    IMPLICIT NONE

    ! Arguments:
    !  Input
    INTEGER, INTENT(IN) :: chem_opt

    !-------- Initialize land cover type --------
    CALL nl_get_mminlu( 1, mminlu_loc )
    CALL wrf_debug(15, 'drydep_simple_init: Using the following land cover type:')
    CALL wrf_debug(15, TRIM(mminlu_loc))

    !-------- Initialize aerosol tracer number for chem_opt --------
    ! This could be recovered from the aerosol scheme data directly
    IF(chem_opt == MOZCART_KPP .OR. chem_opt == GOCART_SIMPLE) THEN
      ! MOZCART
      naerspecies = 14
    ELSE
      CALL wrf_error_fatal ('drydep_simple_init: chem_opt not supported')
    ENDIF

    !---- Allocate aerosol mappings and parameters for dry deposition
    IF (.NOT. ALLOCATED(aer_pointers)) THEN
      ALLOCATE(aer_pointers(naerspecies))
    ENDIF
    IF (.NOT. ALLOCATED(aer_dry_diam)) THEN
      ALLOCATE(aer_dry_diam(naerspecies))
    ENDIF
    IF (.NOT. ALLOCATED(aer_dry_dens)) THEN
      ALLOCATE(aer_dry_dens(naerspecies))
    ENDIF

    !---- Initialize aerosol mappings and parameters for dry deposition
    IF(chem_opt == MOZCART_KPP .OR. chem_opt == GOCART_SIMPLE) THEN
      !-- GOCART
      ! Species that experience dry deposition
      aer_pointers(:) = (/ p_sulf, p_seas_1, p_seas_2, p_seas_3, p_seas_4, &
                           p_bc1, p_oc1, p_bc2, p_oc2, p_dust_1, &
                           p_dust_2, p_dust_3, p_dust_4, p_dust_5 /)
      !TODO I used the effective diameter
      aer_dry_diam(:) = (/ 0.312E-6, 0.6E-6, 2.0E-6, 6.5E-6, 15.0E-6, &
                           0.078E-6, 0.174E-6, 0.078E-6, 0.174E-6, 1.46E-6, &
                           2.8E-6, 4.8E-6, 9.0E-6, 16.0E-6 /)

      aer_dry_dens(:) = (/ 1700., 2200., 2200., 2200., 2200., &
                           1000., 1800., 1000., 1800., 2600., &
                           2600., 2600., 2600., 2600. /)
    ELSE
      CALL wrf_error_fatal ('wetscav_simple_init: chem_opt not supported')
    ENDIF

  END SUBROUTINE aer_drydep_simple_init

!----------------------------------------------------------------------------

END MODULE module_aer_drydep_simple



