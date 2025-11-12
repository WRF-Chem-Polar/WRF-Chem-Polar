! 2025/06 Louis Marelle
!
! Purpose:
!  Perform dry deposition and sedimentation for aerosols at the surface, with
!  dry deposition velocities calculated following Emerson et al. (2020), a
!  modified version of the Zhang et al. (2001) parameterization with updated
!  parameters.
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
!  role in radiative effect estimates, Proc. Natl. Acad. Sci. U.S.A.  117 (42)
!  26076-26082, https://doi.org/10.1073/pnas.2014761117
!  (2) Zhang et al., (2001) - A size-segregated particle dry deposition scheme
!  for an atmospheric aerosol module. Atmospheric environment, 35(3), 549-560,
!  https://doi.org/10.1016/S1352-2310(00)00326-5
!  (3) Gerber (1985) - Relative humidity parameterization of the Navy aerosol
!  model (NAM), NRL Rep. 8956. National Research Laboratory, Washington DC.
!  https://apps.dtic.mil/sti/tr/pdf/ADA163209.pdf

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
         9.999, 9.999, 9.999, 9.999, 9.999, &  ! luc 8
         9.999, 9.999, 9.999, 9.999, 9.999, &  ! luc 9
         0.010, 0.010, 0.010, 0.010, 0.010, &  ! luc 10
         0.010, 0.010, 0.010, 0.010, 0.010, &  ! luc 11
         9.999, 9.999, 9.999, 9.999, 9.999, &  ! luc 12
         9.999, 9.999, 9.999, 9.999, 9.999, &  ! luc 13
         9.999, 9.999, 9.999, 9.999, 9.999, &  ! luc 14
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

CONTAINS

!----------------------------------------------------------------------------

  SUBROUTINE aer_drydep_simple(chem, t_phy, rho_phy, ust, dz8w, rh, &
                               xlat, snowh, xice, ivgtyp,           &
                               res_adyn, julday, dtstep,            &
                               ids,ide, jds,jde, kds,kde,           &
                               ims,ime, jms,jme, kms,kme,           &
                               its,ite, jts,jte, kts,kte            )
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
   !    xlat = latitude (-90 to 90)
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
   USE module_model_constants, ONLY: g,piconst
   USE module_state_description, ONLY: p_sulf, p_seas_1, p_seas_2, p_seas_3, &
     p_seas_4, p_dust_1, p_dust_2, p_dust_3, p_dust_4, p_dust_5, num_chem

   IMPLICIT NONE

   !---- Arguments
   REAL, DIMENSION(ims:ime, kms:kme, jms:jme), INTENT(IN) :: t_phy, rho_phy, &
     dz8w, rh
   REAL, DIMENSION(ims:ime, jms:jme), INTENT(IN) :: ust, xlat, snowh, xice
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
   ! Dry deposition rate (s-1)
   REAL :: drydeprate
   ! Interception collection radius parameter (cm)
   REAL :: aintercept
   ! Deposition efficiencies from brownian motion, impaction, and interception
   REAL :: ebrown, eimpact, eintercept
   ! Particle bounce correction factor
   REAL :: r1fac
   ! x,y,z, and chemical species index, chemical species pointer
   INTEGER :: i, j, k, nv, p_chem
   ! Zhang2001 surface category
   INTEGER :: zhangsurfcat
   ! WRF (USGS or MODIS) surface category
   INTEGER :: wrfsurfcat
   ! Zhang2001 Seasonal category
   INTEGER :: iseason
   ! Flag for vegetated surfaces (true if vegetation is present)
   LOGICAL :: is_vegetated

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
   REAL, PARAMETER :: C3_WET_SEASALT = 2.573E-11
   REAL, PARAMETER :: C4_WET_SEASALT = -1.424
   REAL, PARAMETER :: C1_WET_SULFATE = 0.4809
   REAL, PARAMETER :: C2_WET_SULFATE = 3.082
   REAL, PARAMETER :: C3_WET_SULFATE = 3.110E-11
   REAL, PARAMETER :: C4_WET_SULFATE = -1.428

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

       !---- Initialize the seasonal category parameter
       ! This is only used to set up the interception parameters if the grid
       ! cell is vegetated, so no need to worry about "lush vegetation" at the
       ! poles or in other deserts
       IF(xlat(i,j) > 25.0) THEN
         !-- Northern midlatitudes and the Arctic
         iseason = 2 ! Default: autumn with cropland that has not been harvested.
         IF(julday >= 60 .and. julday <= 150) THEN
           iseason = 5 ! Transitional spring with partially green short annuals.
         ELSEIF(julday > 150 .and. julday < 245) THEN
           iseason = 1 ! Midsummer with lush vegetation.
         ELSEIF(julday >= 245 .and. julday <= 305) THEN
           iseason = 2 ! Autumn with cropland that has not been harvested.
         ELSEIF (julday < 60 .or. julday > 305) THEN
           iseason = 3 ! Late autumn after frost, no snow.
         ENDIF
       ELSEIF(xlat(i,j) > -25.0) THEN
         !-- Tropical regions
         iseason = 1 ! Midsummer with lush vegetation.
       ELSE ! xlat(i,j) <= -25.0
         !-- Southern midlatitudes and Antarctica
         iseason = 2 ! Default: autumn with cropland that has not been harvested.
         IF(julday >= 240 .and. julday <= 330) THEN
           iseason = 5 ! Transitional spring with partially green short annuals.
         ELSEIF(julday > 330 .or. julday < 65) THEN
           iseason = 1 ! Midsummer with lush vegetation.
         ELSEIF(julday >= 65 .and. julday <= 125) THEN
           iseason = 2 ! Autumn with cropland that has not been harvested.
         ELSEIF (julday > 125 .and. julday < 240) THEN
           iseason = 3 ! Late autumn after frost, no snow.
         ENDIF
       ENDIF
       ! If there is snow on the ground or sea ice at the surface, use the
       ! winter/frozen ground seasonal category
       IF(xice(i,j) > 0.15 .OR. snowh(i,j) > 0.1) THEN
          iseason = 4 ! winter, snow on ground
       ENDIF

       !---- Initialize the surface category parameters
       zhangsurfcat = 12 ! default to ice
       IF ( trim(mminlu_loc) == 'USGS' ) THEN
         ! Map from usgs to zhang surface category
         IF ((wrfsurfcat >= 1) .and. (wrfsurfcat <= 25)) &
            zhangsurfcat = ZHANGSURF_FROM_USGS(wrfsurfcat)
       ENDIF
       IF( trim(mminlu_loc) == 'MODIFIED_IGBP_MODIS_NOAH' ) THEN
         ! Map from MODIS to zhang surface category
         IF ((wrfsurfcat >= 1) .and. (wrfsurfcat <= 20)) &
            zhangsurfcat = ZHANG_SURF_FROM_MODIS(wrfsurfcat)
       ENDIF
       IF( trim(mminlu_loc) == 'MODIS' ) THEN
         ! Map from MODIS+lakes to zhang surface category
         IF ((wrfsurfcat >= 1) .and. (wrfsurfcat <= 21)) &
            zhangsurfcat = ZHANG_SURF_FROM_MODIS(wrfsurfcat)
       ENDIF
       ! If there is snow on the ground or sea ice at the surface, use the
       ! snow/ice surface category
       IF(xice(i,j) > 0.15 .OR. snowh(i,j) > 0.1) THEN
          zhangsurfcat = 12 ! ice
       ENDIF
       is_vegetated = zhangsurfcat == 8 &
                      .OR. zhangsurfcat == 9 &
                      .OR. zhangsurfcat == 12 &
                      .OR. zhangsurfcat == 13 &
                      .OR. zhangsurfcat == 14

       !---- Land-surface dependent interception parameter converted to cm
       aintercept = RAD_COLLECTOR_ZHANG(iseason, zhangsurfcat) * 1.E2

       DO nv = 1,naerspecies
         !---- Initialize aerosol parameters and convert units
         vdep = 0.0
         p_chem = aer_pointers(nv)

         !-- Aerosol wet diameter (cm)
         aer_dry_d_cm = aer_dry_diam(nv) * 1.0E2
         ! Calculate wet diameter with eq. 10 in Zhang2001, also see Gerber1985
         ! This should only be done for hydrophilic aerosols, sea salt, sulfate
         ! This could be enabled with oc2, bc2, using 'rural' parameters from Zhang2001
         IF(p_chem == p_seas_1 .OR. p_chem == p_seas_2 .OR. p_chem == p_seas_3 .OR. p_chem == p_seas_4) THEN
           aer_wet_d_cm = ( C1_WET_SEASALT*aer_dry_d_cm**C2_WET_SEASALT &
                         / (C3_WET_SEASALT*aer_dry_d_cm**C4_WET_SEASALT - log(relhum)) &
                         + aer_dry_d_cm**3.0 ) ** (1./3.)
           ! WRITE(*,*) 'i,j,nv,rh,aer_diam_dry,aer_diam_wet: ', i,j,nv,relhum,aer_dry_diam(nv),aer_wet_d_cm
         ELSEIF(p_chem == p_sulf) THEN
           aer_wet_d_cm = ( C1_WET_SULFATE*aer_dry_d_cm**C2_WET_SULFATE &
                         / (C3_WET_SULFATE*aer_dry_d_cm**C4_WET_SULFATE - log(relhum)) &
                         + aer_dry_d_cm**3.0 ) ** (1./3.)
           ! WRITE(*,*) 'i,j,nv,rh,aer_diam_dry,aer_diam_wet: ', i,j,nv,relhum,aer_dry_diam(nv),aer_wet_d_cm
         ENDIF

         !-- Aerosol wet density (g/cm3)
         ! Using the dry density for now
         aer_dry_dens_gcm3 = aer_dry_dens(nv) * 1.E-3

         !---- Calculate the sedimentation velocity at the surface
         ! Settling (Stokes) velocity (cm/s)
         vsettl = (GRAVITY*aer_dry_dens_gcm3*aer_wet_d_cm*aer_wet_d_cm) &
                  / (18.*airkinvisc*rho_gcm3)
         ! Cunningham correction factor for the stokes velocity
         cunningham_fact = 1.0 + 2.0 * freepath/aer_wet_d_cm &
                           * (1.257 + 0.400*EXP(-0.55*aer_wet_d_cm/freepath) )
         vsettl = vsettl*cunningham_fact

         !---- Calculate the surface deposition resistance
         ! g cm2 s-2 K-1 / *** -> g cm2 s-2 / ***
         ! -> g cm2 s-2 / (cm2 s-1 g cm-3 * cm) ->  s-1 / (cm-2) -> cm2/s
         aerodiffus = BOLTZMANN * tk * cunningham_fact &
                      / (3.0*PI*airkinvisc*rho_gcm3*aer_wet_d_cm)
         schmidt = airkinvisc/aerodiffus
         IF(is_vegetated) THEN
           stokes = ust_cms*vsettl/(GRAVITY*aintercept)
           eintercept = CINTERCEPT * (aer_wet_d_cm/aintercept)**NUINTERCEPT
         ELSE
           stokes = ust_cms*ust_cms*vsettl/(GRAVITY*airkinvisc)
           eintercept = 0.0
         ENDIF
         eimpact = CIMPACT * (stokes/(ALPHAIMPACT+stokes))**BETAIMPACT
         ebrown = CBROWN * schmidt**(-0.66666666)
         r1fac = EXP(-(stokes**0.5)) ! Zhang2001
         res_sfc = 1. / (EPSZ * ust(i,j) * (ebrown+eimpact+eintercept) * r1fac) ! s/m

         !---- Calculate the total deposition velocity in m/s
         vsettl = vsettl * 0.01 ! m/s
         vdep = vsettl + 1./(res_adyn(i,j) + res_sfc)

         !---- Perform dry deposition
         ! WRITE(*,*) 'i,j,nv,aerd,aerdens,vdep,vsettl,res_sfc,res_adyn',i,j,nv,aer_wet_d_cm,aer_dry_dens_gcm3,vdep,vsettl,res_sfc,res_adyn(i,j)
         drydeprate = vdep/dz8w(i,kts,j) ! s -1
         aer_mmr_old = chem(i,kts,j,p_chem)
         aer_mmr = aer_mmr_old * EXP(-dtstep*drydeprate)
         chem(i,kts,j,p_chem) = MAX(0.0, aer_mmr)

       END DO !nv
       ! WRITE(*,*) 'i,j,nv,vdep,vsettl,res_sfc,res_adyn',i,j,nv,aer_wet_d_cm,aer_dry_dens_gcm3,vdep,vsettl,res_sfc,res_adyn(i,j)
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

    USE module_state_description, ONLY:p_sulf, p_seas_1, p_seas_2, p_seas_3, &
      p_seas_4, p_bc1, p_oc1, p_bc2, p_oc2, p_dust_1, p_dust_2, p_dust_3, &
      p_dust_4, p_dust_5, MOZCART_KPP, GOCART_SIMPLE

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
      ! Effective aerosol dry diameter
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



