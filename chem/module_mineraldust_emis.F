! Louis Marelle 2025/03
!
! Purpose:
!  Compute size-resolved mineral dust aerosol emissions (mass and number) for the
!  sectional 4-bin or 8-bin MOSAIC aerosols, or for the GOCART 5 dust bins.
!
! Remarks:
!  TODO should include a subgrid wind speed parameterization
!  TODO FLEXDUST ignores MB95 basal surface and soil type
!  TODO check G04 differences with AFWA
!
! References
! ==============================================================================
! Source functions:
!  (1) Ginoux et al. (2004) - doi:10.1016/S1364-8152(03)00114-2
!  (2) Ginoux et al. (2001) - https://doi.org/10.1029/2000JD000053
!  (3) Marticorena and Bergametti (1995) - https://doi.org/10.1029/95JD00690
!  (4) FLEXDUST - Groot Zwaaftink et al. (2016) - http://dx.doi.org/10.1002/2016JD025482
!  (5) AFWA - LeGrand et al. (2019) - https://doi.org/10.5194/gmd-12-131-2019
!      (describes the original wrfgocart parameterizations and its issues)
! Size distributions:
!  (1) Kok (2010) - https://doi.org/10.1073/pnas.1014798108
! Soil moisture correction factor
!  (1) Fécan et al. (1999) - https://doi.org/10.1007/s00585-999-0149-7

MODULE module_mineraldust_emis

  IMPLICIT NONE

  !---- Global module-wide variables
  ! These are initialized in the init subroutine

  ! Allocatable size-resolved aerosol arrays (dim is size bin number and depends
  ! on aerosol mechanism)
  REAL, DIMENSION(:), ALLOCATABLE :: bin_diameters_cen
  ! REAL, DIMENSION(:), ALLOCATABLE :: bin_diameters_low
  ! REAL, DIMENSION(:), ALLOCATABLE :: bin_diameters_high
  REAL, DIMENSION(:), ALLOCATABLE :: emiss_size_massfrac

  ! Number of dust aerosol bins for chem_opt in use
  INTEGER :: naersizes

  ! Aerosol density (kg/m3)
  REAL :: aer_dens

  ! chem_opt aerosol model family
  LOGICAL :: aer_is_gocart
  LOGICAL :: aer_is_mosaic4bin
  LOGICAL :: aer_is_mosaic8bin

  ! Soil properties for calculating the saltation flux, from AFWA
  INTEGER, PARAMETER :: NSOIL=9
  REAL, DIMENSION (NSOIL), PARAMETER :: REFF_SOIL=(/0.71E-6,1.37E-6,2.63E-6,5.00E-6,9.50E-6,18.1E-6,34.5E-6,65.5E-6,125.E-6/) ! m
  REAL, DIMENSION (NSOIL), PARAMETER :: DEN_SOIL=(/2500.,2650.,2650.,2650.,2650.,2650.,2650.,2650.,2650./) ! kg/m3
  REAL, DIMENSION (NSOIL), PARAMETER :: FRACSIZE_SOIL=(/1.,0.2,0.2,0.2,0.2,0.2,0.333,0.333,0.333/) ! fraction, 0-1
  INTEGER, DIMENSION (NSOIL), PARAMETER :: SPOINT_SOIL=(/1,2,2,2,2,2,3,3,3/)  ! 1 Clay, 2 Silt, 3 Sand

  !TODO also ensable for USGS
  INTEGER, PARAMETER :: ISICE_MODIS = 12

  CONTAINS

!------------------------------------------------------------------------------

  SUBROUTINE mineraldust_emis( id, dtstep, dz8w, &
         chem_opt, dust_opt, dust_hl_opt,   &
         rho_phy, chem, e_dust,                                     &
         u10, v10, ust, xland, tsk, xice, snowh, &
         clayfrac, sandfrac, &
         erod, erod_hl, vegfra, smois_top, isltyp, ivgtyp, &
         ids,ide, jds,jde, kds,kde,                                 &
         ims,ime, jms,jme, kms,kme,                                 &
         its,ite, jts,jte, kts,kte                                  )

    ! Purpose:
    !  Compute size-resolved mineral dust aerosol number and mass emissions

    ! Arguments:
    ! Input:
    !  id = domain ID number
    !  dtstep = model time step (s)
    !  dz8w = 3D vertical level depth (m)
    !  chem_opt = chemical mechanism + aerosol mechanism option from config_flags
    !  dust_opt = dust emission option from config_flags
    !  dust_hl_opt = high-latitude dust option from config_flags
    !  rho_phy = 3D air density (kg/m3)
    !  u10 = windpeed at 10m in x (positive west-east) direction (m/s)
    !  v10 = windpeed at 10m in y (positive south-north) direction (m/s)
    !  ust = friction velocity (m/s)
    !  xland = landmask flag (1 land, 2 water)
    !  tsk = surface skin temperature (K)
    !  xice = sea ice or lake ice fraction (0-1)
    !  snowh = snow height (m)
    !  clayfrac = soil clay fraction (0-1)
    !  sandfrac = soil sand fraction (0-1)
    !  erod = soil erodibility fraction (0-1)
    !  erod_hl = soil erodibility including high latitude sources (0-1)
    !  vegfra = vegetation fractional cover in grid cell (0-1)
    !  smois_top = soil moisture in top soil layer (m3/m3)
    !  isltyp = soil type category (int category)
    !  ivgtyp = land type category (int category)
    !  ids,ims,its etc. = domain, memory, and tile start and end indices
    ! Input/output
    !  chem = advected chemical species (gases in ppm, aerosols in microg/kg)
    !  e_dust = accumulated dust emission flux (mg/m2)

    !TODO Should include all explicit ONLY: associations for module_state_description
    !TODO chem_opt, dust_opt, param_first_scalar, num_chem, p_dust_1, p_dust_2, p_dust_3, p_dust_4, p_dust_5
    USE module_state_description
    USE module_model_constants, ONLY: piconst
    USE module_data_gocart_dust, ONLY: porosity
    !TODO Is it possible to only include this for mosaic?
    USE module_data_mosaic_asect, ONLY: numptr_aer, lptr_oin_aer, nsize_aer

    IMPLICIT NONE

    !---- Arguments
    ! Input
    INTEGER, INTENT(IN) :: chem_opt, dust_opt, dust_hl_opt
    INTEGER, INTENT(IN) :: id,                                      &
                           ids,ide, jds,jde, kds,kde,               &
                           ims,ime, jms,jme, kms,kme,               &
                           its,ite, jts,jte, kts,kte
    INTEGER, DIMENSION(ims:ime, jms:jme), INTENT(IN) :: isltyp, ivgtyp
    REAL, INTENT(IN) :: dtstep
    REAL, DIMENSION(ims:ime, kms:kme, jms:jme), INTENT(IN) :: dz8w, rho_phy
    REAL, DIMENSION(ims:ime, jms:jme), INTENT(IN) :: u10, v10, ust, xland
    REAL, DIMENSION(ims:ime, jms:jme), INTENT(IN) :: tsk, xice, snowh
    REAL, DIMENSION(ims:ime, jms:jme), INTENT(IN) :: smois_top, vegfra, clayfrac, sandfrac
    REAL, DIMENSION(ims:ime, jms:jme, 3), INTENT(IN) :: erod
    REAL, DIMENSION(ims:ime, jms:jme), INTENT(IN) :: erod_hl
    ! Input/output
    REAL, DIMENSION(ims:ime, kms:kme, jms:jme, num_chem), INTENT(INOUT) :: chem
    REAL, DIMENSION(ims:ime, jms:jme), INTENT(INOUT) :: e_dust

    !---- Local variables
    ! Dust emission fluxes (microg/m2/s, microm3_m2_s, #/m2/s)
    REAL :: dust_mass_emiss, dust_vol_emiss, dust_num_emiss
    ! Bulk dust mass emission flux (microg/m2/s)
    REAL :: dust_mass_emiss_tot
    ! 10-m wind speed (m/s)
    REAL :: w10
    ! Friction velocity (m/s)
    REAL :: ust_in
    ! Clay and sand fraction in soil (0-1)
    REAL :: clayfrac_in, sandfrac_in
    ! Conversion factor from 1/m2/s to 1/kg_dry_air
    REAL :: conv
    ! Surface wetness (0-1)
    REAL :: gwet
    ! Gravimetric soil water content (kg/kg)
    REAL :: smois_grav
    ! Aerosol diameter (m) and air density (kg/m3)
    REAL :: aer_diam, air_dens
    ! Soil grain diameter (m) and density (kg/m3)
    REAL :: soil_diam, soil_dens
    ! Fraction of soil type (sand, clay or silt) in cell (0-1)
    REAL :: soiltype_frac
    ! Indices
    INTEGER :: i, j, isize, iphase
    INTEGER :: iaerbin, isoilbin
    INTEGER :: l_dust, l_num
    ! Erodibility scaled by bare soil fraction
    REAL, DIMENSION (3) :: erod_in
    ! REAL :: dust_s_dist_mb95, stot_mb95, d_aer_diam
    ! Size distribution of dust in the original wrfgocart scheme
    REAL, DIMENSION(5) :: frac_s_wrfgocart = (/0.1,0.25,0.25,0.25,0.25/)
    ! ! Reduce the erodibility in line with the original WRFGOCART scheme (see
    ! ! AFWA publication)
    ! REAL, PARAMETER :: erod_hl_scaling_factor = 0.25
    ! Disable the scaling
    REAL, PARAMETER :: erod_hl_scaling_factor = 1.
    ! ! Scale the emissions by 50% (what is suggested by test simulations TODO)
    ! ! TODO apply to G04 only
    ! REAL, PARAMETER :: erod_hl_scaling_factor = 0.5

    !-------- Compute dust emissions --------
    ! No emissions to cloud-phase
    ! This should be = aiphase but it might not be initialized
    iphase = 1

    DO i=its,ite
      DO j=jts,jte

        ! Only perform dust emissions over snow-free land. Sea ice is
        ! considered land, so exclude sea ice points explicitly.
        IF (xland(i,j)<1.5 .AND. snowh(i,j)<0.1 .AND. xice(i,j)<=0.01) THEN
          !TODO Could also cut emissions when precipitation occurs, based on
          !e.g. precr+preci+precs+precg at bottom level
          ! Soil erodibility
          IF(dust_hl_opt == 1) THEN
            ! Scale the erodibility by (1.-vegfra): only bare soil can emit
            ! dust. This is not already done for HLdust EROD
            erod_in(1) = erod_hl(i,j) * (100.-vegfra(i,j))/100. &
                         * erod_hl_scaling_factor
            erod_in(2) = erod_in(1)
            erod_in(3) = erod_in(1)
          ELSE
            ! WPS EROD already includes bare soil scaling
            erod_in(1) = erod(i,j,1)
            erod_in(2) = erod(i,j,2)
            erod_in(3) = erod(i,j,3)
          ENDIF

          ! Skip dust emission calculation if the soil is not erodible
          IF(erod_in(1) > 1.E-3) THEN
            !-------- Initialize mineral dust emission parameters
            dust_mass_emiss = 0.0
            dust_num_emiss = 0.0
            w10 = SQRT(u10(i,j)**2.0 + v10(i,j)**2.0)
            ust_in = ust(i,j)
            clayfrac_in = clayfrac(i,j)
            sandfrac_in = sandfrac(i,j)
            air_dens = rho_phy(i,kts,j)
            ! Conversion factor from 1/m2/s to 1/kg_dry_air
            conv = 1.0/(rho_phy(i,kts,j) * dz8w(i,kts,j)) * dtstep
            ! Volumetric soil moisture over porosity
            !TODO is isltyp invariant in all LSM setups? Porosity is mapped to
            !     specific soil categories
            gwet = smois_top(i,j)/porosity(isltyp(i,j))
            gwet = MAX(1.E-3, MIN(1.0, gwet))
            IF(dust_hl_opt == 1) THEN
              !TODO only works with MODIS, add USGS support
              IF(ivgtyp(i,j) == ISICE_MODIS) THEN
                ! Do not use the soil water to limit emissions for land ice
                ! grid points if we want high latitude dust emissions.
                ! If a fraction of a land ice grid point is erodible it should
                ! emit dust if the other conditions on wind speed and
                ! snow cover are met. Emissions will be limited by the Fecan
                ! factor based on smois_grav but will not be cut to 0.
                gwet = 0.0
              ENDIF
            ENDIF
            ! Gravimetric soil water content
            smois_grav = smois_top(i,j)/((1.-porosity(isltyp(i,j)))*(2.65*(1-clayfrac_in)+2.50*clayfrac_in))

            ! !---- Calculate total basal surface for MB95
            !!TODO Check if this is still compatible with the new soil size distribution
            ! IF(dust_opt == DUST_MB95)
            !   stot_mb95 = 0.0
            !   DO iaerbin = 1,naersizes
            !     ! Basal surface fraction from MB1995 eq. 30 (TODO m2/kg/m ?)
            !     !Here I normalize by Mtot - ok because the point is to get
            !     !dS/Stot and dS is also normalized
            !     dust_s_dist_mb95 = emiss_size_massfrac(iaerbin) &
            !                        / (2./3. * aer_diam * aer_dens)
            !     d_aer_diam = bin_diameters_high(iaerbin)*1.0E-6 - bin_diameters_low(iaerbin)*1.0E-6
            !     stot_mb95 = stot_mb95 + dust_s_dist_mb95 * d_aer_diam
            !   ENDDO
            ! ENDIF

            !-------- Compute the total aerosol mass emission flux
            dust_mass_emiss_tot = 0.0
            SELECT CASE(dust_opt)
            CASE(DUST_WRFGOCART)
              !---- Compute dust mass emissions with the original WRF-Chem GOCART scheme (based on 4.3.3 version)
              DO iaerbin = 1,naersizes
                aer_diam = bin_diameters_cen(iaerbin)*1.0E-6
                CALL source_mineraldust_gocart(dust_mass_emiss, w10, gwet, &
                    erod_in, aer_dens, aer_diam, air_dens)
                dust_mass_emiss_tot = dust_mass_emiss_tot &
                    + dust_mass_emiss*frac_s_wrfgocart(iaerbin)
              ENDDO

            CASE(DUST_G04)
              !---- Compute dust mass emissions following Ginoux et al. (2004)
              !TODO should redefine this as a new option, since this is now
              !     quite different from the original implementation
              !     How can I name this option?
              !TODO maybe G04_UPDATED for now
              DO isoilbin = 1,NSOIL
                soil_diam = REFF_SOIL(isoilbin)*2.0
                soil_dens = DEN_SOIL(isoilbin)
                IF(SPOINT_SOIL(isoilbin) == 1) THEN
                  ! Clay
                  soiltype_frac = clayfrac_in
                ELSEIF(SPOINT_SOIL(isoilbin) == 2) THEN
                  ! Silt
                  soiltype_frac = MAX(0., 1.-clayfrac_in-sandfrac_in)
                ELSEIF(SPOINT_SOIL(isoilbin) == 3) THEN
                  ! Sand
                  soiltype_frac = sandfrac_in
                ENDIF
                CALL source_mineraldust_ginoux2004(dust_mass_emiss, w10, &
                    ust_in, smois_grav, gwet, erod_in, clayfrac_in, soil_dens, &
                    soil_diam, air_dens)
                dust_mass_emiss_tot = dust_mass_emiss_tot &
                    + dust_mass_emiss*soiltype_frac*FRACSIZE_SOIL(isoilbin)
              ENDDO ! isoilbin

            CASE(DUST_G01)
              !---- Compute dust mass emissions following Ginoux et al. (2001)
              DO isoilbin = 1,NSOIL
                soil_diam = REFF_SOIL(isoilbin)*2.0
                soil_dens = DEN_SOIL(isoilbin)
                IF(SPOINT_SOIL(isoilbin) == 1) THEN
                  ! Clay
                  soiltype_frac = clayfrac_in
                ELSEIF(SPOINT_SOIL(isoilbin) == 2) THEN
                  ! Silt
                  soiltype_frac = MAX(0., 1.-clayfrac_in-sandfrac_in)
                ELSEIF(SPOINT_SOIL(isoilbin) == 3) THEN
                  ! Sand
                  soiltype_frac = sandfrac_in
                ENDIF
              CALL source_mineraldust_ginoux2001(dust_mass_emiss, w10, gwet, &
                  erod_in, soil_dens, soil_diam, air_dens)
              dust_mass_emiss_tot = dust_mass_emiss_tot &
                  + dust_mass_emiss*soiltype_frac*FRACSIZE_SOIL(isoilbin)
              ENDDO

            CASE(DUST_FLEXDUST)
              !---- Compute dust mass emissions following Groot Zwaaftink et al. (2016)
              soil_dens = aer_dens
              CALL source_mineraldust_flexdust(dust_mass_emiss, w10, ust_in, &
                  smois_grav, erod_in, clayfrac_in, sandfrac_in, soil_dens, &
                  air_dens)
              dust_mass_emiss_tot = dust_mass_emiss

            CASE(DUST_MB95)
              !---- Compute dust mass emissions following Marticorena and
              !---- Bergametti (1995)
              DO isoilbin = 1,NSOIL
                soil_diam = REFF_SOIL(isoilbin)*2.0
                soil_dens = DEN_SOIL(isoilbin)
                IF(SPOINT_SOIL(isoilbin) == 1) THEN
                  ! Clay
                  soiltype_frac = clayfrac_in
                ELSEIF(SPOINT_SOIL(isoilbin) == 2) THEN
                  ! Silt
                  soiltype_frac = MAX(0., 1.-clayfrac_in-sandfrac_in)
                ELSEIF(SPOINT_SOIL(isoilbin) == 3) THEN
                  ! Sand
                  soiltype_frac = sandfrac_in
                ENDIF
                CALL source_mineraldust_mb95(dust_mass_emiss, ust_in, &
                         smois_grav, erod_in, clayfrac_in, soil_dens, &
                         soil_diam, air_dens)
                dust_mass_emiss_tot = dust_mass_emiss_tot &
                    + dust_mass_emiss*soiltype_frac*FRACSIZE_SOIL(isoilbin)
              ENDDO
              !TODO Should integrate on basal surface rather than on soil frac, see MB95
              !!TODO Check if this is still compatible with the new soil size distribution
              ! ! Basal surface fraction from MB1995 eq. 30 (TODO m2/kg/m ?)
              ! dust_s_dist_mb95 = emiss_size_massfrac(iaerbin) &
              !                    / (2./3. * aer_diam * aer_dens)
              ! ! Integrate dust mass emissions over basal surface
              ! d_aer_diam = bin_diameters_high(iaerbin)*1.0E-6 - bin_diameters_low(iaerbin)*1.0E-6
              ! dust_mass_emiss_tot = dust_mass_emiss_tot + dust_mass_emiss *  dust_s_dist_mb95/stot_mb95 * d_aer_diam

            CASE DEFAULT
              !TODO could also implement other schemes, e.g. UOC, Tegen-Vignati from EC-EARTH4
              CALL wrf_error_fatal ('mineraldust_emis: dust_opt not supported')

            END SELECT

            !-------- Emit mineral dust aerosol mass and number for each size bin
            DO iaerbin = 1,naersizes
              !---- Distribute total dust emissions to size bin
              IF(dust_opt == DUST_WRFGOCART) THEN
                ! Special case for WRFGOCART for compatibility with previous version
                dust_mass_emiss = dust_mass_emiss_tot * frac_s_wrfgocart(iaerbin)
              ELSE
                ! For other dust schemes distribute to each bin based on emiss_size_massfrac
                dust_mass_emiss = dust_mass_emiss_tot * emiss_size_massfrac(iaerbin)
              ENDIF
              e_dust(i,j) = e_dust(i,j) + dust_mass_emiss*1.E-3

              !---- Convert dust aerosol mass emission flux to number
              ! emission flux (#/m2/s)
              ! 1.E9 is for mm3 to um3 conversion
              dust_vol_emiss = dust_mass_emiss / aer_dens * 1.E9
              dust_num_emiss = dust_vol_emiss / ((piconst/6.0)*bin_diameters_cen(iaerbin)**3.0)
              dust_num_emiss = MAX(dust_num_emiss, 0.0)

              !-------- Emit mineral dust emissions to chem(i,k,j,l)
              !---- Emit MOSAIC aerosols to chem
              IF(aer_is_mosaic4bin .OR. aer_is_mosaic8bin) THEN
                ! MOSAIC dust aerosol
                l_dust = lptr_oin_aer(iaerbin,1,iphase)
                IF (l_dust >= param_first_scalar) THEN
                  chem(i,kts,j,l_dust) = chem(i,kts,j,l_dust) &
                                         + dust_mass_emiss * conv
                END IF
                ! MOSAIC aerosol number
                l_num = numptr_aer(iaerbin,1,iphase)
                IF (l_num >= param_first_scalar) THEN
                  chem(i,kts,j,l_num) = chem(i,kts,j,l_num) + dust_num_emiss * conv
                END IF
              !---- Emit GOCART aerosols to chem
              ELSEIF(aer_is_gocart) THEN
                ! GOCART dust aerosol
                IF(iaerbin == 1 ) THEN
                  l_dust = p_dust_1
                ELSEIF(iaerbin == 2 ) THEN
                  l_dust = p_dust_2
                ELSEIF(iaerbin == 3 ) THEN
                  l_dust = p_dust_3
                ELSEIF(iaerbin == 4 ) THEN
                  l_dust = p_dust_4
                ELSEIF(iaerbin == 5 ) THEN
                  l_dust = p_dust_5
                ENDIF
                IF (l_dust >= param_first_scalar) THEN
                  chem(i,kts,j,l_dust) = chem(i,kts,j,l_dust) &
                                         + dust_mass_emiss * conv
                ENDIF
              !---- Other chem_opt cases are not supported
              ELSE
                CALL wrf_error_fatal ('mineraldust_emis: chem_opt not supported')
              ENDIF ! aer_is_
            END DO ! iaerbin

          END IF ! erod_in
        END IF ! xland
      END DO ! j
    END DO ! i

  END SUBROUTINE mineraldust_emis

!------------------------------------------------------------------------------

  SUBROUTINE source_mineraldust_gocart(dust_mass_emiss, w10, gwet, erod, aer_dens, aer_diam, air_dens)
    ! Purpose:
    !  Compute dust mass emission flux (ug/m2/s) as in the GOCART emissions scheme in WRF-Chem
    !
    ! Remarks:
    !  Only use this for GOCART mechanisms

    ! Arguments:
    ! Input:
    !  w10 = 10-m wind speed (m/s)
    !  gwet = soil moisture (0-1)
    !  erod = soil erodibility (0-1)
    !  aer_dens = aerosol density (kg/m3)
    !  aer_diam = aerosol diameter (m)
    !  air_dens = air density (kg/m3)
    ! Output:
    !  dust_mass_emiss = dust aerosol mass emission flux (ug/m2/s)

    USE module_model_constants, ONLY: g

    IMPLICIT NONE

    !---- Arguments
    ! Input
    REAL, INTENT(IN) :: w10, gwet, aer_dens, aer_diam, air_dens
    REAL, DIMENSION(3), INTENT(IN) :: erod
    ! Output
    REAL, INTENT(OUT) :: dust_mass_emiss

    !---- Local variables
    ! Ginoux source function term (soil erodibility, 0-1)
    REAL :: source_ginoux
    ! Aerosol density, GOCART values (kg/m3)
    REAL :: aer_dens_gocart
    ! Threshold velocity for wind erosion (m/s)
    REAL :: w_ts

    !---- Parameters
    ! Ginoux dust emissions scaling parameter, ug/s2/m5
    REAL, PARAMETER :: cdim_ginoux = 1.0

    ! Ginoux source function S = erod(2) for bins 2-5, erod(3) for bin
    ! 1 in GOCART module
    IF(aer_diam < 2.0E-6) THEN
      source_ginoux = erod(3)
      aer_dens_gocart = 2500.
    ELSE
      source_ginoux = erod(2)
      aer_dens_gocart = 2650.
    ENDIF
    ! Threshold velocity (m/s)
    w_ts = 0.13*1.0E-2*SQRT((aer_dens_gocart*1.E-3)*(g*1.E2)*(aer_diam*1.E2)/(air_dens*1.E-3))* &
           SQRT(1.0+0.006/(aer_dens_gocart*1.E-3)/(g*1.E2)/(aer_diam*1.E2)**2.5)/ &
           SQRT(1.928*(1331.0*(aer_diam*1.E2)**1.56+0.38)**0.092-1.0)
    ! Correct threshold velocity for soil moisture
    IF (gwet < 0.5) THEN
      w_ts = w_ts * (1.2 + 0.2*LOG10(gwet))
      w_ts = MAX(0.0, w_ts)
    ELSE
      w_ts = 1000.0
    ENDIF
    ! Calculate the dust emissions (ug/m2/s)
    IF(w10 > w_ts) THEN
      ! 0.8 factor is ch_dust from chemics_init
      dust_mass_emiss = 0.8 * cdim_ginoux * source_ginoux * w10**2 * (w10 - w_ts)
    ELSE
      dust_mass_emiss = 0.0
    ENDIF
    dust_mass_emiss = MAX(0.0, dust_mass_emiss)

  END SUBROUTINE source_mineraldust_gocart

!------------------------------------------------------------------------------

  SUBROUTINE source_mineraldust_mb95(dust_mass_emiss, ust, smois_grav, erod, clayfrac, soil_dens, soil_diam, air_dens)
    ! Purpose:
    !  Compute dust mass emission flux (ug/m2/s) following Marticorena and Bergametti (1995)

    ! References:
    !  (1) Marticorena and Bergametti (1995) - https://doi.org/10.1029/95JD00690
    ! Soil moisture correction factor
    !  (1) Fécan et al. (1999) - https://doi.org/10.1007/s00585-999-0149-7

    ! Arguments:
    ! Input:
    !  ust = friction velocity (m/s)
    !  smois_grav = gravimetric soil water content (kg/kg)
    !  erod = soil erodibility (0-1)
    !  clayfrac = soil clay fraction (0-1)
    !  soil_dens = soil density (kg/m3)
    !  soil_diam = soil grain diameter (m)
    !  air_dens = air density (kg/m3)
    ! Output:
    !  dust_mass_emiss = dust aerosol mass emission flux (ug/m2/s)

    USE module_model_constants, ONLY: g

    IMPLICIT NONE

    !---- Arguments
    ! Input
    REAL, INTENT(IN) :: ust, smois_grav, clayfrac, soil_dens, soil_diam, air_dens
    REAL, DIMENSION(3), INTENT(IN) :: erod
    ! Output
    REAL, INTENT(OUT) :: dust_mass_emiss

    !---- Local variables
    ! Source function term (soil erodibility, 0-1)
    REAL :: source_erod
    ! Threshold velocity for wind erosion (m/s)
    REAL :: w_ts
    ! Dust horizontal flux from saltation (kg/cm/s)
    REAL :: dust_horiz_flux_g
    ! Dust vertical emission flux (kg/cm2/s)
    REAL :: dust_vert_flux_f
    ! Ratio of the vertical and horizontal flux (cm-1)
    REAL :: f_g_ratio
    ! Adjusted soil clay fraction (%)
    REAL :: clayfrac_adj
    ! Threshold gravimetric water content for threshold velocity correction (kg/kg)
    REAL :: smois_thresh

    source_erod = erod(2)
    clayfrac_adj = MIN(clayfrac*100., 50.)
    ! Threshold velocity from MB1995, converted to IS units
    w_ts = 1.E-2 * 0.129 * SQRT(soil_dens/air_dens * g * soil_diam) \
           * SQRT(1. + 6.E-7/(soil_dens*g*(soil_diam)**2.5)) \
           / SQRT(1.928 * (1754600.*(soil_diam)**1.56 + 0.38)**0.092 - 1.)
    ! Correct threshold velocity for soil moisture (Fécan et al., 1999)
    !TODO cut emissions completely in very wet soils?
    smois_thresh = 1.E-2 * (0.17*(clayfrac_adj) + 0.0014*(clayfrac_adj)**2.0)
    IF(smois_grav > smois_thresh) THEN
      w_ts = w_ts * SQRT(1. + 1.21*(smois_grav - smois_thresh)**0.68)
    ENDIF
    IF(ust > w_ts) THEN
      ! Dust horizontal flux from MB1995, (G in eq. 28, kg/cm/s)
      dust_horiz_flux_g = 2.61 * air_dens*(1.E-6)/g * (ust*1.E2)**3.0 * (1. + w_ts/ust)\
                          * (1. - w_ts**2.0/ust**2.0)
      ! MB1995 eq. 47 - clayfrac 0% - 20%
      clayfrac_adj = MIN(clayfrac*100., 20.)
      f_g_ratio = 10.**(0.134 * clayfrac_adj - 6.) ! cm-1
      dust_vert_flux_f = dust_horiz_flux_g * f_g_ratio ! kg/cm2/s
      dust_mass_emiss = dust_vert_flux_f*1.E5 ! kg/cm2/s -> ug/m2/s
    ELSE
      dust_mass_emiss = 0.0
    ENDIF
    dust_mass_emiss = MAX(0.0, dust_mass_emiss)

  END SUBROUTINE source_mineraldust_mb95

!------------------------------------------------------------------------------

  SUBROUTINE source_mineraldust_ginoux2001(dust_mass_emiss, w10, gwet, erod, soil_dens, soil_diam, air_dens)
    ! Purpose:
    !  Compute dust mass emission flux (ug/m2/s) following Ginoux et al. (2001)

    ! References:
    !  Ginoux et al. (2001) - https://doi.org/10.1029/2000JD000053

    ! Arguments:
    ! Input:
    !  w10 = 10-m wind speed (m/s)
    !  gwet = soil moisture (0-1)
    !  erod = soil erodibility (0-1)
    !  soil_dens = soil density (kg/m3)
    !  soil_diam = soil grain diameter (m)
    !  air_dens = air density (kg/m3)
    ! Output:
    !  dust_mass_emiss = dust aerosol mass emission flux (ug/m2/s)

    USE module_model_constants, ONLY: g

    IMPLICIT NONE

    !---- Arguments
    ! Input
    REAL, INTENT(IN) :: w10, gwet, soil_dens, soil_diam, air_dens
    REAL, DIMENSION(3), INTENT(IN) :: erod
    ! Output
    REAL, INTENT(OUT) :: dust_mass_emiss

    !---- Local variables
    ! Ginoux source function term (soil erodibility, 0-1)
    REAL :: source_ginoux
    ! Threshold velocity for wind erosion (m/s)
    REAL :: w_ts

    !---- Parameters
    ! Ginoux dust emissions scaling parameter, ug/s2/m5
    REAL, PARAMETER :: cdim_ginoux = 1.0

    source_ginoux = erod(2)
    ! Threshold velocity
    IF (gwet < 0.5) THEN
      w_ts = 6.5 * SQRT((soil_dens-air_dens)/air_dens * g * soil_diam) * (1.2 + 0.2*LOG10(gwet))
      w_ts = MAX(0.0, w_ts)
    ELSE
      w_ts = 1000.0
    ENDIF
    ! Calculate the dust emissions (ug/m2/s)
    IF(w10 > w_ts) THEN
      dust_mass_emiss = cdim_ginoux * source_ginoux * w10**2 * (w10 - w_ts)
    ELSE
      dust_mass_emiss = 0.0
    ENDIF
    dust_mass_emiss = MAX(0.0, dust_mass_emiss)

  END SUBROUTINE source_mineraldust_ginoux2001


!------------------------------------------------------------------------------

  SUBROUTINE source_mineraldust_ginoux2004(dust_mass_emiss, w10, ust, &
      smois_grav, gwet, erod, clayfrac, soil_dens, soil_diam, air_dens)
    ! Purpose:
    !  Compute dust mass emission flux (ug/m2/s) following Ginoux et al. (2004)

    ! References:
    !  (1) Ginoux et al. (2004) - doi:10.1016/S1364-8152(03)00114-2
    !  (2) Marticorena and Bergametti (1995) - https://doi.org/10.1029/95JD00690
    !  (3) LeGrand et al. (2019) - https://doi.org/10.5194/gmd-12-131-2019
    ! Soil moisture correction factor
    !  (1) Fécan et al. (1999) - https://doi.org/10.1007/s00585-999-0149-7

    ! Arguments:
    ! Input:
    !  w10 = 10-m wind speed (m/s)
    !  ust = friction velocity (m/s)
    !  smois_grav = gravimetric soil water content (kg/kg)
    !  gwet = soil moisture (0-1)
    !  erod = soil erodibility (0-1)
    !  clayfrac = soil clay fraction (0-1)
    !  soil_dens = soil density (kg/m3)
    !  soil_diam = soil grain diameter (m)
    !  air_dens = air density (kg/m3)
    ! Output:
    !  dust_mass_emiss = dust aerosol mass emission flux (ug/m2/s)

    USE module_model_constants, ONLY: g

    IMPLICIT NONE

    !---- Arguments
    ! Input
    REAL, INTENT(IN) :: w10, ust, smois_grav, gwet, clayfrac, soil_dens, soil_diam, air_dens
    REAL, DIMENSION(3), INTENT(IN) :: erod
    ! Output
    REAL, INTENT(OUT) :: dust_mass_emiss

    !---- Local variables
    ! Ginoux source function term (soil erodibility, 0-1)
    REAL :: source_ginoux
    ! Threshold velocity for wind erosion (m/s)
    REAL :: w_ts
    ! Adjusted soil clay fraction (%)
    REAL :: clayfrac_adj
    ! Threshold gravimetric water content for threshold velocity correction (kg/kg)
    REAL :: smois_thresh

    !---- Parameters
    ! Ginoux dust emissions scaling parameter, ug/s2/m5
    REAL, PARAMETER :: cdim_ginoux = 1.0

    source_ginoux = erod(2)

    ! Threshold velocity from Marticorena and Bergametti, 1995, converted to IS
    ! units
    w_ts = 1.E-2 * 0.129 * SQRT(soil_dens/air_dens * g * soil_diam) \
           * SQRT(1. + 6.E-7/(soil_dens*g*(soil_diam)**2.5)) \
           / SQRT(1.928 * (1754600.*(soil_diam)**1.56 + 0.38)**0.092 - 1.)
    ! ! Correct threshold velocity for soil moisture
    ! IF (gwet < 0.5) THEN
    !   w_ts = w_ts * (1.2 + 0.2*LOG10(gwet))
    !   w_ts = MAX(0.0, w_ts)
    ! ELSE
    !   w_ts = 1000.0
    ! ENDIF
    clayfrac_adj = MIN(clayfrac*100., 50.)
    ! Correct threshold velocity for soil moisture
    IF (gwet < 0.5) THEN
      ! Correct threshold velocity for soil moisture (Fécan et al., 1999)
      ! This is not the original Ginoux correction factor - the reason why is
      ! described in the AFWA paper
      smois_thresh = 1.E-2 * (0.17*(clayfrac_adj) + 0.0014*(clayfrac_adj)**2.0)
      IF(smois_grav > smois_thresh) THEN
        w_ts = w_ts * SQRT(1. + 1.21*(smois_grav - smois_thresh)**0.68)
      ENDIF
    ELSE
      w_ts = 1000.0
    ENDIF
    ! Correct u* threshold velocity to w10 threshold velocity
    IF(w10 > w_ts) THEN ! To avoid low w10 and ust values
      w_ts = w_ts * w10/ust
    ENDIF
    ! Calculate the dust emissions (ug/m2/s)
    IF(w10 > w_ts) THEN
      ! 0.76 factor from the sum of sp fractions in Ginoux et al. (2004)
      dust_mass_emiss = 0.76 * cdim_ginoux * source_ginoux * w10**2 * (w10 - w_ts)
    ELSE
      dust_mass_emiss = 0.0
    ENDIF
    dust_mass_emiss = MAX(0.0, dust_mass_emiss)

  END SUBROUTINE source_mineraldust_ginoux2004

!------------------------------------------------------------------------------

  SUBROUTINE source_mineraldust_flexdust(dust_mass_emiss, w10, ust, smois_grav, erod, clayfrac, sandfrac, soil_dens, air_dens)
    ! Purpose:
    !  Compute dust mass emission flux (ug/m2/s) as in the FLEXDUST model

    ! References:
    !  (1) Groot Zwaaftink et al. (2016) - http://dx.doi.org/10.1002/2016JD025482
    !  (2) Marticorena and Bergametti (1995) - https://doi.org/10.1029/95JD00690
    ! Soil moisture correction factor
    !  (1) Fécan et al. (1999) - https://doi.org/10.1007/s00585-999-0149-7

    ! Arguments:
    ! Input:
    !  w10 = 10-m wind speed (m/s)
    !  ust = friction velocity (m/s)
    !  smois_grav = gravimetric soil water content (kg/kg)
    !  erod = soil erodibility (0-1)
    !  clayfrac = soil clay fraction (0-1)
    !  sandfrac = soil sand fraction (0-1)
    !  soil_dens = soil grain density (kg/m3)
    !  air_dens = air density (kg/m3)
    ! Output:
    !  dust_mass_emiss = dust aerosol mass emission flux (ug/m2/s)

    USE module_model_constants, ONLY: g

    IMPLICIT NONE

    !---- Arguments
    ! Input
    REAL, INTENT(IN) :: w10, ust, smois_grav, soil_dens, air_dens
    REAL, INTENT(IN) :: clayfrac, sandfrac
    REAL, DIMENSION(3), INTENT(IN) :: erod
    ! Output
    REAL, INTENT(OUT) :: dust_mass_emiss

    !---- Local variables
    ! Ginoux source function term (soil erodibility, 0-1)
    REAL :: source_ginoux
    ! Threshold velocity for wind erosion (m/s)
    REAL :: w_ts
    ! Soil particle diameter (m)
    REAL :: soil_diam
    ! Threshold gravimetric water content for threshold velocity correction (kg/kg)
    REAL :: smois_thresh
    ! Sandblasting efficiency for emissions (m-1)
    REAL :: alpha_sandblast
    ! Adjusted clay and sand fraction (%)
    REAL :: clayfrac_adj, sandfrac_adj

    !---- Parameters
    REAL, PARAMETER :: sandfrac_minval = 1. ! %
    REAL, PARAMETER :: scaling_factor_mb95 = 4.8E-4

    !----
    source_ginoux = erod(2)
    clayfrac_adj = MIN(clayfrac*100., 50.)
    sandfrac_adj = sandfrac*100.
    ! Threshold velocity (GZ2016 eq. 2)
    IF(sandfrac_adj >= sandfrac_minval) THEN
      soil_diam = 75.E-6
    ELSE
      soil_diam = 10.E-6
    ENDIF
    w_ts = SQRT(0.0123 * (soil_dens/air_dens * g * soil_diam + 3E-4/(air_dens*soil_diam)))
    ! Correct threshold velocity for soil moisture (Fécan et al., 1999)
    smois_thresh = 1.E-2 * (0.17*(clayfrac_adj) + 0.0014*(clayfrac_adj)**2.0)
    IF(smois_grav > smois_thresh) THEN
      w_ts = w_ts * SQRT(1. + 1.21*(smois_grav - smois_thresh)**0.68)
    ENDIF
    ! Dust emissions from Marticorena and Bergametti (1995) (ug/m2/s)
    clayfrac_adj = MIN(clayfrac_adj, 20.)
    alpha_sandblast = 100. * 10.**(0.134*clayfrac_adj-6.) ! F/G ratio in MB95
    IF(ust > w_ts) THEN
      ! I could not find 4.8E-4 scaling factor from FLEXDUST scheme in MB95
      dust_mass_emiss = 1.E9 * source_ginoux * scaling_factor_mb95 &
                       * alpha_sandblast * air_dens * ust**3. / g &
                       * (1. - w_ts**2.0/ust**2.0) * (1. + w_ts/ust)
    ELSE
      dust_mass_emiss = 0.0
    ENDIF
    dust_mass_emiss = MAX(0.0, dust_mass_emiss)

  END SUBROUTINE source_mineraldust_flexdust

!------------------------------------------------------------------------------

  SUBROUTINE mineraldust_emis_init(chem_opt, dust_opt)

    ! Purpose:
    !  Initialize the parameters for the dust aerosol source functions used
    !  in mineraldust_emis, which depend on chem_opt and dust_opt

    USE module_state_description ! For chem_opt and dust_opt option values

    IMPLICIT NONE

    ! Arguments:
    !  Input
    INTEGER, INTENT(IN) :: chem_opt, dust_opt

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
                        chem_opt == CRI_MOSAIC_4BIN_AQ_KPP .OR. &
                        chem_opt == SAPRC99_MOSAIC_4BIN_VBS2_KPP

    aer_is_mosaic8bin = chem_opt == CBMZ_MOSAIC_8BIN .OR. &
                        chem_opt == CBMZ_MOSAIC_DMS_8BIN .OR. &
                        chem_opt == CBMZ_MOSAIC_8BIN_AQ .OR. &
                        chem_opt == CBMZ_MOSAIC_DMS_8BIN_AQ .OR. &
                        chem_opt == CRI_MOSAIC_8BIN_AQ_KPP .OR. &
                        chem_opt == SAPRC99_MOSAIC_8BIN_VBS2_AQ_KPP .OR. &
                        chem_opt == SAPRC99_MOSAIC_8BIN_VBS2_KPP

    !-------- Initialize mineral dust aerosol bin number for chem_opt --------
    ! This could be recovered from the aerosol scheme data directly, but
    ! emiss_size_massfrac is already chem_opt-dependent and depends on naersizes so
    ! there is no harm in redefining it here just for this module
    IF(aer_is_gocart) THEN
      ! GOCART with 5 dust bins cases
      naersizes = 5
    ELSEIF(aer_is_mosaic8bin) THEN
      ! MOSAIC-8bin cases
      naersizes = 8
    ELSEIF(aer_is_mosaic4bin) THEN
      ! MOSAIC-4bin cases
      naersizes = 4
    ELSE
      CALL wrf_error_fatal ('mineraldust_emis_init: chem_opt not supported')
    ENDIF

    !-------- Initialize size-resolved aerosol arrays for dust emissions --------
    ! Initialize emiss_size_massfrac (microm3/m2/s) and bin_diameters_cen
    ! (micrometer) for each family of aerosol options selected by
    ! chem_opt

    ! Allocate the size-resolved aerosol arrays
    IF (.NOT. ALLOCATED(bin_diameters_cen)) THEN
      ALLOCATE(bin_diameters_cen(naersizes))
    ENDIF
    IF (.NOT. ALLOCATED(emiss_size_massfrac)) THEN
      ALLOCATE(emiss_size_massfrac(naersizes))
    ENDIF

    !-- Initialize dust source functions for GOCART aerosols (5 dust bins)
    IF(aer_is_gocart) THEN
      ! GOCART dust bin center diameters (micrometer)
      ! These should be log-averaged but I follow the GOCART definition instead
      bin_diameters_cen(:) = (/ 1.46, 2.8, 4.8, 9.0, 16.0 /)
      ! Size-resolved emission scaling factors for GOCART
      ! https://doi.org/10.1073/pnas.1014798108
      emiss_size_massfrac(:) = (/ 0.00281, 0.178, 0.288, 0.349, 0.183 /)
      ! ! Value from original GOCART dust scheme
      ! emiss_size_massfrac = (/ 0.1, 0.25, 0.25, 0.25, 0.25 /)
      aer_dens = 2650.

    !-- Initialize dust source functions for MOSAIC-8bin aerosols
    ELSEIF(aer_is_mosaic8bin) THEN
      ! MOSAIC-8bin bin center diameters (micrometer)
      bin_diameters_cen(:) = (/ 0.0552, 0.1105, 0.2210, 0.4419, &
                                0.8839, 1.7678, 3.5355, 7.0711 /)
      ! Size-resolved emission scaling factors for MOSAIC-8bin
      ! https://doi.org/10.1073/pnas.1014798108
      emiss_size_massfrac(:) = (/ 1.93E-05, 2.40E-04, 2.04E-03, 1.20E-02,&
                                  4.97E-02, 1.46E-01, 3.12E-01, 4.78E-01 /)
      aer_dens = 2600.

    !-- Initialize dust source functions for MOSAIC-4bin aerosols
    ELSEIF(aer_is_mosaic4bin) THEN
      ! MOSAIC-4bin bin center diameters (micrometer)
      bin_diameters_cen(:) = (/ 0.0781, 0.3125, 1.2500, 5.0000 /)
      ! Size-resolved emission scaling factors for MOSAIC-4bin
      ! https://doi.org/10.1073/pnas.1014798108
      emiss_size_massfrac(:) = (/ 5.32E-05, 5.61E-03, 1.37E-01, 8.58E-01 /)
      aer_dens = 2600.

    ELSE ! aer_is_
      CALL wrf_error_fatal ('mineraldust_emis_init: chem_opt not supported')

    ENDIF ! aer_is_

  END SUBROUTINE mineraldust_emis_init

END MODULE module_mineraldust_emis

