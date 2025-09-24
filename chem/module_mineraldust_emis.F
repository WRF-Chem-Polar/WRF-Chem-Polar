! Louis Marelle 2025/03
!
! Purpose:
!  Compute size-resolved mineral dust aerosol emissions (mass and number) for the
!  sectional 4-bin or 8-bin MOSAIC aerosols, or for the GOCART 5 dust bins.
!
! References
! ==============================================================================
! Source functions:
!  (1) [Ginoux2004] Ginoux, P., Prospero, J. M., Torres, O., & Chin, M.
!      (2004). Long-term simulation of global dust distribution with the GOCART
!      model: correlation with North Atlantic Oscillation. Environmental Modelling
!      & Software, 19(2), 113-128 - doi:10.1016/S1364-8152(03)00114-2
!  (2) [MB1995] Marticorena, B., and G. Bergametti (1995), Modeling the
!      atmospheric dust cycle: 1. Design of a soil-derived dust emission
!      scheme, J. Geophys. Res., 100(D8), 16415–16430,
!      https://doi.org/10.1029/95JD00690
!  (3) [GZ2016] FLEXDUST - Groot Zwaaftink, C. D., Grythe, H., Skov, H., &
!      Stohl, A. (2016). Substantial contribution of northern high‐latitude
!      sources to mineral dust in the Arctic. Journal of Geophysical Research:
!      Atmospheres, 121(22), 13-678 - http://dx.doi.org/10.1002/2016JD025482
! Size distribution:
!  (1) [Kok2010] J.F. Kok, A scaling theory for the size distribution of
!      emitted dust aerosols suggests climate models underestimate the size of
!      the global dust cycle, Proc. Natl. Acad. Sci. U.S.A. 108 (3) 1016-1021,
!      https://doi.org/10.1073/pnas.1014798108 (2011).
! Soil moisture and precipitation correction factors
!  (1) [Fécan1999] Fécan, F., Marticorena, B., and Bergametti, G.: Parametrization
!      of the increase of the aeolian erosion threshold wind friction velocity
!      due to soil moisture for arid and semi-arid areas, Ann. Geophys., 17,
!      149–157, https://doi.org/10.1007/s00585-999-0149-7, 1999.
!  (2) [UNCCD2022] United Nations Convention to Combat Desertification
!      (UNCCD). 2022. Sand and Dust Storms Compendium: Information and
!      Guidance on Assessing and Addressing the Risks. Bonn, Germany.
!  (3) [GZ2017] Groot Zwaaftink, C. D., Arnalds, Ó.,
!      Dagsson-Waldhauserova, P., Eckhardt, S., Prospero, J. M., and
!      Stohl, A.: Temporal and spatial variability of Icelandic dust
!      emissions and atmospheric transport, Atmos. Chem. Phys., 17,
!      10865–10878, https://doi.org/10.5194/acp-17-10865-2017, 2017.
! General references
!  (1) [LeGrand2019]AFWA - LeGrand et al. (2019) - https://doi.org/10.5194/gmd-12-131-2019
!      (describes the original wrfgocart parameterizations and its issues)

MODULE module_mineraldust_emis

  IMPLICIT NONE

  !---- Global module-wide variables
  ! These are initialized in the init subroutine

  ! Allocatable size-resolved aerosol arrays (dim is size bin number and depends
  ! on aerosol mechanism)
  REAL, DIMENSION(:), ALLOCATABLE :: bin_diameters_cen
  REAL, DIMENSION(:), ALLOCATABLE :: emiss_size_massfrac

  ! Number of dust aerosol bins for chem_opt in use
  INTEGER :: naersizes

  ! Aerosol density (kg m-3)
  REAL :: aer_dens

  ! chem_opt aerosol model family
  LOGICAL :: aer_is_gocart
  LOGICAL :: aer_is_mosaic4bin
  LOGICAL :: aer_is_mosaic8bin

  ! Soil properties for calculating the saltation flux, from AFWA
  ! These values are identical to the data in phys/module_data_gocart_dust.F,
  ! but they have unclear names, and in case the GOCART module is changed in
  ! the future it is better to have a copy here
  INTEGER, PARAMETER :: NSOIL=9
  REAL, DIMENSION (NSOIL), PARAMETER :: REFF_SOIL=(/0.71E-6,1.37E-6,2.63E-6,5.00E-6,9.50E-6,18.1E-6,34.5E-6,65.5E-6,125.E-6/) ! m
  REAL, DIMENSION (NSOIL), PARAMETER :: DEN_SOIL=(/2500.,2650.,2650.,2650.,2650.,2650.,2650.,2650.,2650./) ! kg m-3
  REAL, DIMENSION (NSOIL), PARAMETER :: FRACSIZE_SOIL=(/1.,0.2,0.2,0.2,0.2,0.2,0.333,0.333,0.333/) ! fraction, 0-1
  INTEGER, DIMENSION (NSOIL), PARAMETER :: SPOINT_SOIL=(/1,2,2,2,2,2,3,3,3/)  ! 1 Clay, 2 Silt, 3 Sand


  CONTAINS

!------------------------------------------------------------------------------

  SUBROUTINE mineraldust_emis( dtstep, dz8w, tsk, &
         chem_opt, dust_opt, dust_hl_opt,   &
         rho_phy, chem, e_dust,                                     &
         u10, v10, ust, xland, xice, snowh, &
         raincv, rainncv, &
         clayfrac, sandfrac, &
         erod, erod_hl, vegfra, smois_top, isltyp, ivgtyp, &
         ids,ide, jds,jde, kds,kde,                                 &
         ims,ime, jms,jme, kms,kme,                                 &
         its,ite, jts,jte, kts,kte                                  )

    ! Purpose:
    !  Compute size-resolved mineral dust aerosol number and mass emissions

    ! Arguments:
    ! Input:
    !  dtstep = model time step (s)
    !  dz8w = 3D vertical level depth (m)
    !  tsk = skin temperature (K)
    !  chem_opt = chemical mechanism + aerosol mechanism option from config_flags
    !  dust_opt = dust emission option from config_flags
    !  dust_hl_opt = high-latitude dust option from config_flags
    !  rho_phy = 3D air density (kg m-3)
    !  u10 = windpeed at 10m in x (positive west-east) direction (m s-1)
    !  v10 = windpeed at 10m in y (positive south-north) direction (m s-1)
    !  ust = friction velocity (m s-1)
    !  xland = landmask flag (1 land, 2 water)
    !  xice = sea ice or lake ice fraction (0-1)
    !  snowh = snow height (m)
    !  raincv = time step rain (mm)
    !  rainncv = time step rain (mm)
    !  clayfrac = soil clay fraction (0-1)
    !  sandfrac = soil sand fraction (0-1)
    !  erod = Ginoux S-function scaled by bare soil fraction ("erodibility") (0-1)
    !  erod_hl = Ginoux S-function including high latitude sources ("erodibility") (0-1)
    !  vegfra = vegetation fractional cover in grid cell (0-1)
    !  smois_top = soil moisture in top soil layer (m3 m-3)
    !  isltyp = soil type category (int category)
    !  ivgtyp = land type category (int category)
    !  ids,ims,its etc. = domain, memory, and tile start and end indices
    ! Input/output
    !  chem = advected chemical species (gases in ppm, aerosols in µg kg-1)
    !  e_dust = accumulated dust emission flux (kg m-2)

    USE module_state_description, ONLY: PARAM_FIRST_SCALAR, NUM_CHEM, &
                           P_DUST_1, P_DUST_2, P_DUST_3, P_DUST_4, P_DUST_5, &
                           DUST_WRFGOCART, DUST_WCP, DUST_FLEXDUST
    USE module_model_constants, ONLY: PICONST
    USE module_data_gocart_dust, ONLY: POROSITY
    USE module_data_mosaic_asect, ONLY: numptr_aer, lptr_oin_aer, nsize_aer

    IMPLICIT NONE

    !---- Arguments
    ! Input
    INTEGER, INTENT(IN) :: chem_opt, dust_opt, dust_hl_opt
    INTEGER, INTENT(IN) :: ids,ide, jds,jde, kds,kde,               &
                           ims,ime, jms,jme, kms,kme,               &
                           its,ite, jts,jte, kts,kte
    INTEGER, DIMENSION(ims:ime, jms:jme), INTENT(IN) :: isltyp, ivgtyp
    REAL, INTENT(IN) :: dtstep
    REAL, DIMENSION(ims:ime, kms:kme, jms:jme), INTENT(IN) :: dz8w, rho_phy
    REAL, DIMENSION(ims:ime, jms:jme), INTENT(IN) :: tsk, u10, v10, ust, xland
    REAL, DIMENSION(ims:ime, jms:jme), INTENT(IN) :: xice, snowh, raincv, rainncv
    REAL, DIMENSION(ims:ime, jms:jme), INTENT(IN) :: smois_top, vegfra, clayfrac, sandfrac
    REAL, DIMENSION(ims:ime, jms:jme, 3), INTENT(IN) :: erod
    REAL, DIMENSION(ims:ime, jms:jme), INTENT(IN) :: erod_hl
    ! Input/output
    REAL, DIMENSION(ims:ime, kms:kme, jms:jme, NUM_CHEM), INTENT(INOUT) :: chem
    REAL, DIMENSION(ims:ime, jms:jme), INTENT(INOUT) :: e_dust

    !---- Local variables
    ! Dust emission fluxes (µg m-2 s-1, µm3 m-2 s-1, m-2 s-1)
    REAL :: dust_mass_emiss, dust_vol_emiss, dust_num_emiss
    ! Bulk dust mass emission flux (µg m-2 s-1)
    REAL :: dust_mass_emiss_tot
    ! 10-m wind speed (m s-1)
    REAL :: w10
    ! Friction velocity (m s-1)
    REAL :: ust_in
    ! Clay and sand fraction in soil (0-1)
    REAL :: clayfrac_in, sandfrac_in
    ! Skin temperature
    REAL :: tsk_in
    ! Conversion factor from m-2 s-1 to kg_dry_air-1
    REAL :: conv
    ! Surface wetness (0-1)
    REAL :: gwet
    ! Gravimetric soil water content in top layer (kg k-1g)
    REAL :: smois_grav
    ! Volumetric soil water content in top layer (m3 m-3)
    REAL :: smois_vol
    ! Aerosol diameter (m) and air density (kg m-3)
    REAL :: aer_diam, air_dens
    ! Soil grain diameter (m) and soil grain density (kg m-3)
    REAL :: soil_diam, soil_dens
    ! Fraction of soil type (sand, clay or silt) in cell (0-1)
    REAL :: soiltype_frac
    ! Soil freezing point (K)
    REAL :: tsk_freeze
    ! Precipitation rate (mm hr-1)
    REAL :: precip_rate
    ! Indices
    INTEGER :: i, j, isize, iphase
    INTEGER :: iaerbin, isoilbin
    INTEGER :: l_dust, l_num
    ! Ginoux S-function scaled by bare soil fraction ("Erodibility" in WRF-Chem)
    REAL :: s_ginoux
    ! Size distribution of dust in the original wrfgocart scheme
    REAL, DIMENSION(5) :: frac_s_wrfgocart = (/0.1,0.25,0.25,0.25,0.25/)

    !-------- Compute dust emissions --------
    ! No emissions to cloud-phase
    ! This should be = aiphase but it might not be initialized
    iphase = 1

    DO i=its,ite
      DO j=jts,jte

        ! Only perform dust emissions over snow-free land. Sea ice is
        ! considered land, so exclude sea ice points explicitly.
        IF (xland(i,j)<1.5 .AND. snowh(i,j)<0.1 .AND. xice(i,j)<=0.01) THEN
          IF(dust_hl_opt == 1) THEN
            ! Scale the S-function by (1.-vegfra): only bare soil can emit
            ! dust. This is not already done for HLdust EROD
            s_ginoux = erod_hl(i,j) * (100.-vegfra(i,j))/100.
          ELSE
            ! WPS EROD already includes bare soil scaling
            s_ginoux = SUM(erod(i,j,:))
          ENDIF

          ! Skip dust emission calculation if there is no potential dust source
          IF(s_ginoux > 1.E-3) THEN
            !-------- Initialize mineral dust emission parameters
            dust_mass_emiss = 0.0
            dust_num_emiss = 0.0
            w10 = SQRT(u10(i,j)**2.0 + v10(i,j)**2.0)
            ust_in = ust(i,j)
            clayfrac_in = clayfrac(i,j)
            sandfrac_in = sandfrac(i,j)
            air_dens = rho_phy(i,kts,j)
            precip_rate = (raincv(i,j)+rainncv(i,j)) / dtstep * 3600.
            tsk_in = tsk(i,j)
            ! Conversion factor from m-2 s-1 to kg_dry_air-1
            conv = 1.0/(rho_phy(i,kts,j) * dz8w(i,kts,j)) * dtstep
            ! Volumetric soil moisture over porosity
            smois_vol = smois_top(i,j)
            gwet = smois_vol/POROSITY(isltyp(i,j))
            gwet = MAX(1.E-3, MIN(1.0, gwet))
            ! Gravimetric soil water content
            smois_grav = smois_vol/((1.-POROSITY(isltyp(i,j)))*(2.65*(1-clayfrac_in)+2.50*clayfrac_in))

            !-------- Compute the total aerosol mass emission flux
            dust_mass_emiss_tot = 0.0
            SELECT CASE(dust_opt)
            CASE(DUST_WRFGOCART)
              !---- Compute dust mass emissions with the original WRF-Chem
              !     GOCART scheme (based on the 4.3.3 version)
              DO iaerbin = 1,naersizes
                aer_diam = bin_diameters_cen(iaerbin)*1.0E-6
                CALL source_mineraldust_gocart(dust_mass_emiss=dust_mass_emiss,  &
                       w10=w10, gwet=gwet, s_ginoux=s_ginoux, &
                       aer_dens=aer_dens, aer_diam=aer_diam, &
                       air_dens=air_dens)
                dust_mass_emiss_tot = dust_mass_emiss_tot &
                    + dust_mass_emiss*frac_s_wrfgocart(iaerbin)
              ENDDO

            CASE(DUST_WCP)
              !---- Compute dust mass emissions, WRF-Chem-Polar source
              ! Based on Ginoux2004, correcting original issues in the scheme
              ! following Legrand2019 and other corrections to enable
              ! emissions from high latitude regions
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
                CALL source_mineraldust_wcp(dust_mass_emiss=dust_mass_emiss, &
                       w10=w10, ust=ust_in, smois_grav=smois_grav, &
                       smois_vol=smois_vol, s_ginoux=s_ginoux, &
                       clayfrac=clayfrac_in, soil_dens=soil_dens, &
                       soil_diam=soil_diam, air_dens=air_dens, &
                       tsk=tsk_in, precip_rate=precip_rate)
                dust_mass_emiss_tot = dust_mass_emiss_tot &
                    + dust_mass_emiss*soiltype_frac*FRACSIZE_SOIL(isoilbin)
              ENDDO ! isoilbin

            CASE(DUST_FLEXDUST)
              !---- Compute dust mass emissions following Groot Zwaaftink et al. (2016)
              soil_dens = aer_dens
              CALL source_mineraldust_flexdust(dust_mass_emiss=dust_mass_emiss, &
                     w10=w10, ust=ust_in, smois_grav=smois_grav, &
                     s_ginoux=s_ginoux, clayfrac=clayfrac_in, &
                     sandfrac=sandfrac_in, soil_dens=soil_dens, &
                     air_dens=air_dens, precip_rate=precip_rate)
              dust_mass_emiss_tot = dust_mass_emiss


            CASE DEFAULT
              CALL wrf_error_fatal ('mineraldust_emis: dust_opt not supported')

            END SELECT

            !-------- Emit mineral dust aerosol mass and number for each size bin
            DO iaerbin = 1,naersizes
              !---- Distribute total dust emissions to size bin
              IF(dust_opt == DUST_WRFGOCART) THEN
                ! For WRFGOCART, use the old (incorrect) frac_s_wrfgocart to
                ! distribute emissions to each bin
                dust_mass_emiss = dust_mass_emiss_tot * frac_s_wrfgocart(iaerbin)
              ELSE
                ! For other dust emission models, distribute total emissions
                ! to each bin based on emiss_size_massfrac
                dust_mass_emiss = dust_mass_emiss_tot * emiss_size_massfrac(iaerbin)
              ENDIF
              e_dust(i,j) = e_dust(i,j) + dust_mass_emiss*dtstep*1.E-9

              !---- Convert dust aerosol mass emission flux to number
              ! emission flux (m-2 s-1)
              ! 1.E9 is for mm3 to um3 conversion
              dust_vol_emiss = dust_mass_emiss / aer_dens * 1.E9
              dust_num_emiss = dust_vol_emiss / ((PICONST/6.0)*bin_diameters_cen(iaerbin)**3.0)
              dust_num_emiss = MAX(dust_num_emiss, 0.0)

              !-------- Emit mineral dust emissions to chem(i,k,j,l)
              !---- Emit MOSAIC aerosols to chem
              IF(aer_is_mosaic4bin .OR. aer_is_mosaic8bin) THEN
                ! MOSAIC dust aerosol
                l_dust = lptr_oin_aer(iaerbin,1,iphase)
                IF (l_dust >= PARAM_FIRST_SCALAR) THEN
                  chem(i,kts,j,l_dust) = chem(i,kts,j,l_dust) &
                                         + dust_mass_emiss * conv
                END IF
                ! MOSAIC aerosol number
                l_num = numptr_aer(iaerbin,1,iphase)
                IF (l_num >= PARAM_FIRST_SCALAR) THEN
                  chem(i,kts,j,l_num) = chem(i,kts,j,l_num) + dust_num_emiss * conv
                END IF
              !---- Emit GOCART aerosols to chem
              ELSEIF(aer_is_gocart) THEN
                ! GOCART dust aerosol
                IF(iaerbin == 1 ) THEN
                  l_dust = P_DUST_1
                ELSEIF(iaerbin == 2 ) THEN
                  l_dust = P_DUST_2
                ELSEIF(iaerbin == 3 ) THEN
                  l_dust = P_DUST_3
                ELSEIF(iaerbin == 4 ) THEN
                  l_dust = P_DUST_4
                ELSEIF(iaerbin == 5 ) THEN
                  l_dust = P_DUST_5
                ENDIF
                IF (l_dust >= PARAM_FIRST_SCALAR) THEN
                  chem(i,kts,j,l_dust) = chem(i,kts,j,l_dust) &
                                         + dust_mass_emiss * conv
                ENDIF
              !---- Other chem_opt cases are not supported
              ELSE
                CALL wrf_error_fatal ('mineraldust_emis: chem_opt not supported')
              ENDIF ! aer_is_
            END DO ! iaerbin

          END IF ! s_ginoux
        END IF ! xland
      END DO ! j
    END DO ! i

  END SUBROUTINE mineraldust_emis

!------------------------------------------------------------------------------

  SUBROUTINE source_mineraldust_gocart(dust_mass_emiss, w10, gwet, s_ginoux, &
                                       aer_dens, aer_diam, air_dens)
    ! Purpose:
    !  Compute dust mass emission flux (µg m-2 s-1) as in the GOCART emissions scheme in WRF-Chem
    !
    ! Remarks:
    !  Only use this for GOCART mechanisms

    ! Arguments:
    ! Input:
    !  w10 = 10-m wind speed (m s-1)
    !  gwet = soil moisture (0-1)
    !  s_ginoux = Ginoux S-function ("erodibility") (0-1)
    !  aer_dens = aerosol density (kg m-3)
    !  aer_diam = aerosol diameter (m)
    !  air_dens = air density (kg m-3)
    ! Output:
    !  dust_mass_emiss = dust aerosol mass emission flux (µg m-2 s-1)

    USE module_model_constants, ONLY: G

    IMPLICIT NONE

    !---- Arguments
    ! Input
    REAL, INTENT(IN) :: w10, gwet, aer_dens, aer_diam, air_dens
    REAL, INTENT(IN) :: s_ginoux
    ! Output
    REAL, INTENT(OUT) :: dust_mass_emiss

    !---- Local variables
    ! Aerosol density, GOCART values (kg m-3)
    REAL :: aer_dens_gocart
    ! Threshold friction velocity for wind erosion (m s-1)
    REAL :: w_ts

    !---- Parameters
    ! Ginoux dust emissions scaling parameter, µg s-2 m-5
    ! 0.8 factor is ch_dust from chemics_init, 0.25 factor is from the
    ! erodibility scaling in the original GOCART emissions
    REAL, PARAMETER :: cdim_ginoux = 1.0 * 0.8 * 0.25

    ! Ginoux source functio
    IF(aer_diam < 2.0E-6) THEN
      aer_dens_gocart = 2500.
    ELSE
      aer_dens_gocart = 2650.
    ENDIF
    ! Threshold velocity (m s-1)
    w_ts = 0.13*1.0E-2*SQRT((aer_dens_gocart*1.E-3)*(G*1.E2)*(aer_diam*1.E2)/(air_dens*1.E-3))* &
           SQRT(1.0+0.006/(aer_dens_gocart*1.E-3)/(G*1.E2)/(aer_diam*1.E2)**2.5)/ &
           SQRT(1.928*(1331.0*(aer_diam*1.E2)**1.56+0.38)**0.092-1.0)
    ! Correct threshold velocity for soil moisture
    IF (gwet < 0.5) THEN
      w_ts = w_ts * (1.2 + 0.2*LOG10(gwet))
      w_ts = MAX(0.0, w_ts)
    ELSE
      w_ts = 1000.0
    ENDIF
    ! Calculate the dust emissions (µg m-2 s-1)
    IF(w10 > w_ts) THEN
      dust_mass_emiss = cdim_ginoux * s_ginoux * w10**2. * (w10 - w_ts)
    ELSE
      dust_mass_emiss = 0.0
    ENDIF
    dust_mass_emiss = MAX(0.0, dust_mass_emiss)

  END SUBROUTINE source_mineraldust_gocart

!------------------------------------------------------------------------------

  SUBROUTINE source_mineraldust_wcp(dust_mass_emiss, w10, ust, &
      smois_grav, smois_vol, s_ginoux, clayfrac, soil_dens, soil_diam, air_dens, &
      tsk, precip_rate)
    ! Purpose:
    !  Compute dust mass emission flux (µg m-2 s-1) using the source function
    !  developped for the WRF-Chem-Polar model, based on Ginoux2004, with bug
    !  corrections, updated parameters, and updates for emissions in high
    !  latitude regions

    ! References:
    !  (1) [Ginoux2004] Ginoux et al. (2004) - doi:10.1016/S1364-8152(03)00114-2
    !  (2) [MB1995] Marticorena and Bergametti (1995) - https://doi.org/10.1029/95JD00690
    !  (3) [Legrand2019] LeGrand et al. (2019) - https://doi.org/10.5194/gmd-12-131-2019
    ! Soil moisture, soil freezing, and precipitaton correction factors:
    !  (4) [Fécan1999] Fécan, F. et al. (1999) - https://doi.org/10.1007/s00585-999-0149-7
    !  (5) [UNCCD2022] UNCCD, 2022. Sand and Dust Storms Compendium: Information and
    !      Guidance on Assessing and Addressing the Risks. Bonn, Germany.
    !  (6) [GZ2017] Groot Zwaaftink C. D. et al. (2017) - https://doi.org/10.5194/acp-17-10865-2017

    ! Arguments:
    ! Input:
    !  w10 = 10-m wind speed (m s-1)
    !  ust = friction velocity (m s-1)
    !  smois_grav = gravimetric soil water content (kg kg-1)
    !  smois_vol = volumetric soil water content (m3 m-3)
    !  s_ginoux = Ginoux S-function ("erodibility") (0-1)
    !  clayfrac = soil clay fraction (0-1)
    !  soil_dens = soil grain density (kg m-3)
    !  soil_diam = soil grain diameter (m)
    !  air_dens = air density (kg m-3)
    !  tsk = skin temperature (K)
    !  precip_rate = instantaneous precipitation rate (mm hr-1)
    ! Output:
    !  dust_mass_emiss = dust aerosol mass emission flux (µg m-2 s-1)

    IMPLICIT NONE

    !---- Arguments
    ! Input
    REAL, INTENT(IN) :: w10, ust, smois_grav, smois_vol, clayfrac, soil_dens
    REAL, INTENT(IN) :: soil_diam, air_dens, s_ginoux, tsk, precip_rate
    ! Output
    REAL, INTENT(OUT) :: dust_mass_emiss

    !---- Local variables
    ! Threshold friction velocity for wind erosion (m s-1)
    REAL :: w_ts
    ! Threshold 10-m velocity for wind erosion (m s-1)
    REAL :: w10_ts
    ! Correction factor for the threshold velocity due to soil moisture
    ! w_ts_wet/w_ts_dry (unitless)
    REAL :: fcorr_w_ts
    ! Soil moisture correction factor (0-1)
    REAL :: soil_moisture_corr
    ! Soil freezing point (K)
    REAL :: tsk_freeze

    !---- Parameters
    ! Ginoux dust emissions scaling parameter, µg s-2 m-5
    ! 0.25*0.8 = 0.2 factor from original GOCART implementation
    REAL, PARAMETER :: CDIM = 0.2
    ! Use soil moisture to limit emissions, by increasing the threshold
    ! velocity and decreasing the total flux. If false, emissions are not
    ! affected by soil moisture directly, but emissions are cut when
    ! precipitation exceeds 1.0 mm hr-1 following GZ2017
    LOGICAL, PARAMETER :: USE_SOIL_MOISTURE = .true.
    REAL, PARAMETER :: PRECIP_RATE_CRIT = 1.0 ! mm hr-1

    dust_mass_emiss = 0.0

    ! Calculate soil freezing point, from UNCCD2022
    IF(smois_vol < 0.2) THEN
      tsk_freeze = 260.
    ELSEIF(smois_vol < 0.4) THEN
      tsk_freeze = 260. + 12.*(smois_vol-0.2)/0.2
    ELSE
      tsk_freeze = 272.
    ENDIF

    ! No emissions over frozen soils
    IF(tsk > tsk_freeze) THEN
      ! Threshold velocity from MB1995
      w_ts = w_ts_mb1995(soil_diam=soil_diam, soil_dens=soil_dens, &
                         air_dens=air_dens)

      ! Soil moisture / precipitation inhibition of emissions
      IF(USE_SOIL_MOISTURE) THEN
        ! Correct threshold velocity for soil moisture
        ! We use Fécan1999 instead of the original Ginoux2004 factor
        fcorr_w_ts = fcorr_f1999(smois_grav=smois_grav, clayfrac=clayfrac)
        w_ts = w_ts * fcorr_w_ts
      ELSEIF(precip_rate > PRECIP_RATE_CRIT) THEN
        ! Disable emissions if it rains (GZ2017)
        w_ts = 1000.0
      ENDIF ! USE_SOIL_MOISTURE

      ! Calculate the dust emissions (µg m-2 s-1)
      IF(ust > w_ts .AND. ust > 1.E-6) THEN
        ! Convert threshold friction velocity to a threshold 10-m velocity using
        ! the local w10/ust ratio
        w10_ts = w_ts * w10/ust
        IF(w10 > w10_ts) THEN
          dust_mass_emiss = CDIM * s_ginoux &
                            * w10**2. * (w10 - w10_ts)
        ENDIF
      ENDIF ! ust > w_ts

      ! Soil moisture reduction of emissions
      IF(USE_SOIL_MOISTURE) THEN
        ! Correct emission flux for soil moisture, from UNCCD2022
        IF(smois_vol > 0.15) THEN
          soil_moisture_corr = 0.15/smois_vol
        ELSE
          soil_moisture_corr = 1.0
        ENDIF
      ENDIF ! USE_SOIL_MOISTURE

      dust_mass_emiss = dust_mass_emiss * soil_moisture_corr
      ! Prevent negative emissions
      dust_mass_emiss = MAX(0.0, dust_mass_emiss)
    ENDIF ! tsk

  END SUBROUTINE source_mineraldust_wcp

!------------------------------------------------------------------------------

  SUBROUTINE source_mineraldust_flexdust(dust_mass_emiss, w10, ust, &
               smois_grav, s_ginoux, clayfrac, sandfrac, soil_dens, &
               air_dens, precip_rate)
    ! Purpose:
    !  Compute dust mass emission flux (µg m-2 s-1) as in the FLEXDUST model
    !  Does not include special treatment for the Icelandic emissions as in
    !  the original FLEXDUST model, but the s_ginoux factor can be updated to
    !  increase emissions from Iceland.

    ! References:
    !  (1) [GZ2016] Groot Zwaaftink et al. (2016) - http://dx.doi.org/10.1002/2016JD025482
    !  (2) [MB1995] Marticorena and Bergametti (1995) - https://doi.org/10.1029/95JD00690
    ! Soil moisture and precipitaton correction factors
    !  (1) [Fécan1999] Fécan et al. (1999) - https://doi.org/10.1007/s00585-999-0149-7
    !  (2) [GZ2017] Groot Zwaaftink C. D. et al. (2017) - https://doi.org/10.5194/acp-17-10865-2017

    ! Arguments:
    ! Input:
    !  w10 = 10-m wind speed (m s-1)
    !  ust = friction velocity (m s-1)
    !  smois_grav = gravimetric soil water content (kg kg-1)
    !  s_ginoux = Ginoux s-function ("erodibility") (0-1)
    !  clayfrac = soil clay fraction (0-1)
    !  sandfrac = soil sand fraction (0-1)
    !  soil_dens = soil grain density (kg m-3)
    !  air_dens = air density (kg m-3)
    !  precip_rate = instantaneous precipitation rate (mm hr-1)
    ! Output:
    !  dust_mass_emiss = dust aerosol mass emission flux (µg m-2 s-1)

    USE module_model_constants, ONLY: G

    IMPLICIT NONE

    !---- Arguments
    ! Input
    REAL, INTENT(IN) :: w10, ust, smois_grav, soil_dens, air_dens, precip_rate
    REAL, INTENT(IN) :: clayfrac, sandfrac
    REAL, INTENT(IN) :: s_ginoux
    ! Output
    REAL, INTENT(OUT) :: dust_mass_emiss

    !---- Local variables
    ! Threshold velocity for wind erosion (m s-1)
    REAL :: w_ts
    ! Soil particle diameter (m)
    REAL :: soil_diam
    ! Sandblasting efficiency for emissions (m-1)
    REAL :: alpha_sandblast
    ! Adjusted clay and sand fraction (%)
    REAL :: clayfrac_adj, sandfrac_adj
    ! Correction factor for the threshold velocity due to soil moisture
    ! w_ts_wet/w_ts_dry (unitless)
    REAL :: fcorr_w_ts

    !---- Parameters
    REAL, PARAMETER :: sandfrac_minval = 1. ! %
    REAL, PARAMETER :: scaling_factor_mb95 = 4.8E-4
    ! Use soil moisture to limit emissions, by increasing the threshold
    ! velocity. If false, emissions are not affected by soil moisture
    ! directly, but emissions are cut when precipitation exceeds 1.0 mm hr-1
    ! following GZ2017
    LOGICAL, PARAMETER :: USE_SOIL_MOISTURE = .true.
    REAL, PARAMETER :: PRECIP_RATE_CRIT = 1.0 ! mm hr-1

    !---- Calculate dust emission flux following GZ2016
    sandfrac_adj = sandfrac*100.
    ! Threshold velocity (GZ2016 eq. 2)
    IF(sandfrac_adj >= sandfrac_minval) THEN
      soil_diam = 75.E-6
    ELSE
      soil_diam = 10.E-6
    ENDIF
    w_ts = SQRT( 0.0123 * (soil_dens/air_dens * G * soil_diam &
                          + 3E-4/(air_dens*soil_diam)) )
    ! Soil moisture / precipitation inhibition of emissions
    IF(USE_SOIL_MOISTURE) THEN
      ! Correct threshold velocity for soil moisture (Fécan et al., 1999)
      fcorr_w_ts = fcorr_f1999(smois_grav=smois_grav, clayfrac=clayfrac)
      w_ts = w_ts * fcorr_w_ts
    ELSEIF(precip_rate > PRECIP_RATE_CRIT) THEN
      ! Disable emissions if it rains (GZ2017)
      w_ts = 1000.0
    ENDIF ! USE_SOIL_MOISTURE
    ! Sandblasting efficiency, from MB95 (m-1)
    clayfrac_adj = MIN(clayfrac*100., 20.)
    alpha_sandblast = 100. * 10.**(0.134*clayfrac_adj-6.)
    IF(ust > w_ts) THEN
      ! Dust emissions following GZ2016 based on MB95
      dust_mass_emiss = 1.E9 * s_ginoux * scaling_factor_mb95 &
                       * alpha_sandblast * air_dens * ust**3. / G &
                       * (1. - w_ts**2./ust**2.) * (1. + w_ts/ust)
    ELSE
      dust_mass_emiss = 0.0
    ENDIF
    dust_mass_emiss = MAX(0.0, dust_mass_emiss)

  END SUBROUTINE source_mineraldust_flexdust


!------------------------------------------------------------------------------

  FUNCTION w_ts_mb1995(soil_diam, soil_dens, air_dens)
    ! Threshold friction velocity from Marticorena and Bergametti, 1995
    !
    ! References:
    !   [MB1995] Marticorena, B., and G. Bergametti (1995), Modeling the
    !   atmospheric dust cycle: 1. Design of a soil-derived dust emission
    !   scheme, J. Geophys. Res., 100(D8), 16415–16430,
    !   https://doi.org/10.1029/95JD00690
    !
    ! Arguments:
    !   soil_diam = soil grain diameter (m)
    !   soil_dens = soil grain density (kg m-3)
    !   air_dens = air density (kg m-3)
    ! Output:
    !   w_ts_mb1995 = Threshold friction velocity for dust (m s-1)
    !
    USE module_model_constants, ONLY: G
    IMPLICIT NONE
    !-- Arguments
    REAL, INTENT(IN) :: soil_diam, soil_dens, air_dens
    REAL :: w_ts_mb1995
    !-- Local variables
    ! K factor from MB1995 eq. 4
    REAL :: k_mb1995
    !-- Calculate the threshold friction velocity
    ! MB1995 eq. 4, soil_diam in m instead of cm, and dens in kg m-3
    ! instead of g cm-3
    k_mb1995 = SQRT(soil_dens/air_dens * G * soil_diam) &
               * SQRT(1. + 6.E-7/(soil_dens*G*(soil_diam)**2.5))
    ! MB1995 eq. 6 and 5., soil_diam in m instead of cm
    w_ts_mb1995 = 1.E-2 * 0.129 * k_mb1995 &
           / SQRT(1.928 * (1754600.*(soil_diam)**1.56 + 0.38)**0.092 - 1.)
    RETURN
  END FUNCTION w_ts_mb1995

!------------------------------------------------------------------------------

  FUNCTION fcorr_f1999(smois_grav, clayfrac)
    ! Correction factor for the threshold friction velocity due to soil
    ! moisture effects, from Fécan et al. (1999)
    !
    ! References:
    !   [Fécan1999] Fécan, F., Marticorena, B., and Bergametti, G.:
    !   Parametrization of the increase of the aeolian erosion threshold wind
    !   friction velocity due to soil moisture for arid and semi-arid areas,
    !   Ann. Geophys., 17, 149–157, https://doi.org/10.1007/s00585-999-0149-7,
    !   1999.
    !
    ! Arguments:
    !   smois_grav = gravimetric surface soil water content (kg k-1g)
    !   clayfrac = clay fraction in soil (0-1)
    ! Output:
    !   fcorr_w_ts_f1999 = Correction factor w_ts_wet/w_ts_dry (unitless)
    !
    IMPLICIT NONE
    !-- Arguments
    REAL, INTENT(IN) :: smois_grav, clayfrac
    REAL :: fcorr_f1999
    !-- Local variables
    ! Clay fraction in soil (%)
    REAL :: clayfrac_adj
    ! Threshold gravimetric soil water content (% mass)
    REAL :: smois_thresh
    !-- Calculate the correction factor (Fécan1999 eq. 15)
    clayfrac_adj = MIN(clayfrac, 0.5)
    smois_thresh = 17.*clayfrac_adj + 14.*(clayfrac_adj)**2.0
    IF(smois_grav*100. > smois_thresh) THEN
      fcorr_f1999 = SQRT(1. + 1.21*(smois_grav*100. - smois_thresh)**0.68)
    ENDIF
    RETURN
  END FUNCTION fcorr_f1999

!------------------------------------------------------------------------------

  SUBROUTINE mineraldust_emis_init(chem_opt)

    ! Purpose:
    !  Initialize the parameters for the dust aerosol source functions used
    !  in mineraldust_emis, which depend on chem_opt
    !
    ! References:
    !  [Kok2010] J.F. Kok, A scaling theory for the size distribution of
    !  emitted dust aerosols suggests climate models underestimate the size of
    !  the global dust cycle, Proc. Natl. Acad. Sci. U.S.A. 108 (3) 1016-1021,
    !  https://doi.org/10.1073/pnas.1014798108 (2011).

    USE module_state_description, ONLY: MOZCART_KPP, T1_MOZCART_KPP, &
          GOCART_SIMPLE, GOCARTRADM2, GOCARTRACM_KPP, &
          CBMZ_MOSAIC_4BIN, CBMZ_MOSAIC_DMS_4BIN, CBMZ_MOSAIC_4BIN_AQ, &
          CBMZ_MOSAIC_DMS_4BIN_AQ, MOZART_MOSAIC_4BIN_KPP, &
          MOZART_MOSAIC_4BIN_AQ_KPP, MOZART_MOSAIC_4BIN_DMS_AQ_KPP, &
          CRI_MOSAIC_4BIN_AQ_KPP, SAPRC99_MOSAIC_4BIN_VBS2_KPP, CBMZ_MOSAIC_8BIN, &
          CBMZ_MOSAIC_DMS_8BIN, CBMZ_MOSAIC_8BIN_AQ, CBMZ_MOSAIC_DMS_8BIN_AQ, &
          CRI_MOSAIC_8BIN_AQ_KPP, SAPRC99_MOSAIC_8BIN_VBS2_AQ_KPP, &
          SAPRC99_MOSAIC_8BIN_VBS2_KPP

    IMPLICIT NONE

    ! Arguments:
    !  Input
    INTEGER, INTENT(IN) :: chem_opt

    ! Parameters:
    ! Flag for the size distribution of dust emissions. If true, include the
    ! dust mass beyond the aerosol scheme size limits into the first and last
    ! bins
    ! LOGICAL, PARAMETER :: INCLUDE_ALL_DUST_MASS = .true.
    LOGICAL, PARAMETER :: INCLUDE_ALL_DUST_MASS = .false.

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
    ! Initialize emiss_size_massfrac (µm3 m-2 s-1) and bin_diameters_cen
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
      ! Kok2010, eq. 6
      IF(INCLUDE_ALL_DUST_MASS) THEN
        emiss_size_massfrac(:) = (/ 0.0435, 0.106, 0.219, 0.485, 0.1465 /)
      ELSE
        emiss_size_massfrac(:) = (/ 0.0435, 0.106, 0.219, 0.485, 0.144 /)
      ENDIF
      aer_dens = 2650.

    !-- Initialize dust source functions for MOSAIC-8bin aerosols
    ELSEIF(aer_is_mosaic8bin) THEN
      ! MOSAIC-8bin bin center diameters (micrometer)
      bin_diameters_cen(:) = (/ 0.0552, 0.1105, 0.2210, 0.4419, &
                                0.8839, 1.7678, 3.5355, 7.0711 /)
      ! Size-resolved emission scaling factors for MOSAIC-8bin
      ! Kok2010, eq. 6
      IF(INCLUDE_ALL_DUST_MASS) THEN
        emiss_size_massfrac(:) = (/ 8.12E-07, 1.48E-05, 1.95E-04, 1.80E-03, &
                                    1.19E-02, 5.73E-02, 2.02E-01, 7.27E-01 /)
      ELSE
        emiss_size_massfrac(:) = (/ 7.83E-07, 1.48E-05, 1.95E-04, 1.80E-03, &
                                    1.19E-02, 5.73E-02, 2.02E-01, 4.56E-01 /)
      ENDIF
      aer_dens = 2600.

    !-- Initialize dust source functions for MOSAIC-4bin aerosols
    ELSEIF(aer_is_mosaic4bin) THEN
      ! MOSAIC-4bin bin center diameters (micrometer)
      bin_diameters_cen(:) = (/ 0.0781, 0.3125, 1.2500, 5.0000 /)
      ! Size-resolved emission scaling factors for MOSAIC-4bin
      ! Kok2010, eq. 6
      IF(INCLUDE_ALL_DUST_MASS) THEN
        emiss_size_massfrac(:) = (/ 1.56E-05, 2.00E-03, 6.92E-02, 9.288E-01 /)
      ELSE
        emiss_size_massfrac(:) = (/ 1.56E-05, 2.00E-03, 6.92E-02, 6.57E-01  /)
      ENDIF
      aer_dens = 2600.

    ELSE ! aer_is_
      CALL wrf_error_fatal ('mineraldust_emis_init: chem_opt not supported')

    ENDIF ! aer_is_

  END SUBROUTINE mineraldust_emis_init

END MODULE module_mineraldust_emis

