! 2025/11 Louis Marelle
!
! Purpose:
! - Perform wet scavenging of aerosols by grid scale rainout/washout
!
! Remarks:
! - Washout and rainout from Luo et al. (2019, 2020)
!
! References:
! ==============================================================================
! (1) Luo, G., Yu, F., and Moch, J. M.: Further improvement of wet process
!     treatments in GEOS-Chem v12.6.0: impact on global distributions of
!     aerosols and aerosol precursors, Geosci. Model Dev., 13, 2879–2903,
!     https://doi.org/10.5194/gmd-13-2879-2020, 2020.
! (2) Luo, G., Yu, F., and Schwab, J.: Revised treatment of wet scavenging
!     processes dramatically improves GEOS-Chem 12.0.0 simulations of surface
!     nitric acid, nitrate, and ammonium over the United States, Geosci. Model
!     Dev., 12, 3439–3447, https://doi.org/10.5194/gmd-12-3439-2019, 2019.


MODULE module_aer_wetscav_simple

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: aer_wetscav_simple, aer_wetscav_simple_init

  !---- Global module-wide variables
  !-- Allocatable aerosol arrays
  ! Dim is number of aerosol species and depends on aerosol mechanism.
  ! These are initialized in aer_wetscav_simple_init
  ! Number of wet-scavenged aerosol species for chem_opt (dim of
  ! aer_pointers)
  INTEGER :: n_aer
  ! List of aerosol species experiencing rainout/washout
  INTEGER, DIMENSION(:), ALLOCATABLE :: aer_pointers
  ! Rainout activation efficiency of species (0-1)
  REAL, DIMENSION(:), ALLOCATABLE :: aq_aer_ratio, ice_aer_ratio
  ! Aerosol diameter (µm) for each species in aer_pointers
  REAL, DIMENSION(:), ALLOCATABLE :: aer_dry_diam
  ! Deposition flux of species during timestep (ug/m2)
  REAL, DIMENSION(:), ALLOCATABLE :: deposition_flux

  !---- Parameters
  ! Minimum values for performing wet scavenging:
  ! Aerosol MMR (µg/kg), cloud/precip water MMR (µg/kg), cloud/precip fraction
  ! of cell area (0-1), loss fraction from wet scavenging process (0-1)
  REAL, PARAMETER :: AER_MINVAL = 1.0E-12, WATER_MINVAL = 1.0E-12, &
                     AERA_FRAC_MINVAL = 0.0001, LOSSFRAC_MINVAL = 1.0E-6
  ! Scaling factor for resuspension rate, to take into account that
  ! evaporation of hydrometeors can be incomplete. Value 0.5 from GEOS-CHEM
  REAL, PARAMETER :: RESUSP_SCALING = 0.5
  ! Flags for disabling in-cloud or below-cloud scavenging, or resuspension
  LOGICAL, PARAMETER :: DO_RAINOUT = .true., DO_WASHOUT = .true., &
                        DO_RESUSPENSION = .true.

  ! chem_opt aerosol model family
  LOGICAL :: aer_is_gocart

CONTAINS

!----------------------------------------------------------------------------

  SUBROUTINE aer_wetscav_simple( chem, qsrflx, dtstep, config_flags,  &
                             t_phy, rho_phy, dz8w,                        &
                             cldfra, rainprod, evapprod,                  &
                             precr, preci, precs, precg,                  &
                             qc, qi, qs, qr, qg,                          &
                             ids,ide, jds,jde, kds,kde,                   &
                             ims,ime, jms,jme, kms,kme,                   &
                             its,ite, jts,jte, kts,kte                    )

  ! Purpose:
  ! Calculates and removes wet scavenging tendencies for aerosols
  !
  ! Arguments:
  ! Input:
  !   dtstep = model time step (s)
  !   config_flags = namelist configuration flags
  !   t_phy = 3D physical temperature (K)
  !   dz8w = 3D vertical level depth (m)
  !   rho_phy = 3D air density (kg/m3)
  !   cldfra = 3D cloud fraction (unitless, 0-1)
  !   rainprod = 3D rain production rate (kg/kg dry air/s)
  !   evapprod = 3D rain evaporation rate (kg/kg dry air/s)
  !   precr = 3D rain precipitation rate through grid cell (kg/m2/s)
  !   preci = 3D ice precipitation rate through grid cell (kg/m2/s)
  !   precs = 3D snow precipitation rate through grid cell (kg/m2/s)
  !   precg = 3D graupel precipitation rate through grid cell (kg/m2/s)
  !   qc = cloud water mixing ratio (kg/kg dry air)
  !   qi = cloud ice mixing ratio (kg/kg dry air)
  !   qs = snow water mixing ratio (kg/kg dry air)
  !   qr = rain water mixing ratio (kg/kg dry air)
  !   qg = graupel water mixing ratio (kg/kg dry air)
  !   ids,ims,its etc. = domain, memory, and tile start and end indices
  ! Input/output:
  !   chem = advected chemical species (gases in ppm, aerosols in ug/kg)
  ! Output:
  !   qsrflx = aerosol column change due to scavenging (µg/m2)

  ! USE associations:
  USE module_configure, ONLY: param_first_scalar, num_chem, grid_config_rec_type
  USE module_state_description, ONLY: &
        p_seas_1, p_seas_2, p_seas_3, p_seas_4, p_dust_1, &
        p_dust_2, p_dust_3, p_dust_4, p_dust_5

  IMPLICIT NONE

  !---- Arguments
  TYPE(grid_config_rec_type), INTENT(IN) :: config_flags
  INTEGER, INTENT(IN) :: ids,ide, jds,jde, kds,kde,    &
                         ims,ime, jms,jme, kms,kme,    &
                         its,ite, jts,jte, kts,kte
  REAL, INTENT(IN) :: dtstep
  REAL, DIMENSION(ims:ime,kms:kme,jms:jme), INTENT(IN) :: t_phy, rho_phy, &
      dz8w, qi, qc, qs, qr, qg, rainprod, evapprod, precr, preci, precs, &
      precg, cldfra
  REAL, DIMENSION(ims:ime,kms:kme,jms:jme,num_chem), INTENT(INOUT) :: chem
  REAL, DIMENSION(ims:ime, jms:jme, num_chem),          &
        INTENT(OUT) :: qsrflx

  !---- Local variables
  ! x,y,z, and chemical species index, chemical species pointer
  INTEGER :: i, j, k, iaer, p_chem
  ! Cell fractions where rainout/washout is applied (See Luo et al., 2020)
  REAL :: frac_rainout, frac_washout
  ! Vertical dimension of cldfra for cloud fraction calculation
  INTEGER :: dimz
  ! Aerosol conc. in cloud phase, aerosol conc. lost to rainout, aerosols
  ! conc. resuspended from precip (ug/kg)
  REAL :: cloudborne_aer, rainout_aer, resusp_aer
  ! Fraction of aerosols activated in clouds (0-1)
  ! Total fraction, liquid clouds only, ice clouds only
  REAL :: cloudborne_frac, cloudborne_frac_liq, cloudborne_frac_ice
  ! Cloud water mixing ratio (kg/kg dry air), precipitation mixing ratio
  ! (kg/kg dry air), liquid fraction of cloud water (0-1)
  REAL :: cloud_water, precip_water, liquid_fraction
  ! Rainout activation efficiency of species (0-1)
  REAL :: aq_aer_ratio_iaer, ice_aer_ratio_iaer
  ! Total 3D precipitation trough cell for rain/snow/ice (kg/m2/hr), rain
  ! production rate at i,k,j (kg/kg/s)
  REAL :: ppr, pps, ppi, rainrate
  ! Temperature (K)
  REAL :: tk
  ! Cloud fraction at i,k,j (0D)
  REAL :: cloud_fraction
  ! Fraction of aerosol species washed out/rained out/resuspended (0-1)
  REAL :: washfrac, rainfrac, resusp_frac
  ! Efficiency of aerosol uptake in cloud water and cloud ice (0-1), for GOCART
  REAL :: cloud_uptk_effic_liq, cloud_uptk_effic_ice
  ! Conversion factor from /kg to /m2
  REAL :: conv
  ! 1D vertical arrays for rainfrac, wahsfrac, resusp_frac, cloud_aer_effic values
  REAL, DIMENSION(kts:kte) :: rainfrac_1d, washfrac_fine_1d, washfrac_coarse_1d, resusp_frac_1d
  REAL, DIMENSION(kts:kte) :: cloud_uptk_effic_liq_1d, cloud_uptk_effic_ice_1d
  ! First-order rainout rate (s-1)
  REAL :: k_rain
  ! Value of chem(i,k,j,p) before deposition step (ug/kg)
  REAL :: chem_spec_before

  !-------- Perform wet scavenging for aerosol species -------
  IF(DO_RAINOUT .OR. DO_WASHOUT) THEN
  DO j=jts,jte
  DO i=its,ite
    !-- Initialize rainout, washout, resuspension parameters
    frac_washout = 0.0
    rainfrac_1d(kts:kte) = 0.0
    washfrac_fine_1d(kts:kte) = 0.0
    washfrac_coarse_1d(kts:kte) = 0.0
    resusp_frac_1d(kts:kte) = 0.0
    cloud_uptk_effic_liq_1d(kts:kte) = 1.0
    cloud_uptk_effic_ice_1d(kts:kte) = 1.0

    !---- Calculate parameters for rainout and washout schemes
    DO k=kte-1,kts,-1 ! Do not change, needs to loop in descending order

      !-- Initialize physical parameters at i,k,j
      tk = t_phy(i, k, j)
      cloud_fraction = cldfra(i,k,j)
      ! Add rainprod back to cloud_water because rainprod was already
      ! removed from cloud water in microphysics
      cloud_water = qc(i,k,j) + qi(i,k,j) + rainprod(i,k,j)*dtstep
      ! Add evapprod back to precip_water because evapprod was already
      ! removed from precip in microphysics
      precip_water = qr(i,k,j) + qs(i,k,j) + qg(i,k,j) &
                     + evapprod(i,k,j)*dtstep
      IF(cloud_water > WATER_MINVAL) THEN
        liquid_fraction = qc(i,k,j) / cloud_water
        liquid_fraction = MIN(1.0, MAX(0.0, liquid_fraction))
      ELSE
        liquid_fraction = 1.0
      ENDIF
      rainrate = rainprod(i,k,j)
      ! Total precipitation through grid cell in kg/m2/hr
      ppr = MAX(0., precr(i,k,j) * 3600.)
      pps = MAX(0., precs(i,k,j) * 3600.)
      ppi = MAX(0., (preci(i,k,j) + precg(i,k,j)) * 3600.)

      !-- Calculate rainout parameters
      ! frac_rainout is also needed to calculate frac_washout for washout
      CALL calc_rainout_parameters( rainfrac, frac_rainout, &
               cloud_uptk_effic_liq, cloud_uptk_effic_ice, &
               tk, cloud_fraction, cloud_water, liquid_fraction, &
               rainrate, dtstep)
      rainfrac_1d(k) = rainfrac
      cloud_uptk_effic_liq_1d(k) = cloud_uptk_effic_liq
      cloud_uptk_effic_ice_1d(k) = cloud_uptk_effic_ice

      !-- Calculate washout parameters
      ! The grid cell fraction experiencing washout, frac_washout, is the
      ! max between frac_rainout at k and frac_washout at k+1 (this is why
      ! we loop over k in descending order)
      frac_washout = MIN(1.0, MAX(frac_washout, frac_rainout))
      IF(DO_WASHOUT .AND. ppr+pps+ppi > 0.0 &
         .AND. frac_washout > AERA_FRAC_MINVAL) THEN
        CALL washfrac_fine_aerosol(washfrac, ppr, pps, ppi, tk, &
                                   frac_washout, dtstep)
        washfrac_fine_1d(k) = washfrac
        CALL washfrac_coarse_aerosol(washfrac, ppr, pps, ppi, tk, &
                                     frac_washout, dtstep)
        washfrac_coarse_1d(k) = washfrac
      ENDIF ! DO_WASHOUT .AND. ...

      !-- Calculate resuspension parameters
      ! We apply resuspension to the total deposition flux from the column
      ! above, so this assumes that the timescale of precipitation to level
      ! k is smaller than dtstep
      IF(DO_RESUSPENSION .AND. precip_water > WATER_MINVAL &
         .AND. evapprod(i,k,j) > 0.) THEN
        resusp_frac = RESUSP_SCALING * evapprod(i,k,j) * dtstep / precip_water
        resusp_frac = MAX(0.0, MIN(1.0, resusp_frac))
        resusp_frac_1d(k) = resusp_frac
      END IF ! DO_RESUSPENSION .AND. ...

    END DO ! k

    !---- Perform rainout and washout
    deposition_flux(:) = 0.0
    DO iaer = 1,n_aer
      ! Calculate rainout and washout by precip falling through the
      ! vertical column
      ! Do not change the loop order, needs to loop in descending order
      DO k=kte-1,kts,-1
        conv = rho_phy(i,k,j) * dz8w(i,k,j)
        IF(aer_is_gocart) THEN
          p_chem = aer_pointers(iaer)
        ENDIF

        !---- Perform rainout ("in-cloud scavenging")
        IF(DO_RAINOUT .AND. rainfrac_1d(k) > LOSSFRAC_MINVAL) THEN
          !-- Rainout: calculate the mass concentration of cloudborne aerosols
          IF(aer_is_gocart) THEN
            !- GOCART: Liquid-phase cloudborne aerosols for rainout
            aq_aer_ratio_iaer = aq_aer_ratio(iaer)
            cloudborne_frac_liq = aq_aer_ratio_iaer*cloud_uptk_effic_liq_1d(k)
            !- GOCART: Ice-phase cloudborne aerosols for rainout
            IF(tk <= 237.0) THEN
              ! Homogeneous nucleation - In cold clouds, all aerosols can
              ! activate, but reduce the max activation to 50% as in Luo
              ! (2019, 2010)
              ice_aer_ratio_iaer = 0.5
            ELSE
              ! Heterogeneous nucleation - In mixed-phase clouds, ice
              ! activation is only possible for INP
              ice_aer_ratio_iaer = ice_aer_ratio(iaer)
            ENDIF
            cloudborne_frac_ice = ice_aer_ratio_iaer*cloud_uptk_effic_ice_1d(k)
          ELSE
            CALL wrf_error_fatal('aer_wetscav_simple: unknown chem_opt')
          ENDIF ! aer_is_gocart

          !-- Appply rainout to chem
          IF(cloudborne_frac_liq > 0.0 .OR. cloudborne_frac_ice > 0.) THEN
            !- GOCART: apply rainout
            IF(aer_is_gocart) THEN
              ! Remove GOCART mass
              IF(p_chem > param_first_scalar) THEN
                cloudborne_frac = cloudborne_frac_liq+cloudborne_frac_ice
                CALL remove_aer(aerconc=chem(i,k,j,p_chem), &
                                dep_flux=deposition_flux(iaer), &
                                lossfrac=rainfrac_1d(k)*cloudborne_frac, &
                                conv=conv)
              ENDIF
            ELSE
              CALL wrf_error_fatal('aer_wetscav_simple: unknown chem_opt')
            ENDIF ! aer_is_gocart
          ENDIF

        ENDIF ! DO_RAINOUT .AND. ...


        !---- Perform washout ("below-cloud scavenging")
        IF(DO_WASHOUT .AND. (washfrac_fine_1d(k) > LOSSFRAC_MINVAL &
                        .OR. washfrac_coarse_1d(k) > LOSSFRAC_MINVAL)) THEN
          !-- Calculate washfrac, the fraction of washed out aerosols
          ! Fine aerosol washout (diam<=2 micrometer)
          IF(aer_dry_diam(iaer) <= 2.0E-6) THEN
            washfrac = washfrac_fine_1d(k)
          ! Coarse aerosol washout (diam>2 micrometer)
          ELSEIF(aer_dry_diam(iaer) > 2.0E-6) THEN
            washfrac = washfrac_coarse_1d(k)
          ELSE
            washfrac = 0.0
          END IF

          !-- Apply washout to chem
          IF(washfrac > 0.0) THEN
            !- GOCART: apply washout
            IF(aer_is_gocart) THEN
              ! Remove GOCART mass
              IF(p_chem > param_first_scalar) THEN
                CALL remove_aer(aerconc=chem(i,k,j,p_chem), &
                                dep_flux=deposition_flux(iaer), &
                                lossfrac=washfrac, conv=conv)
              ENDIF
            ELSE
              CALL wrf_error_fatal('aer_wetscav_simple: unknown chem_opt')
            ENDIF
          ENDIF

        END IF ! DO_WASHOUT .AND. ...


        !---- Perform resuspension of aerosols from evaporating precip. water
        ! We apply resuspension to the total deposition flux from the column
        ! above
        IF(DO_RESUSPENSION .AND. resusp_frac_1d(k) > LOSSFRAC_MINVAL) THEN
          !-- GOCART: apply resuspension
          IF(aer_is_gocart .AND. p_chem > param_first_scalar) THEN
            ! Resuspend GOCART mass
            IF(p_chem > param_first_scalar) THEN
              resusp_aer = resusp_frac_1d(k) * deposition_flux(iaer) / conv
              chem(i,k,j,p_chem) = MAX(0.0, chem(i,k,j,p_chem) + resusp_aer)
              deposition_flux(iaer) = deposition_flux(iaer) * (1.0 - resusp_frac_1d(k))
            END IF
          ELSE
            CALL wrf_error_fatal('aer_wetscav_simple: unknown chem_opt')
          ENDIF ! aer_is_gocart
        END IF ! DO_RESUSPENSION .AND. ...

      END DO ! k

      !---- Output deposition fluxes
      IF(aer_is_gocart) THEN
        p_chem = aer_pointers(iaer)
        IF(p_chem > param_first_scalar) THEN
          ! WRF-Chem convention is negative from atmosphere to the ground
          qsrflx(i,j,p_chem) = -deposition_flux(iaer)
        ENDIF
      ENDIF

    END DO ! iaer

  END DO ! i
  END DO ! j
  ENDIF ! DO_RAINOUT .OR. DO_WASHOUT

  END SUBROUTINE aer_wetscav_simple

!----------------------------------------------------------------------------

  SUBROUTINE calc_rainout_parameters(rainfrac, frac_rainout, &
     cloud_uptk_effic_liq, cloud_uptk_effic_ice, &
     tk, cloud_fraction, cloud_water, &
     liquid_fraction, rainrate, dtstep)

    ! Purpose:
    ! Calculates parameters for aerosol rainout wet scavenging
    !
    ! Arguments:
    ! Input:
    !   tk = temperature (K)
    !   cloud_fraction = cloud fraction (0-1)
    !   liquid_fraction = fraction of cloud water in liquid form (0-1)
    !   rainrate = precipitation formation rate (kg/kg/s)
    !   dtstep = time step (s)
    ! Input/Output
    !   cloud_water = cloud water mixing ratio (kg/kg)
    ! Output:
    !   rainfrac = fraction of aerosols lost to rainout (unitless)
    !   frac_rainout = rainout area fraction (unitless)
    !   cloud_uptk_effic_liq = efficiency of aerosol uptake in cloud water
    !   cloud_uptk_effic_ice = efficiency of aerosol uptake in cloud ice

    IMPLICIT NONE

    !-- Arguments
    REAL, INTENT(IN) :: tk, cloud_fraction, liquid_fraction, rainrate, dtstep
    REAL, INTENT(INOUT) :: cloud_water
    REAL, INTENT(OUT) :: rainfrac, frac_rainout, cloud_uptk_effic_liq, &
                         cloud_uptk_effic_ice

    !-- Calculate frac_rainout, the precipitating fraction of the grid cell as defined
    ! in Luo et al. (2020)
    ! Also calculate rainfrac, the fraction of rained out aerosols
    rainfrac = 0.0
    frac_rainout = 0.0
    IF(cloud_fraction > AERA_FRAC_MINVAL &
       .AND. cloud_water > WATER_MINVAL &
       .AND. rainrate > 0.0) THEN
      CALL rainfrac_aerosol(rainfrac, frac_rainout, rainrate, dtstep, &
                            cloud_water)
    END IF

    !-- Scaling factors for cloudborne aerosol uptake efficiency, for rainout
    cloud_uptk_effic_liq = 1.0
    cloud_uptk_effic_ice = 0.0
    IF(DO_RAINOUT .AND. cloud_fraction > AERA_FRAC_MINVAL) THEN
      !- Liquid-phase cloudborne aerosol uptake efficiency for rainout
      IF(tk > 237.0) THEN
        cloud_uptk_effic_liq = cloud_fraction
        ! In mixed-phase clouds, correct by the liquid water fraction
        IF(tk < 273. .AND. liquid_fraction > 0.0) THEN
          cloud_uptk_effic_liq = MIN(cloud_uptk_effic_liq,  &
                                    cloud_uptk_effic_liq*liquid_fraction)
        ENDIF
      ENDIF ! tk > 237.0
      !- Ice-phase cloudborne aerosol uptake efficiency for rainout
      IF(tk <= 237.0) THEN
        ! Cold clouds - Uptake by homogeneous nucleation = 100%
        cloud_uptk_effic_ice = cloud_fraction
      ELSEIF(tk < 258.0) THEN
        ! Mixed-phase clouds - Uptake by heterogeneous nucleation
        cloud_uptk_effic_ice = cloud_fraction
        ! Reduce ice-phase scavenging efficiency at warmer temperatures
        ! To take into account the temperature dependence of
        ! heterogeneous nucleation
        cloud_uptk_effic_ice = cloud_uptk_effic_ice &
                             * EXP(0.46*(273.16-tk)-11.6)/153.5
        ! Correct by the ice fraction (1.0-liquid_fraction)
        cloud_uptk_effic_ice = MIN(cloud_uptk_effic_ice, &
                                cloud_uptk_effic_ice &
                                * (1.0 - liquid_fraction))
      ELSE
        ! At high temperatures ice nucleation is improbable so
        ! cloudborne aerosol is 0
        cloud_uptk_effic_ice = 0.0
      ENDIF ! tk <= 237.0
    ENDIF ! DO_RAINOUT .AND. ...

  END SUBROUTINE calc_rainout_parameters

!----------------------------------------------------------------------------

  SUBROUTINE rainfrac_aerosol(rainfrac, frac_rainout, rainrate, dtstep, &
                              cloud_water)

    ! Purpose:
    ! Calculates rainfrac, the fraction of cloudborne aerosols lost to rainout.
    !
    ! Arguments:
    ! Input:
    !   rainrate - precipitation formation rate (kg/kg/s)
    !   dtstep - time step (s)
    !   cloud_water - cloud water mixing ratio (kg/kg)
    ! Output:
    !   rainfrac - fraction of aerosols lost to rainout (unitless)
    !   frac_rainout - rainout area fraction (unitless)
    !
    ! Remarks:
    ! From GEOS-Chem, Luo et al (2020)

    IMPLICIT NONE

    ! Arguments
    REAL, INTENT(IN) :: rainrate
    REAL, INTENT(IN) :: dtstep
    REAL, INTENT(INOUT) :: cloud_water
    REAL, INTENT(OUT) :: rainfrac
    REAL, INTENT(OUT) :: frac_rainout
    ! Local variables and parameters
    REAL :: k_rain ! first-order rainout rate (s-1)
    REAL,  PARAMETER :: k_0 = 1.0E-4 ! s-1

    ! Do not allow more than 100% of the cloud water to precipitate
    cloud_water = MAX(rainrate*dtstep, cloud_water)
    k_rain = k_0 + rainrate / cloud_water
    ! Area of grid box undergoing precipitation
    frac_rainout = rainrate / (k_rain * cloud_water)
    frac_rainout = MIN(frac_rainout, 1.0)
    ! Fraction of aerosols lost to rainout
    rainfrac = frac_rainout * (1.0 - EXP(-k_rain*dtstep))

  END SUBROUTINE rainfrac_aerosol

!----------------------------------------------------------------------------

  SUBROUTINE washfrac_fine_aerosol(washfrac, ppr, pps, ppi, tk, &
                                   frac_washout, dtstep)

    ! Purpose:
    ! Calculates washfrac, the fraction of fine (r_eff < 1 µm) aerosols
    ! lost to washout.
    !
    ! Arguments:
    ! Input:
    !   ppr = rain precipitation through cell (kg/m2/hr = mm/hr)
    !   pps = snow precipitation through cell (kg/m2/hr = mm/hr)
    !   ppi = ice precipitation through cell (kg/m2/hr = mm/hr)
    !   tk = temperature (K)
    !   frac_washout = fraction of cell experiencing washout
    !   dtstep - time step (s)
    ! Output:
    !   washfrac - fraction of aerosols lost to washout (unitless)
    !
    ! Comments:
    ! From GEOS-Chem, Luo et al (2020), using GC12 parameters

    IMPLICIT NONE

    ! Arguments
    REAL, INTENT(IN) :: ppr
    REAL, INTENT(IN) :: pps
    REAL, INTENT(IN) :: ppi
    REAL, INTENT(IN) :: tk
    REAL, INTENT(IN) :: frac_washout
    REAL, INTENT(IN) :: dtstep
    REAL, INTENT(OUT) :: washfrac

    ! Local variables
    REAL :: lossrate_rain,lossrate_snow,lossrate_ice,lossrate

    ! Calculate wet removal rate from each precip type
    lossrate_rain = 4.3E-6 * (ppr/frac_washout)**0.61
    lossrate_snow = 8.8E-6 * (pps/frac_washout)**0.96
    lossrate_ice = 8.8E-6/5.0 * (ppi/frac_washout)**0.96

    ! Calculate total wet removal rate and washout fraction
    lossrate = lossrate_rain+lossrate_snow+lossrate_ice
    washfrac = frac_washout * (1.0 - EXP(-lossrate*dtstep))

  END SUBROUTINE washfrac_fine_aerosol

!----------------------------------------------------------------------------

  SUBROUTINE washfrac_coarse_aerosol(washfrac, ppr, pps, ppi, tk, &
                                     frac_washout, dtstep)

    ! Purpose:
    ! Calculates washfrac, the fraction of coarse (r_eff > 1 µm) aerosols lost
    ! to washout.
    !
    ! Arguments:
    ! Input:
    !   ppr = rain precipitation through cell (kg/m2/hr = mm/hr)
    !   pps = snow precipitation through cell (kg/m2/hr = mm/hr)
    !   ppi = ice precipitation through cell (kg/m2/hr = mm/hr)
    !   tk = temperature (K)
    !   frac_washout = fraction of cell experiencing washout
    !   dtstep - time step (s)
    ! Output:
    !   washfrac - fraction of aerosols lost to washout (unitless)
    !
    ! Comments:
    ! From GEOS-Chem, Luo et al (2020), using GC12 parameters

    IMPLICIT NONE

    ! Arguments
    REAL, INTENT(IN) :: ppr
    REAL, INTENT(IN) :: pps
    REAL, INTENT(IN) :: ppi
    REAL, INTENT(IN) :: tk
    REAL, INTENT(IN) :: frac_washout
    REAL, INTENT(IN) :: dtstep
    REAL, INTENT(OUT) :: washfrac

    ! Local variables
    REAL :: lossrate_rain,lossrate_snow,lossrate_ice,lossrate

    ! Calculate wet removal rate from each precip type
    lossrate_rain = 2.6E-4 * (ppr/frac_washout)**0.79
    lossrate_snow = 4.2E-4 * (pps/frac_washout)**0.96
    lossrate_ice = 4.2E-4/5.0 * (ppi/frac_washout)**0.96

    ! Calculate total wet removal rate and washout fraction
    lossrate = lossrate_rain+lossrate_snow+lossrate_ice
    washfrac = frac_washout * (1.0 - EXP(-lossrate*dtstep))

  END SUBROUTINE washfrac_coarse_aerosol
 
!----------------------------------------------------------------------------

  SUBROUTINE remove_aer(aerconc, dep_flux, lossfrac, conv)

    ! Purpose
    !  Remove fraction lossfrac of aerosol concentration aerconc, and
    !  accumulate the resulting deposition flux to dep_flux
    !
    ! Arguments
    !  aerconc = aerosol concentration (mass/kg or number/kg)
    !  dep_flux = accumulated deposition flux (mass/m2 or number/m2)
    !  lossfrac = fraction of aerosol concentration to remove (0-1)
    !  conv = conversion factor from kg-1 to m-2

    IMPLICIT NONE

    ! Arguments
    REAL, INTENT(INOUT) :: aerconc, dep_flux
    REAL, INTENT(IN) :: lossfrac, conv

    ! Local variables
    REAL :: aerconc_before,aerconc_loss

    ! Deposit aerosols and calculate deposition flux
    aerconc_before = aerconc
    aerconc_loss = lossfrac * aerconc
    aerconc = MAX(0.0, aerconc - aerconc_loss)
    dep_flux = dep_flux + (aerconc_before - aerconc) * conv

  END SUBROUTINE remove_aer

!----------------------------------------------------------------------------

  SUBROUTINE aer_wetscav_simple_init(chem_opt)

    ! Purpose
    !  Initialize aer_wetscav_simple parameters, which depend on chem_opt
    !
    ! Arguments
    !  chem_opt = chemistry option from config_flags, from namelist

    USE module_state_description, ONLY: MOZCART_KPP, &
        T1_MOZCART_KPP, GOCART_SIMPLE, GOCARTRADM2,GOCARTRACM_KPP, &
        p_seas_1, p_seas_2, p_seas_3, p_seas_4, p_bc2, p_oc2, p_dust_1, &
        p_dust_2, p_dust_3, p_dust_4, p_dust_5, p_sulf

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

    !-------- Initialize number of aerosols for chem_opt --------
    IF(aer_is_gocart) THEN
      ! GOCART
      n_aer = 12
    ELSE
      CALL wrf_error_fatal ('wetscav_simple_init: chem_opt not supported')
    ENDIF

    !---- Allocate aerosol mappings and parameters for wet scavenging
    IF (.NOT. ALLOCATED(aer_dry_diam)) THEN
      ALLOCATE(aer_dry_diam(n_aer))
    ENDIF
    IF(aer_is_gocart) THEN
      IF (.NOT. ALLOCATED(deposition_flux)) THEN
        ALLOCATE(deposition_flux(n_aer))
      ENDIF
    ENDIF
    IF (.NOT. ALLOCATED(aer_pointers)) THEN
      ALLOCATE(aer_pointers(n_aer))
    ENDIF
    IF (.NOT. ALLOCATED(aq_aer_ratio)) THEN
      ALLOCATE(aq_aer_ratio(n_aer))
    ENDIF
    IF (.NOT. ALLOCATED(ice_aer_ratio)) THEN
      ALLOCATE(ice_aer_ratio(n_aer))
    ENDIF

    !---- Initialize aerosol mappings and parameters for wet scavenging
    IF(aer_is_gocart) THEN
      !-- GOCART
      ! Species that experience wet scavenging (sulfate deposition needs to
      ! be disabled in mozcart_wetscav)
      aer_pointers(:) = (/ p_sulf, p_seas_1, p_seas_2, p_seas_3, p_seas_4, &
                           p_bc2, p_oc2, p_dust_1, p_dust_2, p_dust_3, &
                           p_dust_4, p_dust_5 /)
      ! Aerosol effective diameters (µm) for aer_pointers species
      aer_dry_diam(:) = (/ 0.312E-6, 0.6E-6, 2.0E-6, 6.5E-6, 15.0E-6, &
                           0.078E-6, 0.174E-6, 1.46E-6, 2.8E-6, 4.8E-6, &
                           9.0E-6, 16.0E-6 /)
      ! Fraction of activated aerosols for aer_pointers species
      ! We do not use the values from Luo et al. (2020) because they are
      ! mostly sensitive to hygroscopicity, even though size is more important
      ! for CCN activation
      ! Instead we assume that only small hygrophobic aerosols (e.g. BC1, OC1)
      ! are inefficient CCN
      aq_aer_ratio(:) = (/ 0.5, 0.5, 0.5, 0.5, 0.5, &
                           0.5, 0.5, 0.5, 0.5, 0.5, &
                           0.5, 0.5 /)
      ! Activated fraction in ice clouds is usually very low except at very low
      ! temperatures, even for efficient INP - pick 0.5 for now but the true
      ! value is probably lower
      ice_aer_ratio(:) = (/ 0., 0., 0., 0., 0., &
                            0., 0., 0.5, 0.5, 0.5, &
                            0.5, 0.5 /)
    ELSE
      CALL wrf_error_fatal ('wetscav_simple_init: chem_opt not supported')
    ENDIF

  END SUBROUTINE aer_wetscav_simple_init

!----------------------------------------------------------------------------

END MODULE module_aer_wetscav_simple

