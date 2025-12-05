! 2025/05 Louis Marelle
!
! Purpose:
! Perform wet scavenging of aerosols by grid scale rainout/washout
!
! Remarks:
! Washout and rainout from Luo et al. (2019, 2020)
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
!
!TODO get wet removal fluxes out for validation
!TODO use reasonable thresholds instead of 0.0 in all the module
!TODO remove unused arguments, e.g. dx, dy, config_flags, rho_phy, p_phy
!TODO maybe use double precision
!TODO qc,qi,qr,qs,qg should use b4mp arrays
!

MODULE module_aer_wetscav_simple

   ! USE module_wrf_error, ONLY : wrf_debug

   IMPLICIT NONE
   PRIVATE

   PUBLIC :: aer_wetscav_simple, aer_wetscav_simple_init

   !---- Global module-wide variables
   !-- Allocatable aerosol arrays
   ! Dim is number of aerosol species and depends on aerosol mechanism.
   ! These are initialized in aer_wetscav_simple_init
   REAL, DIMENSION(:), ALLOCATABLE :: mech_aerosol_pointers
   ! Number of wet-scavenged aerosol species for chem_opt (dim of
   ! aer_pointers)
   INTEGER :: naerspecies
   ! List of aerosol species experiencing rainout/washout
   INTEGER, DIMENSION(:), ALLOCATABLE :: aer_pointers
   ! Rainout activation efficiency of species (0-1)
   REAL, DIMENSION(:), ALLOCATABLE :: aq_aer_ratio, ice_aer_ratio
   ! Deposition flux of species during timestep (ug/m2)
   REAL, DIMENSION(:), ALLOCATABLE :: deposition_flux
   ! List of fine aerosol species in mechanism for washout (r<=1 micrometer)
   LOGICAL, DIMENSION(:), ALLOCATABLE :: aer_is_fine
   ! List of coarse aerosol species in mechanism for washout (r>1 micrometer)
   LOGICAL, DIMENSION(:), ALLOCATABLE :: aer_is_coarse

CONTAINS

!----------------------------------------------------------------------------

   SUBROUTINE aer_wetscav_simple( id, dtstep, config_flags,                &
                              t_phy, p_phy, dz8w, qsrflx,                  &
                              chem, rho_phy, cldfra, rainprod, evapprod,   &
                              precr, preci, precs, precg, dx, dy,          &
                              qv, qc, qi, qs, qr, qg,                      &
                              ids,ide, jds,jde, kds,kde,                   &
                              ims,ime, jms,jme, kms,kme,                   &
                              its,ite, jts,jte, kts,kte                    )

   ! Purpose:
   ! Calculates and removes wet scavenging tendencies for aerosols
   !
   ! Arguments:
   ! Input:
   !   id = domain ID
   !   dtstep = model time step (s)
   !   config_flags = namelist configuration flags
   !   t_phy = 3D physical temperature (K)
   !   p_phy = 3D pressure (Pa)
   !   dz8w = 3D vertical level depth (m)
   !   rho_phy = 3D air density (kg/m3)
   !   cldfra = 3D cloud fraction (unitless, 0-1)
   !   rainprod = 3D rain production rate (kg/kg dry air/s)
   !   evapprod = 3D rain evaporation rate (kg/kg dry air/s)
   !   precr = 3D rain precipitation rate through grid cell (kg/m2/s)
   !   preci = 3D ice precipitation rate through grid cell (kg/m2/s)
   !   precs = 3D snow precipitation rate through grid cell (kg/m2/s)
   !   precg = 3D graupel precipitation rate through grid cell (kg/m2/s)
   !   dx, dy = model resolution in x and y (km)
   !   qv = specific humidity = water vapor mixing ratio (kg/kg dry air)
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
   !TODO Use ONLY:
   USE module_configure
   USE module_state_description
   ! p_seas_1, p_seas_2, p_seas_3, p_seas_4, p_bc2, p_bc1, p_oc2, p_dust_1,
   ! p_dust_2, p_dust_3, p_dust_4, p_sulf from module_state_description
   ! USE module_data_gocart_seas, only: reff_seas
   ! USE module_data_gocart_dust, only: reff_dust

   IMPLICIT NONE

   !---- Arguments
   TYPE(grid_config_rec_type), INTENT(IN) :: config_flags
   INTEGER, INTENT(IN) :: ids,ide, jds,jde, kds,kde,    &
                          ims,ime, jms,jme, kms,kme,    &
                          its,ite, jts,jte, kts,kte,    &
                          id
   REAL, INTENT(IN) :: dx, dy
   REAL, INTENT(IN) :: dtstep
   REAL, DIMENSION(ims:ime, kms:kme, jms:jme)         ,        &
       INTENT(IN) :: t_phy, p_phy, rho_phy, dz8w
   REAL, DIMENSION(ims:ime, kms:kme, jms:jme), &
       INTENT(IN) :: qv, qi, qc, qs, qr, qg
   REAL, DIMENSION(ims:ime, kms:kme, jms:jme), &
       INTENT(IN) :: rainprod, evapprod
   REAL, DIMENSION(ims:ime, kms:kme, jms:jme), &
       INTENT(IN) :: precr, preci, precs, precg
   REAL, DIMENSION(ims:ime, kms:kme, jms:jme)         ,        &
       INTENT(IN) :: cldfra
   REAL, DIMENSION(ims:ime, kms:kme, jms:jme, num_chem),          &
       INTENT(INOUT) :: chem
   REAL, DIMENSION(ims:ime, jms:jme, num_chem),          &
         INTENT(OUT) :: qsrflx

   !---- Local variables and parameters
   ! x,y,z, and chemical species index, chemical species pointer
   INTEGER :: i, j, k, nv, p_chem
   ! Cell fractions where rainout/washout is applied (See Luo et al., 2020)
   REAL :: f, fr
   ! Vertical dimension of cldfra for cloud fraction calculation
   INTEGER :: dimz
   ! Aerosol conc. in cloud phase: total cloud, ice-phase only, liquid-phase
   ! only (µg/kg)
   REAL :: cloudborne_aer, cloudborne_aer_ice, cloudborne_aer_liq
   ! Aerosol conc. lost to rainout, aerosols conc. resuspended from precip
   ! (µg/kg)
   REAL :: rainout_aer, resusp_aer
   ! Condensed cloud water mixing ratio (kg/kg dry air), precipitation mixing
   ! ratio (kg/kg dry air), liquid fraction of cloud water (0-1)
   REAL :: cloud_water, precip_water, liquid_fraction
   ! Rainout activation efficiency of species (0-1)
   REAL :: aq_aer_ratio_nv, ice_aer_ratio_nv
   ! Total 3D precipitation trough cell for rain/snow/ice (kg/m2/hr), rain
   ! production rate at i,k,j (kg/kg/s)
   REAL :: ppr, pps, ppi, rainrate
   ! Temperature (K)
   REAL :: tk
   ! Fraction of aerosol species washed out/rained out/resuspended (0-1)
   REAL :: washfrac, rainfrac, resusp_frac
   ! First-order rainout rate (s-1)
   REAL :: k_rain
   ! Value of chem(i,k,j,p) before deposition step (ug/kg)
   REAL :: chem_spec_before
   ! Flags for disabling in-cloud or below-cloud scavenging
   LOGICAL :: DO_RAINOUT, DO_WASHOUT

   DO_RAINOUT = .true.
   DO_WASHOUT = .true.

   !-------- Perform wet scavenging for aerosol species -------
   DO j=jts,jte
     DO i=its,ite
       fr = 0.0
       deposition_flux(:) = 0.0

       DO k=kte-1,kts,-1 ! Do not change, needs to loop in descending order
         !---- Calculate parameters for rainout and washout schemes
         tk = t_phy(i, k, j)
         ! Add rainprod back to cloud_water because rainprod was already
         ! removed from cloud water in microphysics
         cloud_water = qc(i,k,j) + qi(i,k,j) + rainprod(i,k,j)*dtstep
         ! Add evapprod back to precip_water because evapprod was already
         ! removed from precip in microphysics
         precip_water = qr(i,k,j) + qs(i,k,j) + qg(i,k,j) &
                        + evapprod(i,k,j)*dtstep
         IF(cloud_water .GT. 1.0E-12) THEN
           liquid_fraction = MIN(1.0, MAX(0.0, qc(i,k,j)/cloud_water))
         ELSE
           liquid_fraction = 1.0
         ENDIF
         rainrate = rainprod(i,k,j)
         ! Total precipitation through grid cell in kg/m2/hr
         ppr = precr(i,k,j) * 3600.
         pps = precs(i,k,j) * 3600.
         ppi = (preci(i,k,j) + precg(i,k,j)) * 3600.
         ppr = MAX(ppr, 0.)
         pps = MAX(pps, 0.)
         ppi = MAX(ppi, 0.)
         ! Calculate f, the precipitating fraction of the grid cell as defined
         ! in Luo et al. (2020)
         ! Also calculate rainfrac, the fraction of rained out aerosols
         IF(cldfra(i, k, j) .GT. 0.001 .AND. cloud_water .GT. 1.0E-12 &
            .AND. rainrate .GT. 0.0) THEN
           CALL rainfrac_aerosol(rainfrac, f, rainrate, dtstep, cloud_water)
         ELSE
           rainfrac = 0.0
           f = 0.0
         END IF
         ! The grid cell fraction experiencing washout, fr, is the max between f
         ! at k and fr at k+1 (this is why we loop over k in descending order)
         fr = MIN(1.0, MAX(fr, f))

         DO nv = 1,naerspecies
           p_chem = aer_pointers(nv)

           !---- Perform rainout ("in-cloud scavenging")
           IF(DO_RAINOUT .AND. f .GT. 0 .AND. cldfra(i,k,j) .GT. 0.001) THEN

             !-- Calculate the mass concentration of cloudborne aerosols
             ! Liquid-phase cloudborne aerosols
             cloudborne_aer_liq = 0.0
             aq_aer_ratio_nv = aq_aer_ratio(nv)
             IF(tk .GT. 237.0 .AND. aq_aer_ratio_nv .GT. 0.0) THEN
               cloudborne_aer_liq = chem(i,k,j,p_chem) * aq_aer_ratio_nv &
                                * cldfra(i,k,j)
               ! In mixed-phase clouds, correct by the liquid water fraction
               IF(tk .LT. 273.0 .AND. liquid_fraction .GT. 0.0) THEN
                 cloudborne_aer_liq = MIN(cloudborne_aer_liq,  &
                                  chem(i,k,j,p_chem) * aq_aer_ratio_nv &
                                  * cldfra(i,k,j) * liquid_fraction)
               ENDIF
             ENDIF ! aq_aer_ratio_nv .GT. 0.0 .AND. tk .GT. 237.0
             ! Ice-phase cloudborne aerosols
             cloudborne_aer_ice = 0.0
             ice_aer_ratio_nv = ice_aer_ratio(nv)
             IF(tk .LE. 237.0) THEN
               ! 50 % activation efficiency in cold clouds (homogeneous
               ! nucleation)
               ice_aer_ratio_nv = 0.5
               cloudborne_aer_ice = chem(i,k,j,p_chem) * cldfra(i,k,j)
             ELSEIF(tk .LT. 258.0 .AND. ice_aer_ratio_nv .GT. 0.0) THEN
               ! Reduce ice-phase scavenging efficiency at warmer temperatures
               ! To take into account the temperature dependence of
               ! heterogeneous nucleation
               ice_aer_ratio_nv = ice_aer_ratio_nv &
                                  * EXP(0.46*(273.16-tk)-11.6)/153.5
               cloudborne_aer_ice = chem(i,k,j,p_chem) * ice_aer_ratio_nv &
                                    * cldfra(i,k,j)
               ! Correct by the ice fraction (1.0-liquid_fraction)
               cloudborne_aer_ice = MIN(cloudborne_aer_ice,  &
                                cloudborne_aer_ice * (1.0-liquid_fraction))
             ELSE
               ! At high temperatures ice nucleation is improbable so
               ! cloudborne aerosol is 0
               cloudborne_aer_ice = 0.0
             ENDIF ! tk .LE. 237.0

             !-- Appply rainout to chem
             cloudborne_aer = cloudborne_aer_liq + cloudborne_aer_ice
             IF(cloudborne_aer .GT. 0.0) THEN
               rainout_aer = rainfrac * cloudborne_aer
               chem_spec_before = chem(i,k,j,p_chem)
               chem(i,k,j,p_chem) = MAX(0.0, chem(i,k,j,p_chem) - rainout_aer)
               deposition_flux(nv) = deposition_flux(nv) + (chem_spec_before - chem(i,k,j,p_chem)) &
                                     * rho_phy(i,k,j) * dz8w(i,k,j)
             ENDIF
           ENDIF ! f .GT. 0 .AND. cldfra(i, k, j) .GT. 0.001


           !---- Perform washout ("below-cloud scavenging")
           IF(DO_WASHOUT .AND. ppr+pps+ppi .GT. 0.0 .AND. fr .GT. 0.001) THEN
             !-- Calculate washfrac, the fraction of washed out aerosols
             ! Fine aerosol washout (r<=1 micrometer)
             IF(aer_is_fine(nv)) THEN
                CALL washfrac_fine_aerosol(washfrac, ppr, pps, ppi, tk, fr, dtstep)
             ! Coarse aerosol washout (r>1 micrometer)
             ELSEIF(aer_is_coarse(nv)) THEN
                CALL washfrac_coarse_aerosol(washfrac, ppr, pps, ppi, tk, fr, dtstep)
             ! Any other aerosol is not washed out
             ELSE
               washfrac = 0.0
             END IF
             !-- Apply washout to chem
             IF(washfrac .GT. 0.0) THEN
               chem_spec_before = chem(i,k,j,p_chem)
               chem(i,k,j,p_chem) = MAX(0.0, (1.0-washfrac) * chem(i,k,j,p_chem))
               deposition_flux(nv) = deposition_flux(nv) + (chem_spec_before - chem(i,k,j,p_chem)) &
                                     * rho_phy(i,k,j) * dz8w(i,k,j)
             ENDIF
           END IF ! ppr+pps+ppi .GT. 0 .AND. fr .GT. 0.01


           !---- Perform resuspension of aerosols from evaporating precip. water
           ! We apply resuspension to the total deposition flux from the column
           ! above, so this assumes that the timescale of precipitation to level
           ! k is smaller than dtstep
           IF(precip_water .GT. 0. .AND. evapprod(i,k,j) .GT. 0.) THEN
             ! Scale by 0.5 as in GEOS-Chem, to take into account that
             ! evaporation of hydrometeors can be incomplete
             resusp_frac = 0.5 * evapprod(i,k,j) * dtstep / precip_water
             resusp_frac = MAX(0.0, MIN(1.0, resusp_frac))
             resusp_aer = (resusp_frac * deposition_flux(nv) / (rho_phy(i,k,j) * dz8w(i,k,j)))
             chem(i,k,j,p_chem) = MAX(0.0, chem(i,k,j,p_chem) + resusp_aer)
             deposition_flux(nv) = deposition_flux(nv) * (1.0 - resusp_frac)
           END IF

         END DO ! nv
       END DO ! k

       !---- Output deposition fluxes
       DO nv = 1,naerspecies
         p_chem = aer_pointers(nv)
         IF(p_chem > param_first_scalar) THEN
           ! WRF-Chem convention is negative from atmosphere to the ground
           qsrflx(i,j,p_chem) = -deposition_flux(nv)
         ENDIF
       END DO ! nv

     END DO ! i
   END DO ! j

   END SUBROUTINE aer_wetscav_simple

!----------------------------------------------------------------------------

   SUBROUTINE rainfrac_aerosol(rainfrac, f, rainrate, dtstep, cloud_water)

     ! Purpose:
     ! Calculates rainfrac, the fraction of cloudborne aerosols lost to rainout.
     !
     ! Arguments:
     ! Input:
     !   rainrate - precipitation formation rate (kg/kg/s)
     !   dtstep - time step (s)
     !   cloud_water - condensed cloud water mixing ratio (kg/kg)
     ! Output:
     !   rainfrac - fraction of aerosols lost to rainout (unitless)
     !   f - rainout area fraction (unitless)
     !
     ! Remarks:
     ! From GEOS-Chem, Luo et al (2020)

     IMPLICIT NONE

     ! Arguments
     REAL, INTENT(IN) :: rainrate
     REAL, INTENT(IN) :: dtstep
     REAL, INTENT(INOUT) :: cloud_water
     REAL, INTENT(OUT) :: rainfrac
     REAL, INTENT(OUT) :: f
     ! Local variables and parameters
     REAL :: k_rain ! first-order rainout rate (s-1)
     REAL,  PARAMETER :: k_0 = 1.0E-4 ! s-1

     ! Do not allow more than 100% of the cloud water to precipitate
     cloud_water = MAX(rainrate*dtstep, cloud_water)
     k_rain = k_0 + rainrate / cloud_water
     ! Area of grid box undergoing precipitation
     f = rainrate / (k_rain * cloud_water)
     f = MIN(f, 1.0)
     ! Fraction of aerosols lost to rainout
     rainfrac = f * (1.0 - EXP(-k_rain*dtstep))

   END SUBROUTINE rainfrac_aerosol

!----------------------------------------------------------------------------

   SUBROUTINE washfrac_fine_aerosol(washfrac, ppr, pps, ppi, tk, fr, dtstep)

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
     !   fr = fraction of cell experiencing washout
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
     REAL, INTENT(IN) :: fr
     REAL, INTENT(IN) :: dtstep
     REAL, INTENT(OUT) :: washfrac

     ! Local variables
     REAL :: lossrate_rain,lossrate_snow,lossrate_ice,lossrate

     ! Calculate wet removal rate from each precip type
     lossrate_rain = 4.3E-6 * (ppr/fr)**0.61
     lossrate_snow = 8.8E-6 * (pps/fr)**0.96
     lossrate_ice = 8.8E-6/5.0 * (ppi/fr)**0.96

     ! Calculate total wet removal rate and washout fraction
     lossrate = lossrate_rain+lossrate_snow+lossrate_ice
     washfrac = fr * (1.0 - EXP(-lossrate*dtstep))

   END SUBROUTINE washfrac_fine_aerosol

!----------------------------------------------------------------------------

   SUBROUTINE washfrac_coarse_aerosol(washfrac, ppr, pps, ppi, tk, fr, dtstep)

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
     !   fr = fraction of cell experiencing washout
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
     REAL, INTENT(IN) :: fr
     REAL, INTENT(IN) :: dtstep
     REAL, INTENT(OUT) :: washfrac

     ! Local variables
     REAL :: lossrate_rain,lossrate_snow,lossrate_ice,lossrate

     ! Calculate wet removal rate from each precip type
     lossrate_rain = 2.6E-4 * (ppr/fr)**0.79
     lossrate_snow = 4.2E-4 * (pps/fr)**0.96
     lossrate_ice = 4.2E-4/5.0 * (ppi/fr)**0.96

     ! Calculate total wet removal rate and washout fraction
     lossrate = lossrate_rain+lossrate_snow+lossrate_ice
     washfrac = fr * (1.0 - EXP(-lossrate*dtstep))

   END SUBROUTINE washfrac_coarse_aerosol

!----------------------------------------------------------------------------

   SUBROUTINE aer_wetscav_simple_init(chem_opt)

     ! Purpose
     !  Initialize aer_wetscav_simple parameters, which depend on chem_opt
     !
     ! Arguments
     !  chem_opt TODO
     !
     ! Remarks
     !  TODO

     USE module_state_description ! For chem_opt option values

     IMPLICIT NONE

     ! Arguments
     !  Input
     INTEGER, INTENT(IN) :: chem_opt

    !-------- Initialize aerosol tracer number for chem_opt --------
    ! This could be recovered from the aerosol scheme data directly
    IF(chem_opt == MOZCART_KPP .OR. chem_opt == GOCART_SIMPLE) THEN
      ! MOZCART
      naerspecies = 13
    ELSE
      CALL wrf_error_fatal ('wetscav_simple_init: chem_opt not supported')
    ENDIF

    !---- Allocate aerosol mappings and parameters for wet scavenging
    IF (.NOT. ALLOCATED(aer_pointers)) THEN
      ALLOCATE(aer_pointers(naerspecies))
    ENDIF
    IF (.NOT. ALLOCATED(aq_aer_ratio)) THEN
      ALLOCATE(aq_aer_ratio(naerspecies))
    ENDIF
    IF (.NOT. ALLOCATED(ice_aer_ratio)) THEN
      ALLOCATE(ice_aer_ratio(naerspecies))
    ENDIF
    IF (.NOT. ALLOCATED(aer_is_fine)) THEN
      ALLOCATE(aer_is_fine(naerspecies))
    ENDIF
    IF (.NOT. ALLOCATED(aer_is_coarse)) THEN
      ALLOCATE(aer_is_coarse(naerspecies))
    ENDIF
    IF (.NOT. ALLOCATED(deposition_flux)) THEN
      ALLOCATE(deposition_flux(naerspecies))
    ENDIF

    !---- Initialize aerosol mappings and parameters for wet scavenging
    deposition_flux(:) = 0.0
    IF(chem_opt == MOZCART_KPP .OR. chem_opt == GOCART_SIMPLE) THEN
      !-- GOCART
      ! Species that experience wet scavenging (sulfate deposition needs to
      ! be disabled in mozcart_wetscav)
      aer_pointers(:) = (/ p_sulf, p_seas_1, p_seas_2, p_seas_3, p_seas_4, &
                           p_bc2, p_oc2, p_dust_1, p_dust_2, p_dust_3, &
                           p_dust_4, p_dust_5, p_bc1 /)
      ! Fraction of activated aerosols for aer_pointers species
      ! We do not use the values from Luo et al. (2020) because they are
      ! mostly sensitive to hygroscopicity, even though size is more important
      ! for CCN activation
      ! Instead we assume that only small hygrophobic aerosols (e.g. BC1, OC1)
      ! are inefficient CCN
      aq_aer_ratio(:) = (/ 0.5, 0.5, 0.5, 0.5, 0.5, &
                           0.5, 0.5, 0.5, 0.5, 0.5, &
                           0.5, 0.5, 0. /)
      ! Activated fraction in ice clouds is usually very low except at very low
      ! temperatures, even for efficient INP - pick 0.5 for now but it could
      ! be lower
      ice_aer_ratio(:) = (/ 0., 0., 0., 0., 0., &
                            0., 0., 0.5, 0.5, 0.5, &
                            0.5, 0.5, 0. /)
      ! List of washed out fine aerosol for aer_pointers species (r<=1
      ! micrometer)
      aer_is_fine(:) = (/ .true., .true., .true., .false., .false., &
                          .true., .true. , .true., .false., .false., &
                          .false., .false., .false. /)
      ! List of washed out coarse aerosol for aer_pointers species (r>1
      ! micrometer)
      aer_is_coarse(:) = (/ .false., .false., .false., .true., .true., &
                            .false., .false., .false., .true., .true., &
                            .true., .true., .false. /)
    ELSE
      CALL wrf_error_fatal ('wetscav_simple_init: chem_opt not supported')
    ENDIF

   END SUBROUTINE aer_wetscav_simple_init

!----------------------------------------------------------------------------

END MODULE module_aer_wetscav_simple

