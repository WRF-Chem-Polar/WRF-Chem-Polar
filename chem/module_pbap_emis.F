! Louis Marelle, Anderson Da Silva, 2025
!
! Purpose:
!  Compute primary biological aerosol particle (PBAP) emissions from
!  continental sources
!
! Remarks:
!  - Pollen not yet included
!  - Only enabled for MODIS landuse. To enable for USGS too, need to define
!    bactflux_from_usgs_land
!  - Could add separate tracers bact1 and spore1 for bacteria and spores
!  - Could emit to MOSAIC tracers oc_a
!
! References
! ==============================================================================
!  (1) Da Silva, A., Marelle, L., Lapere, R., and Raut, J.-C. Modelling Ice
!      Nucleating Particles from High-Latitude Sources to Reproduce Arctic In
!      Situ Concentrations. In preparation, 2026
!  (2) Chatziparaschos et al. (2024) - https://doi.org/10.5194/egusphere-2024-952
!  (3) Hummel et al. (2015) - https://doi.org/10.5194/acp-15-6127-2015
!  (4) Burrows et al. (2009) - https://doi.org/10.5194/acp-9-9281-2009

MODULE module_pbap_emis

  IMPLICIT NONE

  !---- Global module-wide variables
  ! chem_opt aerosol model family
  LOGICAL :: aer_is_gocart
  LOGICAL :: aer_is_mosaic4bin
  LOGICAL :: aer_is_mosaic8bin
  ! Land cover type name
  CHARACTER*256 :: mminlu_loc

  ! Maps MODIS surface categories to Burrows2009 categories for bacterial
  ! number emissions (m-2 s-1)
  ! https://acp.copernicus.org/articles/9/9281/2009/
  REAL, SAVE :: bactflux_from_modis_land(21)
  DATA bactflux_from_modis_land / &
    ! Burrows2009 bact flux      <---- MODIS category
    0., &! forest                   1 'Evergreen Needleleaf Forest'
    0., &! forest                   2 'Evergreen Broadleaf Forest'
    0., &! forest                   3 'Deciduous Needleleaf Forest'
    0., &! forest                   4 'Deciduous Broadleaf Forest'
    0., &! forest                   5 'Mixed Forests'
  502., &! shrubs                   6 'Closed Shrublands'
  502., &! shrubs                   7 'Open Shrublands'
  648., &! grasslands               8 'Woody Savannas'
  648., &! grasslands               9 'Savannas'
  648., &! grasslands              10 'Grasslands'
  196., &! wetlands                11 'Permanent wetlands'
  704., &! crops                   12 'Croplands'
    0., &! unassigned              13 'Urban and Built-Up'
  704., &! crops                   14 'cropland/natural vegetation mosaic'
   7.7, &! landice                 15 'Snow and Ice'
    0., &! deserts                 16 'Barren or Sparsely Vegetated'
    0., &! seas                    17 'Water'
    0., &! tundra                  18 'Wooded Tundra'
    0., &! tundra                  19 'Mixed Tundra'
    0., &! deserts                 20 'Barren Tundra'
    0./  ! unassigned              21 'Lakes'


  CONTAINS

!------------------------------------------------------------------------------

  SUBROUTINE pbap_emis( id, dtstep, dz8w, &
         chem_opt, pbap_opt, &
         rho_phy, chem, e_bact, e_spore,                            &
         xland, ivgtyp, xice, snowh, lai, t2, q2, &
         ids,ide, jds,jde, kds,kde,                                 &
         ims,ime, jms,jme, kms,kme,                                 &
         its,ite, jts,jte, kts,kte                                  )

    ! Purpose:
    !  Calculate Primary Biological Aerosol Particles (PBAP) emissions from
    !  bacteria and fungal spores

    ! Arguments:
    ! Input:
    !  id = domain ID number
    !  dtstep = model time step (s)
    !  dz8w = 3D vertical level depth (m)
    !  chem_opt = chemical mechanism + aerosol mechanism option from config_flags
    !  pbap_opt = pbap emission option from config_flags
    !  rho_phy = 3D air density (kg m-3)
    !  xland = landmask flag (1 land, 2 water)
    !  ivgtyp = land cover category (int category)
    !  xice = sea ice or lake ice fraction (0-1)
    !  snowh = snow height (m)
    !  lai = leaf aera index (m2 m-2)
    !  t2 = 2-m temperature (K)
    !  q2 = 2-m specific humidity (kg kg-1)
    !  ids,ims,its etc. = domain, memory, and tile start and end indices
    ! Input/output
    !  chem = advected chemical species (gases in ppm, aerosols in ug kg-1)
    !  e_bact = accumulated bacteria emission flux (mg m-2)
    !  e_spore = accumulated spore emission flux (mg m-2)

    !TODO Should include all explicit ONLY: associations for module_state_description
    USE module_state_description, ONLY: NUM_CHEM, P_OC1, P_BACT1, P_SPORE1, &
                                        PARAM_FIRST_SCALAR
    USE module_model_constants, ONLY: PICONST
    USE module_data_gocartchem, only: OC_MFAC

    IMPLICIT NONE

    !---- Arguments
    ! Input
    INTEGER, INTENT(IN) :: chem_opt, pbap_opt
    INTEGER, INTENT(IN) :: id,                                      &
                           ids,ide, jds,jde, kds,kde,               &
                           ims,ime, jms,jme, kms,kme,               &
                           its,ite, jts,jte, kts,kte
    INTEGER, DIMENSION(ims:ime, jms:jme), INTENT(IN) :: ivgtyp
    REAL, INTENT(IN) :: dtstep
    REAL, DIMENSION(ims:ime, kms:kme, jms:jme), INTENT(IN) :: dz8w, rho_phy
    REAL, DIMENSION(ims:ime, jms:jme), INTENT(IN) :: xland
    REAL, DIMENSION(ims:ime, jms:jme), INTENT(IN) :: xice, snowh, lai, t2, q2
    ! Input/output
    REAL, DIMENSION(ims:ime, kms:kme, jms:jme, NUM_CHEM), INTENT(INOUT) :: chem
    REAL, DIMENSION(ims:ime, jms:jme), INTENT(INOUT) :: e_bact, e_spore

    !---- Local variables
    ! PBAP emission fluxes (ug m-2 s-1, um3 m-2 s-1, m-2 s-1)
    REAL :: bact_mass_emiss, bact_vol_emiss, bact_num_emiss
    REAL :: spore_mass_emiss, spore_vol_emiss, spore_num_emiss
    ! PBAP active sites density (#/m²)
    REAL :: inas_dens_bact, inas_dens_spore
    ! aerosols area (m²)
    REAL :: aersurf_spore, aersurf_bact
    ! Conversion factor from 1/m2/s to 1/kg_dry_air
    REAL :: conv
    ! Air density (kg/m3)
    REAL :: air_dens
    ! Indices
    INTEGER :: i, j

    !---- Parameters
    ! Aerosol diameter (um)
    REAL, PARAMETER :: AER_DIAM_BACT = 1.0
    REAL, PARAMETER :: AER_DIAM_SPORE = 3.0
    ! Aerosol density (kg m-3)
    REAL, PARAMETER :: BACT_AER_DENS = 1.0E3
    REAL, PARAMETER :: SPORE_AER_DENS = 1.0E3
    ! Parameters for the Hummel et al. (2015) fungal spore emission
    ! parameterization
    REAL, PARAMETER :: B1_HUMMEL15 = 20.426
    REAL, PARAMETER :: B2_HUMMEL15 = 3.93E4

    !-------- Compute pbap emissions --------
    DO i=its,ite
      DO j=jts,jte

        ! Only perform pbap emissions over land. Sea ice is considered land,
        ! so exclude sea ice points explicitly.
        IF (xland(i,j)<1.5 .AND. xice(i,j)<=0.01) THEN
          !-------- Initialize PBAP emission parameters
          bact_num_emiss = 0.0
          bact_vol_emiss = 0.0
          bact_mass_emiss = 0.0
          spore_num_emiss = 0.0
          spore_vol_emiss = 0.0
          spore_mass_emiss = 0.0
          air_dens = rho_phy(i,kts,j)
          ! Conversion factor from m-2 s-1 to kg_dry_air-1
          conv = 1.0/(rho_phy(i,kts,j) * dz8w(i,kts,j)) * dtstep

          !-------- Compute the bacterial aerosol number emission flux
          IF(pbap_opt > 0 &
             .AND. (TRIM(mminlu_loc) == 'MODIS' &
             .OR. TRIM(mminlu_loc) == 'MODIFIED_IGBP_MODIS_NOAH') ) THEN
            ! Bacterial flux from Burrows2009
            bact_num_emiss = bactflux_from_modis_land(ivgtyp(i,j))
            IF(snowh(i,j) > 0.1) THEN
              bact_num_emiss = 7.7 ! m-2 s-1, value for land ice in Burrows2009
            ENDIF
          ELSE
            bact_num_emiss = 0.
          ENDIF

          !-------- Compute the spore aerosol number emission flux
          IF(pbap_opt > 0) THEN
            ! Spore flux from Hummel2015
            spore_num_emiss = B1_HUMMEL15*(t2(i,j)-275.82) + B2_HUMMEL15*q2(i,j)*lai(i,j)
            IF(snowh(i,j) > 0.1) THEN
              spore_num_emiss = 0.
            ENDIF
          ELSE
            spore_num_emiss = 0.
          ENDIF

          !---- Convert PBAP aerosol number emission flux (m-2 s-1) to mass
          ! emission flux (ug m-2 s-1)
          ! 1E-9 for kg m-3 to ug um-3 conversion
          ! Bacteria
          bact_vol_emiss = bact_num_emiss * (PICONST/6.0)*AER_DIAM_BACT**3.0
          bact_mass_emiss = bact_vol_emiss * BACT_AER_DENS * 1.0E-9
          bact_mass_emiss = MAX(bact_mass_emiss, 0.0)
          e_bact(i,j) = e_bact(i,j) + bact_mass_emiss*1.E-3*dtstep
          ! Fungal spores
          spore_vol_emiss = spore_num_emiss * (PICONST/6.0)*AER_DIAM_SPORE**3.0
          spore_mass_emiss = spore_vol_emiss * SPORE_AER_DENS * 1.0E-9
          spore_mass_emiss = MAX(spore_mass_emiss, 0.0)
          e_spore(i,j) = e_spore(i,j) + spore_mass_emiss*1.E-3*dtstep

          !-------- Emit PBAP emissions to chem(i,k,j,l)
          !---- Emit GOCART aerosols to chem
          IF(aer_is_gocart) THEN
            ! GOCART bacteria aerosols
            IF (P_OC1 >= PARAM_FIRST_SCALAR) THEN
              chem(i,kts,j,P_OC1) = chem(i,kts,j,P_OC1) &
                                     + bact_mass_emiss * conv / OC_MFAC
            ENDIF
            ! GOCART fungal spore aerosols
            IF (P_OC1 >= PARAM_FIRST_SCALAR) THEN
              chem(i,kts,j,P_OC1) = chem(i,kts,j,P_OC1) &
                                     + spore_mass_emiss * conv / OC_MFAC
            ENDIF
            ! PBAP number tracers for the INP diagnostic (module_inp_diag).
            ! Least-invasive: the OC1 mass emission above is unchanged; here we
            ! additionally carry the bacteria/spore NUMBER (m-2 s-1 -> #/kg via
            ! conv) so species-specific bio INP can be diagnosed. These advect
            ! and deposit as standard chem tracers.
            IF (P_BACT1 >= PARAM_FIRST_SCALAR) THEN
              chem(i,kts,j,P_BACT1) = chem(i,kts,j,P_BACT1) &
                                     + bact_num_emiss * conv
            ENDIF
            IF (P_SPORE1 >= PARAM_FIRST_SCALAR) THEN
              chem(i,kts,j,P_SPORE1) = chem(i,kts,j,P_SPORE1) &
                                     + spore_num_emiss * conv
            ENDIF
          !---- Other chem_opt cases are not supported
          ELSE
            CALL wrf_error_fatal ('pbap_emis: chem_opt not supported')
          ENDIF ! aer_is_

        END IF ! xland
      END DO ! j
    END DO ! i

  END SUBROUTINE pbap_emis

!------------------------------------------------------------------------------

  SUBROUTINE pbap_emis_init(chem_opt)

    ! Purpose:
    !  Initialize the parameters for the PBAP aerosol emissions
    !  in pbap_emis, which depend on chem_opt

    USE module_state_description ! For chem_opt option values

    IMPLICIT NONE

    ! Arguments:
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

    !-------- Initialize land cover type --------
    CALL nl_get_mminlu( 1, mminlu_loc )
    IF(TRIM(mminlu_loc) == 'MODIS' &
       .OR. TRIM(mminlu_loc) == 'MODIFIED_IGBP_MODIS_NOAH') THEN
      CALL wrf_debug(15, 'pbap_emis_init: Using the following land cover type:')
      CALL wrf_debug(15, TRIM(mminlu_loc))
    ELSE
      CALL wrf_error_fatal('pabap_emis_init: land cover type not supported')
    END IF

  END SUBROUTINE pbap_emis_init

END MODULE module_pbap_emis

