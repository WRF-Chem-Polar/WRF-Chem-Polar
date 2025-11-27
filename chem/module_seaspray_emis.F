! Louis Marelle, Rémy Lapere, Lefteris Ioannidis, 2025/03
!
! Purpose:
!  Compute size-resolved sea spray aerosol emissions (mass and number) for the
!  sectional 4-bin or 8-bin MOSAIC aerosols, or for the GOCART 4 sea salt bins.
!  Includes marine organic and sea-salt sulfate emission options.
!
! Remarks:
!  -TODO Marine OA fraction should be size-dependent

! References:
! ==============================================================================
! Source functions:
! (1) Grythe et al., 2014: doi:10.5194/acp-14-1277-2014
! (2) Monahan et al., 1986: doi:10.1007/978-94-009-4668-2_16
! (3) Gong et al., 1997: https://doi.org/10.1029/96JD02953
! (4) Salter et al., 2015: doi:10.5194/acp-15-11047-2015
! (5) Ioannidis et al., 2023: https://doi.org/10.5194/acp-23-5641-2023
! (6) O'Dowd et al., 1997: https://doi.org/10.1016/S1352-2310(96)00106-9
! (7) Lapere et al., 2024: https://doi.org/10.5194/acp-24-12107-2024
! Correction factors:
! (1) Jaegle et al., 2011: doi:10.5194/acp-11-3137-2011
! (2) Salisbury et al., 2013: doi:10.1002/2013JC008797
! (3) Nilsson et al., 2001: https://doi.org/10.1029/2000JD900747
! Sea spray SO4 fraction:
! (1) Calhoun et al., 1991, https://doi.org/10.1029/91GL02304
! Sea spray organic fraction:
! (1) Vignati et al., 2010, https://doi.org/10.1016/j.atmosenv.2009.11.013
!

MODULE module_seaspray_emis

  IMPLICIT NONE

  !---- Global module-wide variables
  ! These are initialized in seaspray_emis_init

  ! Allocatable size-resolved aerosol arrays (dim is size bin number and depends
  ! on aerosol mechanism, unit for vol_emiss_dist is um3/m2/s,
  ! vol_emiss_dist_modN is in um3, bin_diameters_cen in um)
  REAL, DIMENSION(:), ALLOCATABLE :: vol_emiss_dist, vol_emiss_dist_mode1, &
    vol_emiss_dist_mode2, vol_emiss_dist_mode3, bin_diameters_cen

  ! Number of sea spray aerosol bins for chem_opt in use
  INTEGER :: naersizes

  ! chem_opt aerosol model family
  LOGICAL :: aer_is_gocart
  LOGICAL :: aer_is_mosaic4bin
  LOGICAL :: aer_is_mosaic8bin

  CONTAINS

!------------------------------------------------------------------------------

  SUBROUTINE seaspray_emis( id, dtstep, dz8w, chem_opt, seas_opt,    &
         seas_leads_opt, &
         seas_so4_opt, seas_oa_opt,                                 &
         rho_phy, chem, e_nacl_ocean, e_nacl_leads,                 &
         u10, v10, xland, tsk, xice, lakemask, leadfrac, chloroa,   &
         ids,ide, jds,jde, kds,kde,                                 &
         ims,ime, jms,jme, kms,kme,                                 &
         its,ite, jts,jte, kts,kte                                  )

    ! Purpose:
    !  Compute size-resolved sea spray aerosol number and mass emissions

    ! Arguments:
    ! Input:
    !  id = domain ID number
    !  dtstep = model time step (s)
    !  dz8w = 3D vertical level depth (m)
    !  chem_opt = chemical mechanism + aerosol mechanism option from config_flags
    !  seas_opt = sea spray emission option from config_flags
    !  seas_leads_opt = sea spray emission from sea ice leads option from config_flags
    !  seas_so4_opt = sea salt-sulfate emission option from config_flags
    !  seas_oa_opt = marine organic emission option from config_flags
    !  rho_phy = 3D air density (kg/m3)
    !  u10 = windpeed at 10m in x (positive west-east) direction (m/s)
    !  v10 = windpeed at 10m in y (positive south-north) direction (m/s)
    !  xland = landmask flag (1 land, 2 water)
    !  tsk = surface skin temperature (K)
    !  xice = sea ice or lake ice fraction (0-1)
    !  lakemask = lake mask (1 for lakes, 0 for rest)
    !  chloroa = cholorphyll-a concentration at sea surface (mg m-3)
    !  leadfrac = fraction of grid-cell covered by sea ice leads (0-1)
    !  ids,ims,its etc. = domain, memory, and tile start and end indices
    ! Input/output
    !  chem = advected chemical species (gases in ppm, aerosols in microg/kg)
    !  e_nacl_ocean = accumulated NaCl emission flux from sea spray (microg/m2) !TODO should be mg/m2
    !  e_nacl_leads = accumulated NaCl emission flux from sea ice leads (microg/m2) !TODO should be mg/m2

    !TODO Should include all explicit ONLY: associations for module_state_description
    ! chem_opt, seas_opt, param_first_scalar, num_chem, p_seas_1, p_seas_2, p_seas_3, p_seas_4, p_oc1, p_sulf
    USE module_state_description
    USE module_model_constants, ONLY: mwdry, piconst
    !TODO Is it possible to only include this for mosaic?
    USE module_data_mosaic_asect, ONLY: numptr_aer, lptr_cl_aer, &
          lptr_na_aer, lptr_so4_aer, lptr_oc_aer, nsize_aer

    IMPLICIT NONE

    !---- Arguments
    ! Input
    INTEGER, INTENT(IN) :: chem_opt, seas_opt
    INTEGER, INTENT(IN) :: seas_leads_opt, seas_so4_opt, seas_oa_opt
    INTEGER, INTENT(IN) :: id,                                      &
                           ids,ide, jds,jde, kds,kde,               &
                           ims,ime, jms,jme, kms,kme,               &
                           its,ite, jts,jte, kts,kte
    REAL, INTENT(IN) :: dtstep
    REAL, DIMENSION(ims:ime, kms:kme, jms:jme), INTENT(IN) :: dz8w, &
                                                              rho_phy
    REAL, DIMENSION(ims:ime, jms:jme), INTENT(IN) :: u10, v10, xland, &
                                          tsk, xice, lakemask, leadfrac, chloroa
    ! Input/output
    REAL, DIMENSION(ims:ime, kms:kme, jms:jme, num_chem), &
        INTENT(INOUT) :: chem
    REAL, DIMENSION(ims:ime, jms:jme), INTENT(INOUT) :: e_nacl_ocean
    REAL, DIMENSION(ims:ime, jms:jme), INTENT(INOUT) :: e_nacl_leads

    !---- Local variables
    ! Sea spray emission fluxes, microm3/m2/s, #/m2/s
    REAL :: seaspray_vol_emiss, seaspray_num_emiss
    ! Sea salt, sulfate, organic emission fluxes, microg/m2/s
    REAL :: nacl_mass_emiss, so4_mass_emiss, org_mass_emiss
    ! 10-m wind (m/s), sea surface temperature (celsius)
    REAL :: w10, sst
    ! Open ocean fraction in grid cell (0-1)
    REAL :: open_ocean_frac
    ! Whitecap fraction in grid cell (0-1), Salter2015 entrainment factor (0-1)
    REAL :: whitecap_frac, entrainment_factor
    ! Lead fraction in grid cell (0-1)
    REAL :: lead_frac
    ! SST scaling factor (unitless), SST constrained to valid range (celsius)
    REAL :: sst_dependence, sst_valid
    ! Salter emissions number production fluxes for emission modes 1,2,3 (1/m2/s)
    REAL :: salter_n1, salter_n2, salter_n3
    ! aerosol bin center diameter (micro m)
    REAL :: dcen
    ! Conversion factor from 1/m2/s to 1/kg_dry_air
    REAL :: conv
    ! Conversion factor from ug/m2/s to ppmv(sulf) for GOCART
    REAL :: conv_sulf_ppmv
    ! Volume fraction of NaCl, SO4, and OA in sea spray (0-1, m3_spec/m3_ssa)
    REAL :: nacl_vol_frac,so4_vol_frac,org_vol_frac
    ! Conversion factor from sea spray volume to nacl mass (ug(nacl)/um3(seaspray))
    REAL :: nacl_mass_to_ss_vol_ratio
    ! Conversion factor from sea spray volume to sulfate mass (ug(so4)/um3(seaspray))
    REAL :: so4_mass_to_ss_vol_ratio
    ! Conversion factor from sea spray volume to organic aerosol mass (ug(oa)/um3(seaspray))
    REAL :: org_mass_to_ss_vol_ratio
    ! Sea-salt-sulfate/NaCl mass ratio (kg/kg) and Marine-organics/NaCl mass ratio (kg/kg)
    REAL :: so4_nacl_mass_ratio, org_nacl_mass_ratio
    ! Vignati org/(org+nacl) mass ratio (kg/kg)
    REAL :: org_mass_ratio
    ! Sea ice lead emission parameterization
    ! Flux ratio between sea-spray from open leads in sea-ice and sea-spray
    ! from the ice-free open ocean, from Nilsson et al. (2001)
    REAL :: leads_to_oo_ratio

    ! Indices
    INTEGER :: i, j, isize, iphase
    INTEGER :: iaerbin
    INTEGER :: l_na, l_cl, l_so4, l_num, l_seasalt, l_oc

    !---- Parameters
    !TODO should be defined in init and/or data module
    ! REAL, PARAMETER :: NACL_AERDEN = 2.17 * 1.0E6 * 1.0E-12 ! g/cm3 -> microg/cm3 -> microg/um3
    ! REAL, PARAMETER :: SO4_AERDEN = 1.8 * 1.0E6 * 1.0E-12 ! g/cm3 -> microg/cm3 -> microg/um3, MOSAIC value
    REAL, PARAMETER :: NACL_AERDEN = 2.2E-6 ! microg/um3
    REAL, PARAMETER :: SO4_AERDEN = 1.8E-6 ! microg/um3, MOSAIC value
    REAL, PARAMETER :: ORG_AERDEN = 1.8E-6 ! microg/um3, MOSAIC value
    REAL, PARAMETER :: MW_NA = 22.99, MW_CL = 35.45, MW_SO4 = 96.06, MW_ORG = 12.011
    REAL, PARAMETER :: GOCART_OC_MFAC = 1.8 ! kg(OA)/kg(OC)
    ! Adjustment factor for Salter emissions mode 3
    ! The base Salter et al.(2015) parameterization overestimates SSA PM10
    ! by a significant amount (x3 to x5).
    ! The difference is not too high for small particles (as shown in S2015 Figure
    ! 5b) so the issue is with coarse mode 3.
    ! The supermicron PM ratio between Salter and Gong goes between ~1.98 at
    ! low temperatures to ~5.28 at high temperatures, so we chose to adjust
    ! by 1/2.5, which is the ratio in the 10C/20C range
    REAL, PARAMETER :: salter_mode3_corr = 2.5
    ! Polynomial parameters for the Salter2015 emissions SST factors
    REAL, PARAMETER :: A1=-5.2168E5, B1=3.31725E7, C1=-6.95275E8, D1=1.0684E10
    REAL, PARAMETER :: A2=0.0, B2=7.374E5, C2=-2.4803E7, D2=7.7373E8
    REAL, PARAMETER :: A3=0.0, B3=1.4210E4, C3=1.4662E7, D3=1.7075E8
    ! Logicals
    LOGICAL :: DO_SEAS_SO4
    LOGICAL :: DO_SEAS_ORG

    !-------- Compute SSA emissions --------
    ! Activate Marine sulfate emissions based on seas_so4_opt
    IF(seas_so4_opt > 0.0) THEN
       DO_SEAS_SO4 = .true.
    ELSE
       DO_SEAS_SO4 = .false.
    ENDIF
    ! Activate Marine organic emissions based on seas_oa_opt
    IF(seas_oa_opt > 0.0) THEN
       DO_SEAS_ORG = .true.
    ELSE
       DO_SEAS_ORG = .false.
    ENDIF

    ! No emissions to cloud-phase
    !TODO Should be = aiphase but might not be initialized
    iphase = 1

    DO i=its,ite
      DO j=jts,jte

        ! Only perform sea spray emissions over open sea water (xland>1.5,
        ! lakemask(i,j)<0.5)
        ! Sea ice is considered land, so include fractional sea ice points
        ! explicitly. We cannot use xice(i,j) .LT. 0.0 on its own because this
        ! is true for all land points.
        IF ( (xland(i,j) .GT. 1.5 &
            .OR. (xice(i,j) .GT. 0.01 .AND. xice(i,j) .LT. 1.0)) &
            .AND. lakemask(i,j) .LT. 0.5) THEN

          !-------- Initialize sea spray emission parameters --------
          nacl_mass_emiss = 0.0
          so4_mass_emiss = 0.0
          org_mass_emiss = 0.0
          whitecap_frac = 0.0
          entrainment_factor = 0.0
          sst_dependence = 1.0
          ! 10-m wind speed (m/s)
          w10 = SQRT(u10(i, j)**2.0 + v10(i, j)**2.0)
          ! Conversion factor from 1/m2/s to 1/kg_dry_air
          conv = 1.0/(rho_phy(i,kts,j) * dz8w(i,kts,j)) * dtstep
          !---- Open ocean fraction
          ! When fractional sea ice is present, scale the emission flux by the
          ! open ocean fraction
          open_ocean_frac = 1.0 - xice(i,j)
          !---- Whitecap fraction, used later in sea spray source functions
          IF(seas_opt == SEASIOANNIDIS) THEN
            ! Whitecap fraction from Salisbury et al. (2013), recommended by
            ! Ioannidis et al. (2023)
            whitecap_frac = 4.60E-5 * w10**2.26
          ELSEIF(seas_opt == SEASSALTER) THEN
            entrainment_factor = 2.0E-8 * w10**3.41
          ELSE
            ! Whitecap fraction from Monahan et al. (1986)
            whitecap_frac = 3.84E-6 * w10**3.41
          ENDIF ! seas_opt
          !---- SST correction factors
          IF(seas_opt == SEASMONAHAN .OR. seas_opt == SEASIOANNIDIS) THEN
            ! SST correction from Jaegle et al. (2011), valid between -2 C and
            ! 30 C, recommended by Ioannidis et al. (2023)
            sst_valid = MIN(MAX(tsk(i,j)-273.15, -2.0), 30.0)
            sst_dependence = 0.3 + 0.1*sst_valid - 0.0076*sst_valid**2.0 + 0.00021*sst_valid**3.0
          ELSEIF(seas_opt == SEASSALTER) THEN
            ! Salter emissions number fluxes for modes 1,2,3 are SST-dependent
            sst_valid = MIN(MAX(tsk(i,j)-273.15, 2.0), 30.0)
            salter_n1 = A1*sst_valid**3.0 + B1*sst_valid**2.0 + C1*sst_valid + D1
            salter_n2 = A2*sst_valid**3.0 + B2*sst_valid**2.0 + C2*sst_valid + D2
            salter_n3 = A3*sst_valid**3.0 + B3*sst_valid**2.0 + C3*sst_valid + D3
          ENDIF ! seas_opt
          !---- Correction factor for the flux from sea ice leads
          ! Following Lapere et al. (2024)
          leads_to_oo_ratio = 1.0
          !- seas_leads_opt=0 sea_spray*(1-sea_ice) everywhere
          !- seas_leads_opt=1 Nilsson ratio average in leads
          lead_frac = leadfrac(i,j)
          IF ((seas_leads_opt == 1 .AND. lead_frac > 0.0)) THEN
             ! oo flux versus lead flux in Nilsson 2001
             leads_to_oo_ratio = EXP(0.11*w10-1.93)/EXP(0.20*w10-1.71)
             ! Recalculate the open ocean fraction to remove the contribution
             ! from leads
             open_ocean_frac = MIN(1., MAX(0., open_ocean_frac - lead_frac))
          ENDIF

          !---- Volume fractions of organics, so4 and nacl in sea spray
          IF(DO_SEAS_SO4) THEN
            ! so4/na mass ratio = 0.252, Ioannidis et al.( 2023)
            so4_nacl_mass_ratio = 0.252 * MW_NA / (MW_NA + MW_CL)
          ELSE
            so4_nacl_mass_ratio = 0.0
          ENDIF
          IF(DO_SEAS_ORG) THEN
            ! Organic mass fraction from Vignati et al. 2010
            org_mass_ratio = MAX(2.0, MIN(43.5 * chloroa(i,j) + 13.805, 76.0)) / 100.0
            ! Convert from organic mass fraction to organic/nacl mass fraction
            org_nacl_mass_ratio = org_mass_ratio / (1. - org_mass_ratio)
          ELSE
            org_nacl_mass_ratio = 0.0
          ENDIF
          ! Calculate volume fractions (vol_spec/vol_ssa)
          nacl_vol_frac = 1./(1 + so4_nacl_mass_ratio*NACL_AERDEN/SO4_AERDEN &
                              + org_nacl_mass_ratio*NACL_AERDEN/ORG_AERDEN)
          org_vol_frac = nacl_vol_frac*org_nacl_mass_ratio*NACL_AERDEN/ORG_AERDEN
          so4_vol_frac = nacl_vol_frac*so4_nacl_mass_ratio*NACL_AERDEN/SO4_AERDEN

          !-------- Compute sea-spray aerosol mass/num emissions for each aer
          !         size bin --------
          DO iaerbin = 1,naersizes
            !---- Compute sea-spray aerosol  mass and number emission fluxes
            ! Calculate the volume emission flux seaspray_vol_emiss(iaerbin)
            ! (um3/m2/s), from size-resolved volume emission scaling factors
            ! vol_emiss_dist(iaerbin) (unit is source function dependent,
            ! um3/m2/s for most source functions, um3 for Salter et al. 2015)
            IF(seas_opt == SEASSALTER) THEN
              ! Salter2015 emissions have 3 lognormal modes scaled independently
              ! as a function of SST
              seaspray_vol_emiss = entrainment_factor &
                          *(salter_n1*vol_emiss_dist_mode1(iaerbin) &
                          + salter_n2*vol_emiss_dist_mode2(iaerbin) &
                          + salter_n3*vol_emiss_dist_mode3(iaerbin))
            ELSE
              ! Other source functions have the same w10 and SST scaling for all
              ! sizes and only one size distribution is needed and precalculated in
              ! vol_emiss_dist
              seaspray_vol_emiss = sst_dependence*whitecap_frac &
                          *vol_emiss_dist(iaerbin)
            ENDIF ! seas_opt
            !---- Correct sea-spray flux by open ocean fraction or open leads
            !  fraction
            IF(seas_leads_opt == 1 .AND. lead_frac > 0.) THEN
              seaspray_vol_emiss = seaspray_vol_emiss*open_ocean_frac &
                              + seaspray_vol_emiss*leads_to_oo_ratio*lead_frac
            ELSE
              seaspray_vol_emiss = seaspray_vol_emiss*open_ocean_frac
            ENDIF
            ! Do not allow negative emissions
            seaspray_vol_emiss = MAX(seaspray_vol_emiss, 0.0)
            !---- Convert sea spray aerosol volume emission flux (microm3/m2/s)
            !  to speciated aerosol mass emission flux (microg/m2/s)
            nacl_mass_emiss = seaspray_vol_emiss * nacl_vol_frac * NACL_AERDEN
            IF(DO_SEAS_SO4) THEN
              so4_mass_emiss = seaspray_vol_emiss * so4_vol_frac * SO4_AERDEN
            ENDIF
            IF(DO_SEAS_ORG) THEN
               org_mass_emiss = seaspray_vol_emiss * org_vol_frac * ORG_AERDEN
            ENDIF
            !---- Convert sea spray aerosol volume emission flux to number
            ! emission flux (#/m2/s), assuming aerosol volume is the volume at
            ! bin log-center diameter dcen
            dcen = bin_diameters_cen(iaerbin)
            seaspray_num_emiss = seaspray_vol_emiss / ((piconst/6.0)*dcen**3.0)
            seaspray_num_emiss = MAX(seaspray_num_emiss, 0.0)
            !---- Accumulate sea-spray nacl emissions from open ocean, output diagnosis
            e_nacl_ocean(i, j) = e_nacl_ocean(i, j) + nacl_mass_emiss * dtstep


            !-------- Emit sea spray emissions (add them to chem(i,k,j,l)) --------
            !---- Emit MOSAIC aerosols to chem
            IF(aer_is_mosaic4bin .OR. aer_is_mosaic8bin) THEN
              ! MOSAIC sodium aerosol
              l_na = lptr_na_aer(iaerbin,1,iphase)
              IF (l_na >= param_first_scalar) THEN
                chem(i,kts,j,l_na) = chem(i,kts,j,l_na) + nacl_mass_emiss * MW_NA/(MW_NA+MW_CL) * conv
              END IF
              ! MOSAIC chlorine aerosol
              l_cl = lptr_cl_aer(iaerbin,1,iphase)
              IF (l_cl >= param_first_scalar) THEN
                chem(i,kts,j,l_cl) = chem(i,kts,j,l_cl) + nacl_mass_emiss * MW_CL/(MW_NA+MW_CL) * conv
              END IF
              ! MOSAIC sulfate aerosol
              IF(DO_SEAS_SO4) THEN
                l_so4 = lptr_so4_aer(iaerbin,1,iphase)
                IF (l_cl >= param_first_scalar) THEN
                  chem(i,kts,j,l_so4) = chem(i,kts,j,l_so4) + so4_mass_emiss * conv
                END IF
              ENDIF
              ! MOSAIC oc aerosol
              IF(DO_SEAS_ORG) THEN
                l_oc = lptr_oc_aer(iaerbin,1,iphase)
                IF (l_oc >= param_first_scalar) THEN
                   chem(i,kts,j,l_oc) = chem(i,kts,j,l_oc) + org_mass_emiss * conv
                END IF
              END IF
              ! MOSAIC aerosol number
              l_num = numptr_aer(iaerbin,1,iphase)
              IF (l_num >= param_first_scalar) THEN
                chem(i,kts,j,l_num) = chem(i,kts,j,l_num) + seaspray_num_emiss * conv
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
              ! GOCART OC aerosol
              IF(DO_SEAS_ORG) THEN
                IF (p_oc1 >= param_first_scalar) THEN
                  chem(i,kts,j,p_oc1) = chem(i,kts,j,p_oc1) + org_mass_emiss / GOCART_OC_MFAC * conv
                ENDIF
              ENDIF
              ! GOCART sulf aerosol
              IF(DO_SEAS_SO4) THEN
                ! Convert ug/m2/s -> ug/kg -> mg/kg (1E-6 g/g) -> ppmv (1E-6
                ! mol/mol)
                conv_sulf_ppmv =  conv * 1E-3 * mwdry/MW_SO4
                IF (p_sulf >= param_first_scalar) THEN
                  chem(i,kts,j,p_sulf) = chem(i,kts,j,p_sulf) + so4_mass_emiss * conv_sulf_ppmv
                ENDIF
              ENDIF
            !---- Other chem_opt cases are not supported
            ELSE
              CALL wrf_error_fatal ('seaspray_emis: chem_opt not supported')
           ENDIF ! aer_is_


          END DO ! iaerbin

        END IF ! xland
      END DO ! j
    END DO ! i

  END SUBROUTINE seaspray_emis

!------------------------------------------------------------------------------

  SUBROUTINE seaspray_emis_init(chem_opt, seas_opt)

    ! Purpose:
    !  Initialize the parameters for the sea spray aerosol source functions used
    !  in seaspray_emis, which depend on chem_opt and seas_opt
    !
    ! Remarks:
    ! - vol_emiss_dist was precalculated by integrating the source functions
    !   size distributions from the original publications over the WRF-Chem
    !   aerosol size ranges. For this reason if the aerosol bin limits are changed
    !   vol_emiss_dist needs to be recalculated.
    ! - TODO Center diameters and naersizes should be compared to the ones used at
    !   runtime (for MOSAIC, dcen_sect(iaerbin,itype) and nsize_aer(itype) from
    !   module_data_mosaic_asect), and model should crash if there is a mismatch
    !   - TODO If comparing to MOSAIC nsizes and dcen_sect, careful the MOSAIC top
    !     level initializations should have been done before.
    ! - TODO Allocatables here are never deallocated
    ! - aer_is and data declaration must be updated if a new mechanism is defined

    !TODO only use with only:, specify all variables here
    USE module_state_description ! For chem_opt and seas_opt option values

    IMPLICIT NONE

    ! Arguments:
    !  Input
    INTEGER, INTENT(IN) :: chem_opt, seas_opt

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

    !-------- Initialize sea spray aerosol bin number for chem_opt --------
    ! This could be recovered from the aerosol scheme data directly, but
    ! vol_emiss_dist is already chem_opt-dependent and depends on naersizes so
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
      CALL wrf_error_fatal ('seaspray_emis_init: chem_opt not supported')
    ENDIF

    !-------- Initialize size-resolved aerosol arrays for SSA emissions --------
    ! Initialize vol_emiss_dist (microm3/m2/s) and bin_diameters_cen
    ! (micrometer) for each family of aerosol options selected by
    ! chem_opt

    ! Allocate the size-resolved aerosol arrays
    IF (.NOT. ALLOCATED(vol_emiss_dist)) THEN
      ALLOCATE(vol_emiss_dist(naersizes))
    ENDIF
    IF (.NOT. ALLOCATED(vol_emiss_dist_mode1)) THEN
      ALLOCATE(vol_emiss_dist_mode1(naersizes))
    ENDIF
    IF (.NOT. ALLOCATED(vol_emiss_dist_mode2)) THEN
      ALLOCATE(vol_emiss_dist_mode2(naersizes))
    ENDIF
    IF (.NOT. ALLOCATED(vol_emiss_dist_mode3)) THEN
      ALLOCATE(vol_emiss_dist_mode3(naersizes))
    ENDIF
    IF (.NOT. ALLOCATED(bin_diameters_cen)) THEN
      ALLOCATE(bin_diameters_cen(naersizes))
    ENDIF

    !-- Initialize SSA source functions for GOCART aerosols (4 sea salt bins)
    IF(aer_is_gocart) THEN
      ! GOCART sea salt bin center diameters (micrometer)
      ! These should be log-averaged but I follow the GOCART definition instead
      bin_diameters_cen(:) = (/ 0.60, 2.00, 6.5, 15.0 /)
      ! Size-resolved volume emission scaling factors for GOCART
      IF(seas_opt == SEASMONAHAN) THEN
        ! Monahan et al. (1986) SSA source function
        vol_emiss_dist(:) = (/ 367221., 3534527., 2057189., 3732253. /)
      ELSEIF(seas_opt == SEASIOANNIDIS .OR. seas_opt == SEASGONG) THEN
        ! Gong et al. (1997) + O'Dowd et al. (1997) SSA source function
        ! Same as Monahan above 100 nm
        vol_emiss_dist(:) = (/ 367221., 3534527., 2057189., 3732253. /)
      ELSEIF(seas_opt == SEASSALTER) THEN
        ! Salter et al. (2015) SSA source function (3 lognormal modes)
        vol_emiss_dist_mode1(:) = (/ 3.83446671e-03, 8.78064403e-04, 4.04923045e-05, 1.35064728e-07 /)
        vol_emiss_dist_mode2(:) = (/ 1.04739718e-01, 2.81827227e-01, 3.81633973e-02, 7.82888780e-05 /)
        vol_emiss_dist_mode3(:) = (/ 0.05502078, 2.45584876, 2.24370337, 0.02050615 /)
      ELSE
        CALL wrf_error_fatal ('seaspray_emis_init: seas_opt not supported')
      ENDIF

    !-- Initialize SSA source functions for MOSAIC-8bin aerosols
    ELSEIF(aer_is_mosaic8bin) THEN
      ! MOSAIC-8bin bin center diameters (micrometer)
      bin_diameters_cen(:) = (/ 0.0552, 0.1105, 0.2210, 0.4419, &
                             0.8839, 1.7678, 3.5355, 7.0711 /)
      ! Size-resolved volume emission scaling factors for MOSAIC-8bin
      IF(seas_opt == SEASMONAHAN) THEN
        ! Monahan et al. (1986) SSA source function
        vol_emiss_dist(:) = (/ 0., 10673., 29885., 72482., &
                            786274., 2695160., 977044., 1406417. /)
      ELSEIF(seas_opt == SEASIOANNIDIS .OR. seas_opt == SEASGONG) THEN
        ! Gong et al. (1997) + O'Dowd et al. (1997) SSA source function
        vol_emiss_dist(:) = (/ 5037., 14721., 29885., 72482., &
                            786274., 2695160., 977044., 1406417. /)
      ELSEIF(seas_opt == SEASSALTER) THEN
        ! Salter et al. (2015) SSA source function (3 lognormal modes)
        vol_emiss_dist_mode1(:) = (/ 3.25430278e-05, 2.86242430e-04, 1.10845502e-03, 1.90106720e-03, &
                                  1.44815452e-03, 4.89275363e-04, 7.29402899e-05, 4.76220834e-06 /)
        vol_emiss_dist_mode2(:) = (/ 1.52923262e-08, 8.46400525e-06, 9.80949397e-04, 2.46583115e-02, &
                                  1.40982898e-01, 1.91288946e-01, 6.21746746e-02, 4.69120329e-03 /)
        vol_emiss_dist_mode3(:) = (/ 3.32393538e-14, 1.16966419e-09, 4.91681378e-06, 2.53850868e-03, &
                                  1.69813539e-01, 1.61016157e+00, 2.39700335e+00, 5.75051048e-01 /)
      ELSE
        CALL wrf_error_fatal ('seaspray_emis_init: seas_opt not supported')
      ENDIF

    !-- Initialize SSA source functions for MOSAIC-4bin aerosols
    ELSEIF(aer_is_mosaic4bin) THEN
      ! MOSAIC-4bin bin center diameters (micrometer)
      bin_diameters_cen(:) = (/ 0.0781, 0.3125, 1.2500, 5.0000 /)
      ! Size-resolved volume emission scaling factors for MOSAIC-4bin
      IF(seas_opt == SEASMONAHAN) THEN
        ! Monahan et al. (1986) SSA source function
        vol_emiss_dist(:) = (/ 10673., 102367., 3481434., 2383461. /)
      ELSEIF(seas_opt == SEASIOANNIDIS .OR. seas_opt == SEASGONG) THEN
        ! Gong et al. (1997) + O'Dowd et al. (1997) SSA source function
        vol_emiss_dist(:) = (/ 19757., 102367., 3481434., 2383461./)
      ELSEIF(seas_opt == SEASSALTER) THEN
        ! Salter et al. (2015) SSA source function (3 lognormal modes)
        vol_emiss_dist_mode1(:) = (/ 3.18830776e-04, 3.00952222e-03, 1.93742988e-03, 7.77024982e-05 /)
        vol_emiss_dist_mode2(:) = (/ 8.48255180e-06, 2.56392609e-02, 3.32271843e-01, 6.68658779e-02 /)
        vol_emiss_dist_mode3(:) = (/ 1.17045932e-09, 2.54342549e-03, 1.77997511e+00, 2.97205440e+00 /)
      ELSE
        CALL wrf_error_fatal ('seaspray_emis_init: seas_opt not supported')
      ENDIF

    ELSE ! aer_is_
      CALL wrf_error_fatal ('seaspray_emis_init: chem_opt not supported')

    ENDIF ! aer_is_

  END SUBROUTINE seaspray_emis_init

END MODULE module_seaspray_emis

