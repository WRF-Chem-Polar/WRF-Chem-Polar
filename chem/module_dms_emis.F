! 2026 Louis Marelle, Rémy lapere
!
! Purpose:
!  Compute DMS emissions
!
! Remarks:
!  - This requires that DMS_OCEAN is present in the wrflowinp file and that
!    the DMS tracer is defined in chem for this mechanism
!
! References:
! ==============================================================================
! Air-sea exchange functions:
!   Nightingale et al., 2000: https://doi.org/10.1029/1999GB900091
!   Liss and Merlivat, 1986: https://doi.org/10.1007/978-94-009-4738-2_5
!   Wanninkhof, 2014: https://doi.org/10.4319/lom.2014.12.351
! Schmidt number for DMS
!   Saltzman et al., 1993: https://doi.org/10.1029/93JC01858
!   Wanninkhof, 2014: https://doi.org/10.4319/lom.2014.12.351
! Henry's Law Coefficient for DMS
!   Sander, 2015: https://doi.org/10.5194/acp-15-4399-2015
! Sea ice correction factor
!   Loose et al., 2009: doi:10.1029/2008GL036318
! Implementation and evaluation:
!   Marelle et al., 2017: doi:10.5194/gmd-10-3661-2017
!   Lapere et al., 2025: https://dx.doi.org/10.22541/au.175466602.22534100/v1

MODULE module_dms_emis

  IMPLICIT NONE

  !---- Parameters
  ! Schmidt number for co2 at 20°C = 600 (reference value)
  REAL, PARAMETER :: SC_CO2 = 600.

  CONTAINS

!------------------------------------------------------------------------------

  SUBROUTINE dms_emis( id, dtstep, dz8w, dms_opt,                   &
         rho_phy, p_phy, chem, e_dms_ocean,                         &
         u10, v10, xland, tsk, xice, lakemask, dms_ocean,           &
         ids,ide, jds,jde, kds,kde,                                 &
         ims,ime, jms,jme, kms,kme,                                 &
         its,ite, jts,jte, kts,kte                                  )

    ! Purpose:
    !  Compute DMS emissions

    ! Arguments:
    ! Input:
    !  id = domain ID number
    !  dtstep = model time step (s)
    !  dz8w = 3D vertical level depth (m)
    !  dms_opt = dms emission option from config_flags
    !  rho_phy = 3D air density (kg/m3)
    !  p_phy = 3D air pressure (Pa)
    !  u10 = windpeed at 10m in x (positive west-east) direction (m/s)
    !  v10 = windpeed at 10m in y (positive south-north) direction (m/s)
    !  xland = landmask flag (1 land, 2 water)
    !  tsk = surface skin temperature (K)
    !  xice = sea ice or lake ice fraction (0-1)
    !  lakemask = lake mask (1 for lakes, 0 for rest)
    !  dms_ocean = oceanic concentration of DMS in micro mol/m3
    !  ids,ims,its etc. = domain, memory, and tile start and end indices
    ! Input/output
    !  chem = advected chemical species (gases in ppm, aerosols in microg/kg)
    !  e_dms_ocean = accumulated DMS emissions (mol/m2)

    USE module_state_description

    IMPLICIT NONE

    !---- Arguments
    ! Input
    INTEGER, INTENT(IN) :: dms_opt
    INTEGER, INTENT(IN) :: id,                                      &
                           ids,ide, jds,jde, kds,kde,               &
                           ims,ime, jms,jme, kms,kme,               &
                           its,ite, jts,jte, kts,kte
    REAL, INTENT(IN) :: dtstep
    REAL, DIMENSION(ims:ime, kms:kme, jms:jme), INTENT(IN) :: dz8w, &
                                                    rho_phy, p_phy
    REAL, DIMENSION(ims:ime, jms:jme), INTENT(IN) :: u10, v10, xland, &
                                          tsk, xice, lakemask, dms_ocean
    ! Input/output
    REAL, DIMENSION(ims:ime, kms:kme, jms:jme, num_chem), &
        INTENT(INOUT) :: chem
    REAL, DIMENSION(ims:ime, jms:jme), INTENT(INOUT) :: e_dms_ocean

    !---- Local variables
    ! DMS emission flux (mol km-2 hr-1)
    REAL :: dms_emiss
    ! Oceanic DMS content in water (mol m-3)
    REAL :: c_dms_water
    ! Oceanic DMS content in water-air interface (mol m-3)
    REAL :: c_dms_interface
    ! Henry's Law coefficient for DMS (mol m-3 Pa-1)
    REAL :: h_dms
    ! DMS partial pressure in air (Pa)
    REAL :: press_dms
    ! Schmidt number for DMS (dimensionless)
    REAL :: sc_dms
    ! Sea-air exchange coefficient for DMS (cm hour-1)
    REAL :: kw_dms
    ! 10-m wind (m/s), sea surface temperature (celsius)
    REAL :: w10, sst
    ! Conversion factor from mol km-2 hr-1 to ppmv
    REAL :: conv
    ! Sea ice fraction
    REAL :: icefrac
    ! Indices
    INTEGER :: i, j

    !-------- Compute DMS emissions --------

    DO i=its,ite
      DO j=jts,jte

        ! Only perform DMS emissions over open water (xland>1.5),
        ! Sea ice is considered land, so include fractional sea ice points
        ! explicitly. We cannot use xice(i,j) .LT. 0.0 on its own because this
        ! is true for all land points.
        icefrac = xice(i,j)
        IF ( (xland(i,j) .GT. 1.5 &
            .OR. (icefrac .GT. 1.E-9 .AND. icefrac .LT. 1.0)) ) THEN

          !-------- Initialize DMS emission parameters --------
          dms_emiss = 0.0
          ! DMS concentration in water - nmol L-1 to mol m-3
          c_dms_water = dms_ocean(i,j)*1.E-6

          IF(c_dms_water > 0.) THEN
            ! 10-m wind speed (m/s)
            w10 = SQRT(u10(i, j)**2.0 + v10(i, j)**2.0)
            ! Conversion factor from mol km-2 hr-1 to ppmv
            conv = 4.828e-4 / rho_phy(i,kts,j) * dtstep / (dz8w(i,kts,j) * 60.)
            ! Sea surface temperature, from skin temperature
            sst = tsk(i,j)


            !---- Calculate the sea-air exchange coefficient
            kw_dms = 0.
            sc_dms = 0.
            ! Sea-air exchange coefficient for DMS from Nightingale et al. (2000)
            IF(dms_opt == DMS_NIGHTINGALE) THEN
               CALL schmidt_saltzman1993(sc_dms=sc_dms, sst=sst)
               CALL k_emis_n2000(kw_dms=kw_dms, sc_dms=sc_dms, w10=w10)

            ! Sea-air exchange coefficient for DMS from Liss and Merlivat (1986)
            ELSEIF(dms_opt == DMS_LM1986) THEN
               CALL schmidt_saltzman1993(sc_dms=sc_dms, sst=sst)
               CALL k_emis_lm1986(kw_dms=kw_dms, sc_dms=sc_dms, w10=w10)

            ! Sea-air exchange coefficient for DMS from Wanninkhof (2014)
            ELSEIF(dms_opt == DMS_W2014) THEN
               CALL schmidt_wanninkhof2014(sc_dms=sc_dms, sst=sst)
               CALL k_emis_w2014(kw_dms=kw_dms, sc_dms=sc_dms, w10=w10)

            ELSE
              !---- Other dms_opt cases are not supported
              CALL wrf_error_fatal ('dms_emis: dms_opt not supported')

            ENDIF ! dms_opt

            !---- Calculate DMS emissions
            !-- Calculate the DMS concentration in water at the interface
            ! Henry's Law coefficient for DMS (Sander2015 eq. 19 and Table 6)
            h_dms = 5.6E-3 * EXP(3500.* (1./sst - 1/298.15))
            press_dms = 1.E-6 * chem(i,kts,j,p_dms) * p_phy(i,kts,j)
            c_dms_interface =  press_dms * h_dms
            !-- Calculate the DMS emission flux
            ! In practice c_interface is often << c_water - could maybe
            ! neglect
            IF(c_dms_water > c_dms_interface) THEN
              ! e.g. Wanninkhof2014 eq. 1
              dms_emiss = kw_dms/100. * (c_dms_water-c_dms_interface)
            ELSE
              dms_emiss = 0.
            END IF
            ! Account for fractional sea ice (Loose2009 Fig. 3)
            IF(icefrac > 0.01) THEN
              dms_emiss = dms_emiss*(1.0-icefrac)**0.4
            ENDIF
            ! Convert to  WRF-Chem units (mol km-2 hour-1)
            dms_emiss = dms_emiss * 1.0E6

            !---- Add emissions to chem, at the surface
            IF(p_dms > param_first_scalar) THEN
              chem(i,kts,j,p_dms) = chem(i,kts,j,p_dms) + dms_emiss*conv
            ENDIF

            ! Accumulate DMS emissions for output (mol/m2)
            e_dms_ocean(i,j) = e_dms_ocean(i,j) &
                               + dms_emiss * 1.E-6 / 3600. * dtstep

          ENDIF ! c_dms_water
        END IF ! xland
      END DO ! j
    END DO ! i

  END SUBROUTINE dms_emis

!------------------------------------------------------------------------------

  SUBROUTINE schmidt_wanninkhof2014(sc_dms, sst)
    !
    ! Purpose:
    !   Calculate the schmidt coefficient of DMS based on Wanninkhof (2014)

    ! References:
    !   Wanninkhof (2014): https://doi.org/10.4319%2Flom.2014.12.351

    ! Arguments:
    !   Input:
    !     sst = Sea surface temperature (K)
    !   Output:
    !     sc_dms = Schmidt number for DMS (dimensionless)

    !---- Arguments
    REAL, INTENT(OUT) :: sc_dms
    REAL, INTENT(IN) :: sst

    !---- Local variables
    ! Sea surface temperature, corrected (°C)
    REAL :: sst_corr

    ! SST - Schmidt number parameterization valid -2 to 40°C
    sst_corr = sst - 273.15
    sst_corr = MAX(MIN(sst_corr, 40.0), -2.0)
    ! Schmidt number for DMS (Wanninkhof2014 Table 2)
    sc_dms = 2855.7 - 177.63*sst_corr + 6.0438*sst_corr**2.0 &
             - 0.11645*sst_corr**3. + 0.00094743*sst_corr**4.
    ! Different equation for fresh water, could be used for lakes
    ! sc_dms = 2595. - 163.12*sst_corr + 5.5902*sst_corr**2.0 &
    !          - 0.10817*sst_corr**3. + 0.00088204*sst_corr**4.

  END SUBROUTINE schmidt_wanninkhof2014

!------------------------------------------------------------------------------

  SUBROUTINE schmidt_saltzman1993(sc_dms, sst)
    !
    ! Purpose:
    !   Calculate the schmidt coefficient of DMS based on Saltzman et al., 1993

    ! References:
    !   Saltzman et al., 1993: https://doi.org/10.1029/93JC01858

    ! Arguments:
    !   Input:
    !     sst = Sea surface temperature (K)
    !   Output:
    !     sc_dms = Schmidt number for DMS (dimensionless)

    !---- Arguments
    REAL, INTENT(OUT) :: sc_dms
    REAL, INTENT(IN) :: sst

    !---- Local variables
    ! Sea surface temperature, corrected (°C)
    REAL :: sst_corr

    ! SST - Schmidt number parameterization valid 5 to 30°C
    sst_corr = sst - 273.15
    sst_corr = MAX(MIN(sst_corr, 30.0), 5.0)
    ! Schmidt number for DMS (Saltzman1993)
    sc_dms = 2674.0 - 147.12*sst_corr + 3.726*sst_corr**2.0 &
                - 0.038*sst_corr**3.

  END SUBROUTINE schmidt_saltzman1993

!------------------------------------------------------------------------------

  SUBROUTINE k_emis_w2014(kw_dms, sc_dms,  w10)
    !
    ! Purpose:
    !  Calculate DMS exchange coefficent from sea surface to the atmosphere
    !  using the parameterization of Wanninkhof (2014)

    ! References:
    !   Wanninkhof (2014): https://doi.org/10.4319%2Flom.2014.12.351

    ! Arguments:
    !   Input:
    !     sc_dms = Schmidt number for DMS (dimensionless)
    !     w10 = 10-m wind speed (m s-1)
    !   Output:
    !     kw_dms = Sea-air exchange coefficient for DMS (cm hour-1)

    IMPLICIT NONE

    !---- Arguments
    REAL, INTENT(OUT) :: kw_dms
    REAL, INTENT(IN) :: sc_dms
    REAL, INTENT(IN) :: w10

    ! Sea-air exchange coefficient (Wanninkhof2014 eq. 11)
    kw_dms = 0.251 * w10**2 * (sc_dms/660.)**(-0.5)

  END SUBROUTINE k_emis_w2014

!------------------------------------------------------------------------------

  SUBROUTINE k_emis_n2000(kw_dms, sc_dms,  w10)
    !
    ! Purpose:
    !  Calculate DMS exchange coefficent from sea surface to the atmosphere
    !  using the parameterization of Nightingale et al. (2000)

    ! References:
    !   Nightingale et al., 2000: https://doi.org/10.1029/1999GB900091

    ! Arguments:
    !   Input:
    !     sc_dms = Schmidt number for DMS (dimensionless)
    !     w10 = 10-m wind speed (m s-1)
    !   Output:
    !     kw_dms = Sea-air exchange coefficient for DMS (cm hour-1)

    IMPLICIT NONE

    !---- Arguments
    REAL, INTENT(OUT) :: kw_dms
    REAL, INTENT(IN) :: sc_dms
    REAL, INTENT(IN) :: w10

    ! Sea-air exchange coefficient, Nightingale2000 Sect. 5.2 and eq. 2
    kw_dms = (0.222*w10**2 + 0.333*w10) * (sc_dms/SC_CO2)**(-0.5)

  END SUBROUTINE k_emis_n2000

!------------------------------------------------------------------------------

  SUBROUTINE k_emis_lm1986(kw_dms, sc_dms, w10)
    !
    ! Purpose:
    !  Calculate DMS exchange coefficient from sea surface to the atmosphere
    !  using the parameterization of Liss and Merlivat (1986)

    ! References:
    !   Liss and Merlivat, 1986: https://doi.org/10.1007/978-94-009-4738-2_5

    ! Arguments:
    !   Input:
    !     sc_dms = Schmidt number for DMS (dimensionless)
    !     w10 = 10-m wind speed (m s-1)
    !   Output:
    !     kw_dms = Sea-air exchange coefficient for DMS (cm hour-1)

    IMPLICIT NONE

    !---- Arguments
    REAL, INTENT(OUT) :: kw_dms
    REAL, INTENT(IN) :: sc_dms
    REAL, INTENT(IN) :: w10

    ! Sea-air exchange coefficient (LissMerlivat1986 eqs. 7, 8, 9)
    IF(w10 <= 3.6) THEN
      kw_dms = 0.17*w10 * (sc_dms/SC_CO2)**(-0.66666666)
    ELSE IF(w10 <= 13.) THEN
      kw_dms = (2.85*w10 - 9.65) * (sc_dms/SC_CO2)**(-0.5)
    ELSE
      kw_dms = (5.9*w10 - 49.3) * (sc_dms/SC_CO2)**(-0.5)
    END IF

  END SUBROUTINE k_emis_lm1986

!------------------------------------------------------------------------------

END MODULE module_dms_emis
