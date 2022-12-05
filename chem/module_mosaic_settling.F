    module module_mosaic_settling

! 2018/04/19 Louis Marelle
! Calculates aerosol sedimentation velocities and performs the sedimentation for
! MOSAIC aerosols (mass, number and aerosol water). This routine only performs
! sedimentation  above the first level. Deposition by sedimentation in the first
! level is treated in mixactivate.
!
! For each i j, this routine calculates the column of deposition velocities for 
! each aerosol size, taking into account predicted aerosol wet radius (based on 
! water_* aerosols) and aerosol density. The method used is based on the Stokes 
! velocity + cunningham correction factor.
! The routine then updates aerosol mixing ratios at level k to account for the 
! loss of mass from settling in level k to level k-1 and increase of mass from 
! settling in level k+1 to level k (using an explicit method)
! In theory it would be better to perform sedimentation with a more complex
! advection scheme than a simple explicit method, since the explicit method
! introduces a lot of numerical downward diffusion for large particles.
!
! For now, this routine only takes into accound the settling in the 
! interstitial (*_a01 etc.) phase. Activated aerosols do not settle here, 
! although cloud droplet settling might already be treated in another module.
!
! The sedimentation is computed for all MOSAIC size bins, even though it's
! probably not significant for submicronic aerosols. This can certainly be
! optimized later.
!
! This code was developped for the following publication, also showing tests in
! the Arctic:
! Marelle, L., J.-C. Raut, K. S. Law, L. K. Berg, J. D. Fast, R. C. Easter, M.
! Shrivastava, and J. L. Thomas (2017), Improvements to the WRF-Chem 3.5.1 model
! for quasi- hemispheric simulations of aerosols and ozone in the Arctic,
! Geoscientific Model Devel- opment, 10(10), 3661–3677,
! doi:10.5194/gmd-10-3661-2017.
!
!LMarelle No sedimentation done for cw phase and hysteresis water. Should it be
!   done?
!LMarelle This would be more efficient (one less k loop) if the settling velocities
!   were calculated from top to bottom. The code is a bit more readable as is
!   however.


    contains

!-----------------------------------------------------------------------
    subroutine mosaic_settling_driver(                                 &
        t_phy, rho_phy, dz8w, chem, dtstep,                            &
        ids,ide, jds,jde, kds,kde,                                     &
        ims,ime, jms,jme, kms,kme,                                     &
        its,ite, jts,jte, kts,kte                                      )

! Import mosaic data (aerosol density by chem composition, molecular weights,
! bin diameters, bin volumes, etc.). dens_* is in (g/cm3) d*_sect etc. are 
! in (cm), vol_* in cm3
! Also import pointers to number, water and aerosol arrays
    USE module_data_mosaic_asect, ONLY:  &
        dens_water_aer, mw_water_aer, dens_aer, & 
        dlo_sect, dcen_sect, dhi_sect, volumlo_sect, &
        volumhi_sect, ncomp_aer, nsize_aer, ntype_aer, &
        numptr_aer, waterptr_aer, massptr_aer, lptr_so4_aer
    USE module_configure, ONLY:  num_chem
    USE module_data_mosaic_other, ONLY: pi 
    USE module_wrf_error !for wrf-debug
    implicit none

!   subroutine arguments
    integer, intent(in) ::   &
        ids, ide, jds, jde, kds, kde,   &
        ims, ime, jms, jme, kms, kme,   &
        its, ite, jts, jte, kts, kte

    real, intent(in),   &
        dimension( ims:ime, kms:kme, jms:jme ) :: &
        t_phy, rho_phy, dz8w

    real, intent(inout),   &
        dimension( ims:ime, kms:kme, jms:jme, 1:num_chem ) :: &
        chem

    real, intent(in) :: dtstep

!   local variables
    integer it, jt, kt, idiagaa, ntdt, ndt_settl, nsettl
    integer iphase, itype, n, ll, l1 
    integer, parameter :: r8=8
    real(r8) :: temp
    real(r8) :: airdens
    real(r8) :: dumfact, dummass, rnum
    real(r8) :: dtmax, aer_mass_above, aer_mass, aer_water_above, aer_water, aer_num
    real(r8) :: dt_settl, air_volume_vert_ratio
    real(r8) :: aer_mass_change, aer_vol_change, aer_num_change
    real(r8) :: wetdp, wetdens, wetmass, wetvol, &
         drydp, drydens, drymass, dryvol 
    real(r8) :: vsettl
    real(r8),  dimension(kts:kte) :: vsettl_col, dz_ij
!   diameters are in (cm), volumes are (cm^3/kg air) masses are (g/kg air)
!   parameters
    real(r8), parameter :: densdefault = 2.0 ! g/cm3
!    real, parameter :: smallmassaa = 1.0e-20 !g/mole air
!    real, parameter :: smallmassbb = 1.0e-30
    real(r8), parameter :: smallmassaa = 1.0e-19 ! g/kg air
    real(r8), parameter :: smallmassbb = 1.0e-29
    real(r8), parameter :: piover6 = pi/6.0
    real(r8), parameter :: onethird = 1.0/3.0
    CHARACTER (LEN=250) :: message

!   for mass checks
    real(r8),  dimension(kts:kte) :: check_mass_after
    real(r8),  dimension(kts:kte) :: check_mass_before
    real(r8) :: aer_mass_before


! ************************ Variables description *******************************
! *** Subr arguments ***
!    t_phy: 3D air temperature (K)
!    rho_phy: 3D air density (kg m-3)
!    dz8w: 3D grid height (m)
!    chem: advected chemical/aer species (ppm or µg kg-1)
!    dtstep: time step (s)
!
! *** Local variables ***
!    it, jt, kt, iphase, itype, n, ll, : loop indices for wrf i dimension, wrf j 
!        dimension, wrf k dimensuon, aerosol phase, aerosol type, aerosol size, 
!        aerosol chemical composition
! 
!    idiagaa : switch for diagnostic output
!
!    temp : air temperature (K)
!    airdens: air density (g cm-3)
!
!    ndt_settl : number of sub time loops with time step dt_settl 
!                to use if the model time step is too high to comply with the CFL condition
!    dt_settl : time step used in this routine to comply with the CFL condition.
!               In most cases, this will be equal to model dtstep (s)
! 
!    vsettl : settling velocity (cm s-1) at i,j,k,n
!    vsettl_col : settling velocity column (cm s-1) at i,j,n
!    dz_ij : height of vertical levels at i,j (cm ; column) 
!
!    aer_mass_above, aer_mass: aerosol mass mixing ratios (µg kg-1) at (i,j,k+1) and (i,j,k)
!    aer_water_above, aer_water: aerosol water mass mixing ratios (µg kg-1) at
!                                (i,j,k+1) and (i,j,k)
!    aer_num: aerosol number concentration (#/kg) at (i,j,k)
!    aer_mass_change, aer_vol_change, aer_num_change : mass, volume, number
!                     change in grid i,j,k from sedimentation
!    air_volume_vert_ratio : ratio [air_mass in grid (i,j,k+1)]/[air mass in grid i,j,k]
!  
!    rnum : aerosol number (# kg-1 air)
!    wetdp, wetdens, wetmass, wetvol: aerosol wet diameter (cm), density (g cm-3),
!                          mass (g kg-1 air), volume (cm-3 kg-1 air)
!    drydp, drydens, drymass, dryvol: same for dry aerosol
!    ***************************************************************************

!   loop over i and j
    do jt = jts, jte
    do it = its, ite

!   Get aerosol wet density and wet diameter

!   *** for now, just calc vdep_aer for iphase = 1
    iphase = 1

    do itype = 1, ntype_aer
    do n = 1, nsize_aer(itype)

    do kt = kts, kte

      dryvol = 0.0_r8
      drymass = 0.0_r8
  
      do ll = 1, ncomp_aer(itype) 
! This loop is "summing" chemical compounds; e.g. for bin 1, drymass = bc_a01+so4_a01+... etc and 
! dryvol = vol_bc_a01+... etc.  to compute bin diameters ((dryvol/(piover6*rnum))**onethird) 
          l1 = massptr_aer(ll,n,itype,iphase)
          dummass = chem(it,kt,jt,l1) * 1.0E-6_r8 ! chem is in ug/kg, convert to g/kg air
          drymass = drymass + dummass
          dryvol = dryvol + dummass/dens_aer(ll,itype) ! cm3/kg air
      end do
  
      l1 = waterptr_aer(n,itype)
      dummass = chem(it,kt,jt,l1) * 1.0E-6_r8 ! chem is in ug/kg, convert to g/kg air
      wetmass = drymass + dummass
      wetvol = dryvol + dummass/dens_water_aer ! cm3/kg air, dens is in g/cm3
  
      l1 = numptr_aer(n,itype,iphase)
      rnum = chem(it,kt,jt,l1) ! #/kg
  
      if (drymass .le. smallmassbb) then
          drydp = dcen_sect(n,itype)
          drydens = densdefault
          wetdp = drydp
          wetdens = drydens
  
      else ! (drymass > smallmassbb)
          if (drymass .le. smallmassaa) then 
              wetmass = drymass
              wetvol = dryvol
          end if
          drydens = drymass/dryvol ! g/cm3
          wetdens = wetmass/wetvol ! g/cm3
  
          if (rnum .ge. dryvol/volumlo_sect(n,itype)) then 
          ! This is capping the dp diameters to [dlo - dhi],
          ! This routine is estimating one sedimentation velocity per bin
          ! (assuming internal mixing in the bin)
              drydp = dlo_sect(n,itype)
          else if (rnum .le. dryvol/volumhi_sect(n,itype)) then 
              drydp = dhi_sect(n,itype) 
          else
              drydp = (dryvol/(piover6*rnum))**onethird
          end if
  
          if(abs(wetvol).gt.(1000.*abs(dryvol))) then
              dumfact=10.0
          else
              dumfact=abs(wetvol/dryvol)**onethird
              dumfact=max(1.0,min(dumfact,10.0))
          endif
          wetdp = drydp*dumfact
  
      endif !if (drymass .le. smallmassbb)
  
!   Get air temperature & density at it, kt, jt
      temp = t_phy(it,kt,jt)
      airdens = rho_phy(it,kt,jt)*1.0e-3_r8 !kg/m3 -> g/cm3
  
!   Calculate settling velocity for mass
      call aerosol_depvel(   &
           wetdp, wetdens,temp, airdens, vsettl)
  
!     Keep the vsettl column
      vsettl_col(kt)=vsettl

    end do !kt

! Sedimentation

!   vertical level depths in cm at i,j
    dz_ij(kts:kte) = dz8w(it,kts:kte,jt)*1.0E2_r8 !m -> cm
!   Compute the mass column for mass checks, only for one species for now
    l1 = lptr_so4_aer(n,itype,iphase)
    !µg/kg -> µg/m3 -> µg/m2
    check_mass_before(kts:kte)=(dz8w(it,kts:kte,jt)*chem(it,kts:kte,jt,l1)*rho_phy(it,kts:kte,jt))

    do kt = kts, kte-1

      ! Determine the maximum time-step satisying the CFL condition:
      ! dt <= (dz)_min / v_settl
      ntdt=INT(dtstep)
      dtmax = dz_ij(kt+1) / vsettl_col(kt+1) 
      dtmax = max( dtmax, dz_ij(kt) / vsettl_col(kt) ) 
      !LMarelle It might be simpler to just use one dtmax for the whole column like in settl vash
      ndt_settl = MAX( 1, INT( ntdt /dtmax) )
      dt_settl = REAL(ntdt) / REAL(ndt_settl)
      ! limit maximum number of iterations to 12. In this case, settling will be
      ! underestimated (this can only happen with very big time steps, big
      ! particles and very thin vertical levels) !LMarelle add a warning, or
      ! call to wrf debug
      IF (ndt_settl > 12) ndt_settl = 12

      ! This ratio takes into account the change of air mass in grids
      ! between levels kt+1 and kt. If this is not taken into account, mass will
      ! not be conserved during settling (we need to transfer mass, not mixing
      ! ratios, between levels kt+1 and kt)
      air_volume_vert_ratio = ( rho_phy(it,kt+1,jt) * dz_ij(kt+1) ) / (rho_phy(it,kt,jt) * dz_ij(kt) )

      ! aerosol number at kt
      l1 = numptr_aer(n,itype,iphase)
      aer_num = chem(it,kt,jt,l1)

      ! loop on aerosol chemical composition (does not contain aerosol water)
      do ll = 1, ncomp_aer(itype)

        ! aerosol mass at kt and kt+1 for species ll
        l1 = massptr_aer(ll,n,itype,iphase)
        aer_mass = chem(it,kt,jt,l1)
        aer_mass_above = chem(it,kt+1,jt,l1)
        aer_mass_before = aer_mass

!----------------------------------------------------------------------
! Do sedimentation for aerosol mass and number using an explicit scheme
!----------------------------------------------------------------------
        do nsettl = 1, ndt_settl
           ! Surface case, no loss to the level below (deposition to the surface
           ! is already done elsewhere):
           if (kt .eq. kts) then 
             aer_mass =  aer_mass + aer_mass_above * ( vsettl_col(kt+1) * dt_settl/dz_ij(kt+1) ) * &
                         air_volume_vert_ratio 
             aer_mass_change = aer_mass_above * ( vsettl_col(kt+1) * dt_settl/dz_ij(kt+1) ) * & 
                               air_volume_vert_ratio * 1.0E-6_r8 ! 1E-6 = µg/kg to g/kg
             aer_vol_change = aer_mass_change / dens_aer(ll,itype) ! cm3/kg
             ! LMarelle: Calculate the number change using the center diameter,
             ! to be consistent with the number concentration update in
             ! mosaic_addemiss
             aer_num_change = 6.0_r8 * aer_vol_change * 1.0_r8/pi * 1.0_r8 / dcen_sect(n,itype) ** 3.0_r8  ! #/kg 
             aer_num = aer_num + aer_num_change
           ! Top level case, no mass gain from level k+1:
           elseif (kt .eq. kte) then 
             aer_mass =  aer_mass * (1.0_r8 - vsettl_col(kt)*dt_settl/dz_ij(kt)) 
             aer_mass_change = - aer_mass * vsettl_col(kt)*dt_settl/dz_ij(kt) * 1.0E-6_r8
             aer_vol_change = aer_mass_change / dens_aer(ll,itype)
             aer_num_change = 6.0_r8 * aer_vol_change * 1.0_r8/pi * 1.0_r8 / dcen_sect(n,itype) ** 3.0_r8
             aer_num = aer_num + aer_num_change
           ! Default case, aerosols enter from level k+1 and exit to level k-1:
           else 
             aer_mass =  aer_mass * ( 1.0_r8 - vsettl_col(kt)*dt_settl / dz_ij(kt) ) + &
                         aer_mass_above * ( vsettl_col(kt+1) * dt_settl / dz_ij(kt+1) ) * &
                         air_volume_vert_ratio
             aer_mass_change = (- aer_mass * vsettl_col(kt) * dt_settl / dz_ij(kt) + & 
                               aer_mass_above * ( vsettl_col(kt+1) * dt_settl / dz_ij(kt+1) ) * &
                               air_volume_vert_ratio) * 1.0E-6_r8
             aer_vol_change = aer_mass_change / dens_aer(ll,itype)
             aer_num_change = 6.0_r8 * aer_vol_change * 1.0_r8/pi * 1.0_r8 / dcen_sect(n,itype) ** 3.0_r8
             aer_num = aer_num + aer_num_change
           end if
             
        end do !dtsettl 

        ! Diagnostic output if needed
        idiagaa = 0
        if (idiagaa>0) print 9310, it, jt, kt, n, ll,   &
             vsettl_col(kt), dt_settl, dz_ij(kt), aer_mass, aer_mass_above, aer_mass_change, aer_mass-aer_mass_before
             9310  format( 'aerdep', 4i4, i3, 1p, 7e10.2,   &
             2x, 0p )

        ! Update aerosol mass mixing ratios
        if (aer_mass < 0.0_r8 ) aer_mass = 0.0_r8
        l1 = massptr_aer(ll,n,itype,iphase) 
        chem(it,kt,jt,l1) = aer_mass
 
      end do !ll, aer chem comp

      ! Update number concentrations, done here because it is calculated based
      ! on the mass change for all species, calculated in the ll loop above
      if (aer_num < 0.0_r8) aer_num = 0.0_r8
      l1 = numptr_aer(n,itype,iphase)
      chem(it,kt,jt,l1) = aer_num

      ! Now do sedimentation for aerosol water (LMarelle: using a subroutine for this and
      ! the mass change would be cleaner) 
      l1 = waterptr_aer(n,itype)
      aer_water = chem(it,kt,jt,l1)
      aer_water_above = chem(it,kt+1,jt,l1)
      do nsettl = 1, ndt_settl
        if (kt == kts) then ! surface
          aer_water =  aer_water + aer_water_above * ( vsettl_col(kt+1) * dt_settl/dz_ij(kt+1) ) * &
                       air_volume_vert_ratio
        elseif (kt == kte) then ! top
          aer_water =  aer_water * (1.0_r8 - vsettl_col(kt)*dt_settl/dz_ij(kt))
        else ! elsewhere in the column
          aer_water =  aer_water * ( 1.0_r8 - vsettl_col(kt)*dt_settl / dz_ij(kt) )+ &
                       aer_water_above * ( vsettl_col(kt+1) * dt_settl /dz_ij(kt+1) ) * &
                       air_volume_vert_ratio
        end if
      end do ! do nsettl
      ! Update the aerosol water mixing ratio
      if (aer_water < 0.0_r8) aer_water = 0.0_r8
      l1 = waterptr_aer(n,itype)
      chem(it,kt,jt,l1)=aer_water

    end do ! do kt

!   Compute the mass column for mass checks, only for sulfate
    l1 = lptr_so4_aer(n,itype,iphase)
    !µg/kg -> µg/m3 -> µg/m2
    check_mass_after(kts:kte)=dz8w(it,kts:kte,jt)*chem(it,kts:kte,jt,l1)*rho_phy(it,kts:kte,jt)
!   Diagnostic output
    WRITE(message, * ) ' aer_settling it,jt,masscol_before,masscol_after = ',it, jt, sum(check_mass_before), sum(check_mass_after)
    IF (idiagaa >0) CALL wrf_debug ( 15, message )

    end do !n, size
    end do !itype, type

    end do !it
    end do !jt


    end subroutine mosaic_settling_driver




!------------------------------------------------------------------------
    subroutine aerosol_depvel(   &
              dgnum, aerodens,   &
              temp, airdens,            &
              vsettl )

! 2018/10/19, Louis Marelle: This routine is based on the subroutine in
! module_mosaic_drydep.F
!
!   Computes the gravitational settling velocity (stokes velocity + cunningham
!   correction factor)
!
!   input parameters
!    dgnum - geometric mean diameter for aerosol number (cm)
!    aerodens - aerosol density (dgnum and aerodens are for the
!        actual wet distribution)
!    temp - temperature (K)
!    airdens - air density (g/cm^3)
!   output parameters
!    vsettl - gravitational settling velocity (cm/s)
!    
    implicit none

    integer, parameter :: r8=8
!   subroutine arguments
    real(r8), intent(in) :: temp, airdens
    real(r8), intent(in) :: dgnum, aerodens
    real(r8), intent(out) :: vsettl

!   local variables
!    airkinvisc - air kinematic viscosity (cm^2/s)
    real airkinvisc
!   air molecular freepath (cm)
    real freepath
    real vsettl_dgnum
!   Cunningham factor (adimensional)
    real cunningham_fact
    real pi
    parameter (pi = 3.1415926536)
!   gravity = gravitational acceleration in cm/s^2
    real gravity
    parameter (gravity = 980.616)
    
!   Calculate the kinetic viscosity and molecular freepath of air
    airkinvisc = ( 1.8325e-4 * (416.16/(temp+120.0)) *   &
                      ((temp/296.16)**1.5) ) / airdens
    freepath = 7.39758e-4 * airkinvisc / sqrt(temp)

!   Cunningham correction factor for the stokes velocity
    cunningham_fact = 1.0 + 2.0 * freepath/dgnum * ( 1.257+ 0.400 * exp( -0.55 * dgnum/freepath ) )

!   Calculate the sedimentation (Stokes) velocity
    vsettl_dgnum = (gravity*aerodens*dgnum*dgnum)/   &
              (18.*airkinvisc*airdens)

!   Apply the correction factor
    vsettl = vsettl_dgnum * cunningham_fact

    return
    end subroutine aerosol_depvel

!------------------------------------------------------------------------
    end module module_mosaic_settling

