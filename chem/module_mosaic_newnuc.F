!**********************************************************************************  
! This computer software was prepared by Battelle Memorial Institute, hereinafter
! the Contractor, under Contract No. DE-AC05-76RL0 1830 with the Department of 
! Energy (DOE). NEITHER THE GOVERNMENT NOR THE CONTRACTOR MAKES ANY WARRANTY,
! EXPRESS OR IMPLIED, OR ASSUMES ANY LIABILITY FOR THE USE OF THIS SOFTWARE.
!
! MOSAIC module: see module_mosaic_driver.F for references and terms of use
!**********************************************************************************  
	module module_mosaic_newnuc



	use module_peg_util



	implicit none



	contains



!-----------------------------------------------------------------------
	subroutine mosaic_newnuc_1clm( istat_newnuc,   &
		it, jt, kclm_calcbgn, kclm_calcend,   &
		idiagbb_in,   &
		dtchem, dtnuc_in, rsub0,   &
		id, ktau, ktauc, its, ite, jts, jte, kts, kte, newnuc_method, nuc_msa_opt, nuc_msa_fac)
!
!   calculates new particle nucleation for grid points in
!	the i=it, j=jt vertical column over timestep dtnuc_in
!   works with mosaic sectional aerosol packages
!
!   when newnuc_method = 1 (internal parameter), uses
!	h2so4-nh3-h2o ternary nucleation from napari et al., 2002
!       *** this option is currently disabled because the ternary
!           parameterization was recently shown to be invalid
!   when newnuc_method = 2 (internal parameter), uses
!	h2so4-h2o binary nucleation from wexler et al., 1994
!
!       R. Lapere 2025
!       newnuc_method is now a namelist option corresponding to 4 possible schemes:
!       1 = Napari et al., 2002 - currently disabled because invalid
!       2 = Wexler et al., 1994 
!       NEW 2025
!       3 = Kulmala et al., 1998
!       4 = Vehkamaki et al., 2002
!		
!	Note: the size of new particles is the lower-bound size of the smallest size bin, 
!	i.e. the mass of the new particles is calculated assuming they have the size of the smallest bin, 
!	leading to overestimations in depletion of H2SO4, which will feedback to calculated rates. 
!	Loss of particles to coagulation between cluster and smallest bin is also not accounted for, 
!	leading to overestimates of particle number in the smallest bin.


	use module_data_mosaic_asect
	use module_data_mosaic_other
	use module_state_description, only:  param_first_scalar

!   subr arguments
	integer, intent(inout) :: istat_newnuc    ! =0 if no problems
	integer, intent(in) ::   &
		it, jt, kclm_calcbgn, kclm_calcend,   &
		idiagbb_in,   &
		id, ktau, ktauc, its, ite, jts, jte, kts, kte
	real, intent(in) :: dtchem, dtnuc_in
	real, intent(in) :: rsub0(l2maxd,kmaxd,nsubareamaxd)
	integer, intent(in) :: nuc_msa_opt ! toggles MSA nucleation on/off - namelist option
	real, intent(in) :: nuc_msa_fac ! tuning factor applied to nucleation rate if nuc_msa_opt = 1

!       rsub0 holds mixrat values before the aerchemistry calcs

!   NOTE - much information is passed via arrays in 
!		module_data_mosaic_asect and module_data_mosaic_other
!
!   rsub (inout) - trace gas and aerosol mixing ratios
!   aqmassdry_sub, aqvoldry_sub (inout) - 
!		aerosol dry-mass and dry-volume mixing ratios
!   adrydens_sub (inout) - aerosol dry density
!   rsub(ktemp,:,:), cairclm, relhumclm (in) - 
!		air temperature, molar density, and relative humidity

	integer, intent(in) :: newnuc_method

	integer :: k, l, ll, m
	integer :: isize, itype, iphase
	integer :: iconform_numb
	integer :: idiagbb
	integer :: nsize, ntau_nuc
	integer, save :: ncount(10)
	integer :: p1st

	real, parameter :: densdefault = 2.0
	real, parameter :: smallmassbb = 1.0e-30

! nh4hso4 values for a_zsr and b_zsr
    real, parameter :: a_zsr_xx1 =  1.15510
    real, parameter :: a_zsr_xx2 = -3.20815
    real, parameter :: a_zsr_xx3 =  2.71141
    real, parameter :: a_zsr_xx4 =  2.01155
    real, parameter :: a_zsr_xx5 = -4.71014
    real, parameter :: a_zsr_xx6 =  2.04616
    real, parameter :: b_zsr_xx  = 29.4779

	real :: aw
	real :: cair_box
	real :: dens_nh4so4a, dtnuc
	real :: duma, dumb, dumc, duma_msa ! dummy variables for intermediate computations
	real :: rh_box
	real :: qh2so4_avg, qh2so4_cur, qh2so4_del 
	real :: qnh3_avg, qnh3_cur, qnh3_del 
	real :: qnuma_del, qso4a_del, qnh4a_del
	real :: temp_box, qv_box
	real :: xxdens, xxmass, xxnumb, xxvolu

    ! mixing ratio of gas-phase msa - mole/mole_air
    ! cur is current, avg is average between cur and rsub0, del is the variation due to nucleation
	real :: qmsa_avg, qmsa_cur, qmsa_del
	real :: qmsaa_del ! change in aerosol msa mixing ratio - mole/mole_air
	real, parameter :: dens_msaa = 1.48 ! density of msa aerosol - g/cm3
	real :: j2_msa ! msa nucleation rate - #/cm3

	real,save :: dumveca(10), dumvecb(10), dumvecc(10), dumvecd(10), dumvece(10)
	real :: volumlo_nuc(maxd_asize), volumhi_nuc(maxd_asize)

	character(len=100) :: msg

    p1st = PARAM_FIRST_SCALAR

!   check newnuc_method
	istat_newnuc = 0
	if (newnuc_method .le. 1) then
	    if ((it .eq. its) .and. (jt .eq. jts))   &
		call peg_message( lunerr,   &
		'*** mosaic_newnuc_1clm -- illegal newnuc_method' )
	    istat_newnuc = -1
	    return
	end if


!   when dtnuc_in > dtchem, do not perform nucleation calcs 
!   on every chemistry time step
	ntau_nuc = nint( dtnuc_in/dtchem )
	ntau_nuc = max( 1, ntau_nuc )
	if (mod(ktau,ntau_nuc) .ne. 0) return
	dtnuc = dtchem*ntau_nuc


!   set variables that do not change
	idiagbb = idiagbb_in

	itype = 1
	iphase = ai_phase
	nsize = nsize_aer(itype)
	volumlo_nuc(1:nsize) = volumlo_sect(1:nsize,itype)
	volumhi_nuc(1:nsize) = volumhi_sect(1:nsize,itype)


!   loop over subareas (currently only 1) and vertical levels
	do 2900 m = 1, nsubareas

	do 2800 k = kclm_calcbgn, kclm_calcend


!   initialize diagnostics
	if ((it .eq. its) .and.   &
	    (jt .eq. jts) .and. (k .eq. kclm_calcbgn)) then
	    dumveca(:) = 0.0         ! current grid param values
	    dumvecb(:) = +1.0e35     ! param minimums
	    dumvecc(:) = -1.0e35     ! param maximums
	    dumvecd(:) = 0.0         ! param averages
	    dumvece(:) = 0.0         ! param values for highest qnuma_del
	    ncount(:) = 0
	end if


	ncount(1) = ncount(1) + 1
	if (afracsubarea(k,m) .lt. 1.e-4) goto 2700

	cair_box = cairclm(k)
	temp_box = rsub(ktemp,k,m)
	rh_box = relhumclm(k)
	qv_box = rsub(kh2o,k,m)

	qh2so4_cur = max(0.0,rsub(kh2so4,k,m))
	qnh3_cur   = max(0.0,rsub(knh3,k,m))
	qh2so4_avg = 0.5*( qh2so4_cur + max(0.0,rsub0(kh2so4,k,m)) )
	qnh3_avg   = 0.5*( qnh3_cur   + max(0.0,rsub0(knh3,k,m)) )

	qh2so4_del = 0.0
	qnh3_del = 0.0
	qnuma_del = 0.0
	qso4a_del = 0.0
	qnh4a_del = 0.0

	if (nuc_msa_opt .eq. 1 .and. kmsa .ge. p1st) then
	  qmsa_cur = max(0.0,rsub(kmsa,k,m))
	  qmsa_avg   = 0.5*( qmsa_cur   + max(0.0,rsub0(kmsa,k,m)) )
	  qmsa_del = 0.0
	  j2_msa = 0.0
	endif

	dens_nh4so4a = dens_so4_aer

	isize = 0

!   make call to nucleation routine
	if (newnuc_method .eq. 1) then
            call ternary_nuc_mosaic_1box(   &
               dtnuc, temp_box, rh_box, cair_box,   &
               qh2so4_avg, qh2so4_cur, qnh3_avg, qnh3_cur,   &
               nsize, maxd_asize, volumlo_nuc, volumhi_nuc,   &
               isize, qnuma_del, qso4a_del, qnh4a_del,   &
               qh2so4_del, qnh3_del, dens_nh4so4a )
	else if (newnuc_method .eq. 2) then
            call wexler_nuc_mosaic_1box(   &
               dtnuc, temp_box, rh_box, cair_box,   &
               qh2so4_avg, qh2so4_cur, qnh3_avg, qnh3_cur,   &
               nsize, maxd_asize, volumlo_nuc, volumhi_nuc,   &
               isize, qnuma_del, qso4a_del, qnh4a_del,   &
               qh2so4_del, qnh3_del, dens_nh4so4a )
    else if (newnuc_method .eq. 3) then
	    	call kulmala1998(   &
	       	   dtnuc, temp_box, rh_box, cair_box, qv_box,  &
               qh2so4_avg, qh2so4_cur, &
               isize, maxd_asize, volumlo_sect, volumhi_sect,   &
               qnuma_del, qh2so4_del, qso4a_del, dens_so4_aer)
    else if (newnuc_method .eq. 4) then
            call vehkamaki2002(   &
               dtnuc, temp_box, rh_box, cair_box,  &
               qh2so4_avg, qh2so4_cur, &
               isize, maxd_asize, volumlo_sect, volumhi_sect,   &
               qnuma_del, qh2so4_del, qso4a_del, dens_so4_aer)
	else
	    istat_newnuc = -1
	    return
	end if

    ! if namelist option for msa nucleation is activated
    ! compute nucleation rate following Riccobono et al. 2014
	! and modify concentrations of MSA gas, MSA aerosol and aerosol number accordingly
	if (nuc_msa_opt .eq. 1 .and. kmsa .ge. p1st) then
	   call riccobono_nuc_msa(   &
	      dtnuc, temp_box, rh_box, cair_box, qh2so4_avg, qh2so4_cur, &
	      qnh3_avg, maxd_asize, volumlo_nuc, volumhi_nuc, &
	      qnuma_del, qh2so4_del, qmsa_avg, qmsa_cur, qmsa_del, qmsaa_del, &
	      j2_msa, nuc_msa_fac )
	end if


	

!   temporary diagnostics
	dumveca(1) = temp_box
	dumveca(2) = rh_box
	dumveca(3) = rsub(kso2,k,m)
	dumveca(4) = qh2so4_avg
	dumveca(5) = qnh3_avg
	dumveca(6) = qnuma_del
	do l = 1, 6
	    dumvecb(l) = min( dumvecb(l), dumveca(l) )
	    dumvecc(l) = max( dumvecc(l), dumveca(l) )
	    dumvecd(l) = dumvecd(l) + dumveca(l)
	    if (qnuma_del .gt. dumvece(6)) dumvece(l) = dumveca(l)
	end do


!   check for zero new particles
	if (qnuma_del .le. 0.0) goto 2700

!   check for valid isize
	if (isize .ne. 1) ncount(3) = ncount(3) + 1
	if ((isize .lt. 1) .or. (isize .gt. nsize)) then
	    write(msg,93010) 'newnucxx bad isize_nuc' , it, jt, k,   &
		isize, nsize
	    call peg_message( lunerr, msg )
	    goto 2700
	end if
93010	format( a, 3i3, 1p, 9e10.2 )


	ncount(2) = ncount(2) + 1

!   update gas and aerosol so4, nh3/nh4 and msa mixing ratios
	rsub(kh2so4,k,m) = max( 0.0, rsub(kh2so4,k,m) + qh2so4_del )
	rsub(knh3,  k,m) = max( 0.0, rsub(knh3,  k,m) + qnh3_del )
	if (nuc_msa_opt .eq. 1 .and. kmsa .ge. p1st) then
	  rsub(kmsa,  k,m) = max( 0.0, rsub(kmsa,  k,m) + qmsa_del )
	endif
	
	l = lptr_so4_aer(isize,itype,iphase)
	if (l .ge. p1st) then
	    rsub(l,k,m) = rsub(l,k,m) + qso4a_del
	end if
	l = lptr_nh4_aer(isize,itype,iphase)
	if (l .ge. p1st) then
	    rsub(l,k,m) = rsub(l,k,m) + qnh4a_del
	end if

	if (nuc_msa_opt .eq. 1 .and. kmsa .ge. p1st) then
	  l = lptr_msa_aer(isize,itype,iphase)
	  if (l .ge. p1st) then
	    rsub(l,k,m) = rsub(l,k,m) + qmsaa_del
	  end if
	endif

	l = numptr_aer(isize,itype,iphase)
	rsub(l,k,m) = rsub(l,k,m) + qnuma_del
	xxnumb = rsub(l,k,m)

!   update aerosol water, using mosaic parameterizations for nh4hso4
!     duma = (mole-salt)/(mole-salt+water)
!     dumb = (mole-salt)/(kg-water)
!     dumc = (mole-water)/(mole-salt)
	l = waterptr_aer(isize,itype)
	if ((rh_box .gt. 0.10) .and. (l .ge. p1st)) then
	    aw = min( rh_box, 0.98 )
	    if (aw .lt. 0.97) then
		duma =       a_zsr_xx1 +   &
		        aw*( a_zsr_xx2 +   &
		        aw*( a_zsr_xx3 +   &
		        aw*( a_zsr_xx4 +   &
		        aw*( a_zsr_xx5 +   &
		        aw*  a_zsr_xx6 ))))
	    else
		dumb = -b_zsr_xx*log(aw)
	        dumb = max( dumb, 0.5 )
		duma = 1.0/(1.0 + 55.509/dumb)
	    end if
	    duma = max( duma, 0.01 )
	    dumc = (1.0 - duma)/duma
	    rsub(l,k,m) = rsub(l,k,m) + qso4a_del*dumc
	end if


!
!   update dry mass, density, and volume,
!   and check for mean dry-size within bounds
!
	xxmass = aqmassdry_sub(isize,itype,k,m)
	xxdens = adrydens_sub( isize,itype,k,m)
	iconform_numb = 1

	if ((xxdens .lt. 0.1) .or. (xxdens .gt. 20.0)) then
!   (exception) case of drydensity not valid
	    continue
	else
!   (normal) case of drydensity valid (which means drymass is valid also)
!   so increment mass and volume with the so4, nh4 and msa deltas, then calc density
	    xxvolu = xxmass/xxdens
	    duma = qso4a_del*mw_so4_aer + qnh4a_del*mw_nh4_aer
	    if (nuc_msa_opt .eq. 1 .and. kmsa .ge. p1st) then
	      duma_msa = qmsaa_del*mw_msa_aer
	    else
	      duma_msa = 0.0
	    endif
	    xxmass = xxmass + duma + duma_msa
	    xxvolu  = xxvolu  + duma/dens_nh4so4a + duma_msa/dens_msaa
	    if (xxmass .le. smallmassbb) then
!   do this to force calc of dry mass, volume from rsub
		xxdens = 0.001
	    else if (xxmass .gt. 1000.0*xxvolu) then
!   in this case, density is too large.  setting density=1000 forces
!   next IF block while avoiding potential divide by zero or overflow
		xxdens = 1000.0
	    else 
		xxdens = xxmass/xxvolu
	    end if
	end if

	if ((xxdens .lt. 0.1) .or. (xxdens .gt. 20.0)) then
!   (exception) case of drydensity not valid (or drymass extremely small), 
!   so compute from dry mass, volume from rsub
	    ncount(4) = ncount(4) + 1
	    xxmass = 0.0
	    xxvolu  = 0.0
	    do ll = 1, ncomp_aer(itype)
		l = massptr_aer(ll,isize,itype,iphase)
		if (l .ge. p1st) then
		    duma = max( 0.0, rsub(l,k,m) )*mw_aer(ll,itype)
		    xxmass = xxmass + duma + duma_msa
		    xxvolu = xxvolu + duma/dens_aer(ll,itype)
		end if
	    end do
	end if

	if (xxmass .le. smallmassbb) then
!   when drymass extremely small, use default density and bin center size,
!   and zero out water
	    ncount(5) = ncount(5) + 1
	    xxdens = densdefault
	    xxvolu = xxmass/xxdens
	    xxnumb = xxmass/(volumcen_sect(isize,itype)*xxdens)
	    iconform_numb = 0
	    l = waterptr_aer(isize,itype)
	    if (l .ge. p1st) rsub(l,k,m) = 0.0
	    l = hyswptr_aer(isize,itype)
	    if (l .ge. p1st) rsub(l,k,m) = 0.0
	else
	    xxdens = xxmass/xxvolu
	end if

	if (iconform_numb .gt. 0) then
!   check for mean dry-size within bounds, and conform number if not
	    if (xxnumb .gt. xxvolu/volumlo_sect(isize,itype)) then
		ncount(6) = ncount(6) + 1
		xxnumb = xxvolu/volumlo_sect(isize,itype)
	    else if (xxnumb .lt. xxvolu/volumhi_sect(isize,itype)) then
		ncount(7) = ncount(7) + 1
		xxnumb = xxvolu/volumhi_sect(isize,itype)
	    end if
	end if

!   load dry mass, density, volume, and (possibly conformed) number
	l = numptr_aer(isize,itype,iphase)
	rsub(l,k,m) = xxnumb
	adrydens_sub( isize,itype,k,m) = xxdens
	aqmassdry_sub(isize,itype,k,m) = xxmass
	aqvoldry_sub( isize,itype,k,m) = xxvolu


2700	continue

!   temporary diagnostics
	if ((idiagbb .ge. 100) .and.   &
	    (it .eq. ite) .and.    &
	    (jt .eq. jte) .and. (k .eq. kclm_calcend)) then
	    if (idiagbb .ge. 110) then
	      write(msg,93020) 'newnucbb mins ', dumvecb(1:6)
	      call peg_message( lunerr, msg )
	      write(msg,93020) 'newnucbb maxs ', dumvecc(1:6)
	      call peg_message( lunerr, msg )
	      duma = max( 1, ncount(1) ) 
	      write(msg,93020) 'newnucbb avgs ', dumvecd(1:6)/duma
	      call peg_message( lunerr, msg )
	      write(msg,93020) 'newnucbb hinuc', dumvece(1:6)
	      call peg_message( lunerr, msg )
	      write(msg,93020) 'newnucbb dtnuc', dtnuc
	      call peg_message( lunerr, msg )
	    end if
	    write(msg,93030) 'newnucbb ncnt ', ncount(1:7)
	    call peg_message( lunerr, msg )
	end if
93020	format( a, 1p, 10e10.2 )
93030	format( a, 1p, 10i10 )


2800	continue	! k levels

2900	continue	! subareas


	return
	end subroutine mosaic_newnuc_1clm


!-----------------------------------------------------------------------
! note - subrs ternary_nuc_mosaic_1box, ternary_nuc_napari, wexler_nuc_mosaic_1box
!    taken from ternucl04.f90 on 22-jul-2006
!-----------------------------------------------------------------------
        subroutine ternary_nuc_mosaic_1box(   &
           dtnuc, temp_in, rh_in, cair,   &
           qh2so4_avg, qh2so4_cur, qnh3_avg, qnh3_cur,   &
           nsize, maxd_asize, volumlo_sect, volumhi_sect,   &
           isize_nuc, qnuma_del, qso4a_del, qnh4a_del,   &
           qh2so4_del, qnh3_del, dens_nh4so4a )
!.......................................................................
!
! calculates new particle production from h2so4-nh3-h2o ternary nucleation
!    over timestep dtnuc, using nucleation rates from the
!    napari et al. (2002) parameterization
!
! the new particles are "grown" to the lower-bound size of the host code's 
!    smallest size bin.  (this "growth" is purely ad hoc, and would not be
!    necessary if the host code's size bins extend down to ~1 nm.)
!    if the h2so4 and nh3 mass mixing ratios (mixrats) of the grown new 
!    particles exceed the current gas mixrats, the new particle production
!    is reduced so that the new particle mass mixrats match the gas mixrats.
!
! revision history
!    coded by rc easter, pnnl, 20-mar-2006
!
! key routines called: subr ternary_nuc_napari
!
! references:
!    napari, i., m. noppel, h. vehkamaki, and m. kulmala,
!       parameterization of ternary nucleation rates for
!       h2so4-nh3-h2o vapors.  j. geophys. res., 107, 4381, 2002.
!
!.......................................................................
        implicit none

! subr arguments (in)
        real, intent(in) :: dtnuc             ! nucleation time step (s)
        real, intent(in) :: temp_in           ! temperature, in k
        real, intent(in) :: rh_in             ! relative humidity, as fraction
        real, intent(in) :: cair              ! dry-air molar density (mole/cm3)

        real, intent(in) :: qh2so4_avg, qh2so4_cur   ! gas h2so4 mixing ratios (mole/mole-air)
        real, intent(in) :: qnh3_avg, qnh3_cur       ! gas nh3 mixing ratios (mole/mole-air)
             ! qxxx_cur = current value (at end of condensation)
             ! qxxx_avg = average value (from start to end of condensation)

        integer, intent(in) :: nsize                    ! number of aerosol size bins
        integer, intent(in) :: maxd_asize               ! dimension for volumlo_sect, ...
        real, intent(in) :: volumlo_sect(maxd_asize)    ! dry volume at lower bnd of bin (cm3)
        real, intent(in) :: volumhi_sect(maxd_asize)    ! dry volume at upper bnd of bin (cm3)

! subr arguments (out)
        integer, intent(out) :: isize_nuc     ! size bin into which new particles go
        real, intent(out) :: qnuma_del        ! change to aerosol number mixing ratio (#/mole-air)
        real, intent(out) :: qso4a_del        ! change to aerosol so4 mixing ratio (mole/mole-air)
        real, intent(out) :: qnh4a_del        ! change to aerosol nh4 mixing ratio (mole/mole-air)
        real, intent(out) :: qh2so4_del       ! change to gas h2so4 mixing ratio (mole/mole-air)
        real, intent(out) :: qnh3_del         ! change to gas nh3 mixing ratio (mole/mole-air)
                                              ! aerosol changes are > 0; gas changes are < 0

! subr arguments (inout)
        real, intent(inout) :: dens_nh4so4a   ! dry-density of the new nh4-so4 aerosol mass (g/cm3)
                                              ! use 'in' value only if it is between 1.6-2.0 g/cm3

! subr arguments (out) passed via common block  
!    these are used to duplicate the outputs of yang zhang's original test driver
! they are not really needed in wrf-chem, so just make them local
        real :: ratenuclt        ! j = ternary nucleation rate from napari param. (cm-3 s-1)
        real :: rateloge         ! ln (j)
        real :: cnum_h2so4       ! number of h2so4 molecules in the critical nucleus
        real :: cnum_nh3         ! number of nh3 molecules   in the critical nucleus
        real :: cnum_tot         ! total number of molecules in the critical nucleus
        real :: radius_cluster   ! the radius of cluster (nm)

!       common / ternary_nuc_yzhang_cmn01 /   &
!         ratenuclt, rateloge,   &
!         cnum_h2so4, cnum_nh3, cnum_tot, radius_cluster


! local variables
        integer i
        integer, save :: icase = 0, icase_reldiffmax = 0

        real, parameter :: pi = 3.1415926536
        real, parameter :: avogad = 6.022e23   ! avogadro number (molecules/mole)
        real, parameter :: mw_air = 28.966     ! dry-air mean molecular weight (g/mole)

! dry densities (g/cm3) molecular weights of aerosol 
! ammsulf, ammbisulf, and sulfacid (from mosaic  dens_electrolyte values)
        real, parameter :: dens_ammsulf   = 1.769
        real, parameter :: dens_ammbisulf = 1.78
        real, parameter :: dens_sulfacid  = 1.841

! molecular weights (g/mole) of aerosol ammsulf, ammbisulf, and sulfacid
!    for ammbisulf and sulfacid, use 114 & 96 here rather than 115 & 98
!    because we don't keep track of aerosol hion mass
        real, parameter :: mw_ammsulf   = 132.0
        real, parameter :: mw_ammbisulf = 114.0
        real, parameter :: mw_sulfacid  =  96.0
! molecular weights of aerosol sulfate and ammonium
        real, parameter :: mw_so4a      =  96.0
        real, parameter :: mw_nh4a      =  18.0

        real, save :: reldiffmax = 0.0

        real dens_part                ! "grown" single-particle dry density (g/cm3)
        real duma, dumb, dumc, dume
        real dum_m1, dum_m2, dum_m3, dum_n1, dum_n2, dum_n3
        real fogas, foso4a, fonh4a, fonuma
        real freduce                  ! reduction factor applied to nucleation rate
                                      ! due to limited availability of h2so4 & nh3 gases
        real freducea, freduceb
        real gramaero_per_moleso4a    ! (g dry aerosol)/(mole aerosol so4)
        real mass_part                ! "grown" single-particle mass (g)
        real molenh4a_per_moleso4a    ! (mole aerosol nh4)/(mole aerosol so4)
        real nh3conc_in               ! mixing ratio of nh3 for nucl. calc., pptv
        real so4vol_in                ! concentration of h2so4 for nucl. calc., molecules cm-3
        real qmolnh4a_del_max         ! max production of aerosol nh4 over dtnuc (mole/mole-air)
        real qmolso4a_del_max         ! max production of aerosol so4 over dtnuc (mole/mole-air)
        real vol_cluster              ! critical-cluster volume (cm3)
        real vol_part                 ! "grown" single-particle volume (cm3)


!
! if h2so4 vapor < 4.0e-16 mole/moleair ~= 1.0e4 molecules/cm3, 
! exit with new particle formation = 0
!
        isize_nuc = 1
        qnuma_del = 0.0
        qso4a_del = 0.0
        qnh4a_del = 0.0
        qh2so4_del = 0.0
        qnh3_del = 0.0
        if (qh2so4_avg .le. 4.0e-16) return
        if (qh2so4_cur .le. 4.0e-16) return


!
! make call to napari parameterization routine
!

! convert nh3   from mole/mole-air to ppt
!       qnh3_cur = nh3conc(m)*1.0e-12
        nh3conc_in = qnh3_avg/1.0e-12

! convert h2so4 from mole/mole-air to molecules/cm3
!       qh2so4_cur = (so4vol(k)/avogad) / cair
        so4vol_in  = (qh2so4_avg) * cair * avogad

        call ternary_nuc_napari(   &
            temp_in, rh_in, nh3conc_in, so4vol_in,   &
            ratenuclt, rateloge,   &
            cnum_h2so4, cnum_nh3, cnum_tot, radius_cluster )

! if nucleation rate is less than 0.1 particle/cm3/day,
! exit with new particle formation = 0
        if (ratenuclt .le. 1.0e-6) return


! determine size bin into which the new particles go
! (probably it will always be bin #1, but ...)
        vol_cluster = (pi*4.0/3.0)* (radius_cluster**3) * 1.0e-21
        isize_nuc = 1
        vol_part = volumlo_sect(1)
        if (vol_cluster .le. volumlo_sect(1)) then
           continue
        else if (vol_cluster .ge. volumhi_sect(nsize)) then
           isize_nuc = nsize
           vol_part = volumhi_sect(nsize)
        else
           do i = 1, nsize
              if (vol_cluster .lt. volumhi_sect(i)) then
                 isize_nuc = i
                 vol_part = vol_cluster
                 vol_part = min( vol_part, volumhi_sect(i) )
                 vol_part = max( vol_part, volumlo_sect(i) )
                 exit
              end if
           end do
        end if


!
! determine composition and density of the "grown particles"
! the grown particles are assumed to be liquid
!    (since critical clusters contain water)
!    so any (nh4/so4) molar ratio between 0 and 2 is allowed
! assume that the grown particles will have 
!    (nh4/so4 molar ratio) = min( 2, (nh3/h2so4 gas molar ratio) )
!
        if (qnh3_cur .ge. qh2so4_cur) then
! combination of ammonium sulfate and ammonium bisulfate
! dum_n1 & dum_n2 = mole fractions of the ammsulf & ammbisulf
           dum_n1 = (qnh3_cur/qh2so4_cur) - 1.0
           dum_n1 = max( 0.0, min( 1.0, dum_n1 ) )
           dum_n2 = 1.0 - dum_n1
           dum_n3 = 0.0
        else
! combination of ammonium bisulfate and sulfuric acid
! dum_n2 & dum_n3 = mole fractions of the ammbisulf & sulfacid
           dum_n1 = 0.0
           dum_n2 = (qnh3_cur/qh2so4_cur)
           dum_n2 = max( 0.0, min( 1.0, dum_n2 ) )
           dum_n3 = 1.0 - dum_n2
	end if

        dum_m1 = dum_n1*mw_ammsulf
        dum_m2 = dum_n2*mw_ammbisulf
        dum_m3 = dum_n3*mw_sulfacid
        dens_part = (dum_m1 + dum_m2 + dum_m3)/   &
           ((dum_m1/dens_ammsulf) + (dum_m2/dens_ammbisulf)   &
                                  + (dum_m3/dens_sulfacid))
! 25-jul-2006 - use 'in' value only if it is between 1.6-2.0 g/cm3
	if (abs(dens_nh4so4a-1.8) .le. 0.2) then
	    dens_part = dens_nh4so4a
	else
            dens_nh4so4a = dens_part
	end if
        mass_part  = vol_part*dens_part 
        molenh4a_per_moleso4a = 2.0*dum_n1 + dum_n2  ! (mole aerosol nh4)/(mole aerosol so4)
        gramaero_per_moleso4a = dum_m1 + dum_m2 + dum_m3  ! (g dry aerosol)/(mole aerosol so4)


! max production of aerosol dry mass (g-aero/cm3-air)
        duma = max( 0.0, (ratenuclt*dtnuc*mass_part) )
! max production of aerosol so4 (mole-so4a/cm3-air)
        dumc = duma/gramaero_per_moleso4a
! max production of aerosol so4 (mole-so4a/mole-air)
        dume = dumc/cair
! max production of aerosol so4 (mole/mole-air)
! based on ratenuclt and mass_part
        qmolso4a_del_max = dume

! check if max production exceeds available h2so4 vapor
        freducea = 1.0
        if (qmolso4a_del_max .gt. qh2so4_cur) then
           freducea = qh2so4_cur/qmolso4a_del_max
        end if

! check if max production exceeds available nh3 vapor
        freduceb = 1.0
        if (molenh4a_per_moleso4a .ge. 1.0e-10) then
! max production of aerosol nh4 (ppm) based on ratenuclt and mass_part
           qmolnh4a_del_max = qmolso4a_del_max*molenh4a_per_moleso4a
           if (qmolnh4a_del_max .gt. qnh3_cur) then
              freduceb = qnh3_cur/qmolnh4a_del_max
           end if
        end if
        freduce = min( freducea, freduceb )

! if adjusted nucleation rate is less than 0.1 particle/cm3/day,
! exit with new particle formation = 0
        if (freduce*ratenuclt .le. 1.0e-6) return


! note:  suppose that at this point, freduce = 1.0 (no gas-available 
!    constraints) and molenh4a_per_moleso4a < 2.0
! then it would be possible to condense "additional" nh3 and have
!    (nh3/h2so4 gas molar ratio) < (nh4/so4 aerosol molar ratio) <= 2 
! one could do some additional calculations of 
!    dens_part & molenh4a_per_moleso4a to realize this
! however, the particle "growing" is a crude approximate way to get
!    the new particles to the host code's minimum particle size,
!    so are such refinements worth the effort?


! changes to h2so4 & nh3 gas (in mole/mole-air), limited by amounts available
        duma = 0.9999
        qh2so4_del = min( duma*qh2so4_cur, freduce*qmolso4a_del_max )
        qnh3_del   = min( duma*qnh3_cur, qh2so4_del*molenh4a_per_moleso4a )
        qh2so4_del = -qh2so4_del
        qnh3_del   = -qnh3_del

! changes to so4 & nh4 aerosol (in mole/mole-air)
        qso4a_del = -qh2so4_del
        qnh4a_del =   -qnh3_del
! change to aerosol number (in #/mole-air)
        qnuma_del = (qso4a_del*mw_so4a + qnh4a_del*mw_nh4a)/mass_part


        return
        end subroutine ternary_nuc_mosaic_1box



!-----------------------------------------------------------------------
        subroutine ternary_nuc_napari(   &
           temp_in, rh_in, nh3conc_in, so4vol_in,   &
           ratenuclt, rateloge,   &
           cnum_h2so4, cnum_nh3, cnum_tot, radius_cluster )
!*********************************************************************
!  purpose: calculate ternary nucleation rates based on              *
!           napari et al. (2002) parameterization                    *
!                                                                    *
!  revision history:                                                 *
!     coded by yang zhang,  ncsu, nov. 25, 2004                      *
!                                                                    *
!     14-mar-2006 rce - yang sent this on 10-mar-2005;               *
!         converted to lowercase;                                    *
!         moved the nucleation rate calcs                            *
!             from the main program unit to this subr;               *
!         added temporary variables log_rh, log_nh3conc, log_so4vol; *
!         in the nuc subr, apply limits to the input parameters;     *
!         converted to f90;                                          *
!                                                                    *
!  reference:                                                        *
!     napari, i., m, noppel, h. vehkamaki, . and . m. kulmala,       *
!        parameterization of ternary nucleation rates for            *
!        h2so4-nh3-h2o vapors.  j. geophys. res., 107, 4381, 2002b.  *
!                                                                    *
!  note:                                                             *
!     advantage of this parameterization                             *
!         this parameterization reproduces the ternary nucleation    *
!         rates obtained from the full model within the range of     *
!         one order of magnitude, with a cpu saving by a factor      *
!         of 10e5.                                                   *
!                                                                    *
!     limitations / assumptions of this parameterization             *
!         1. the limits of validity are                              *
!                t = 240-300 k, rh=0.05-0.95                         *
!                [h2so4]=10e4-10e9 cm-3 (4e-4 - 40 ppt at stp)       *
!                [nh3]=0.1-100 ppt,                                  *
!                j=10e-5 - 10e6 cm-3 s-1                             *
!         2. it cannot be used to obtain binary h2o-h2so4 or h2o-nh3 *
!            limit, due to the logarithmic dependencies on rh,       *
!            [h2so4], and [nh3]                                      *
!         3. the fit is the worst at high temperatures ( > 298 k)    *
!            and low nucleation rates ( < 0.01 cm-3 s-1), but it is  *
!            more accurate at significant nucleation rates (0.01-0.1 *
!            cm-3 s-1), thus adequate for simulating ternary         *
!            nucleation in the atmosphere                            *
!                                                                    *
!*********************************************************************

      implicit none

! subr arguments (in)
        real, intent(in) :: temp_in           ! temperature, in k
        real, intent(in) :: rh_in             ! relative humidity, as fraction
        real, intent(in) :: nh3conc_in        ! mixing ratio of nh3, pptv
        real, intent(in) :: so4vol_in         ! concentration of h2so4, cm-3

! subr arguments (out)
        real, intent(out) :: ratenuclt        ! ternary nucleation rate, j, cm-3 s-1
        real, intent(out) :: rateloge         ! ln (j)

        real, intent(out) :: cnum_h2so4       ! number of h2so4 molecules
                                              ! in the critical nucleus
        real, intent(out) :: cnum_nh3         ! number of nh3 molecules
                                              ! in the critical nucleus
        real, intent(out) :: cnum_tot         ! total number of molecules
                                              ! in the critical nucleus
        real, intent(out) :: radius_cluster   ! the radius of cluster, nm

! local variables
        integer ncoeff
        parameter ( ncoeff = 4 )      ! total fitting coefficients at each t
                                      ! corresponds to 3rd order polynomial
        integer npoly
        parameter ( npoly = 20 )      ! total number of polynomial functions

        integer n                     ! loop index for functions of polynomial

! polynomial functions
        real f ( npoly )

! coefficients of polynomials
        real a  ( ncoeff, npoly )

        real temp, rh, nh3conc, so4vol     ! bounded values of the input args
        real log_rh, log_nh3conc, log_so4vol

! coefficients of polynomials fi(t), i = 1,  20
        data a  / -0.355297,  -33.8449,      0.34536,    -8.24007e-4,   &
                   3.13735,    -0.772861,    5.61204e-3, -9.74576e-6,   &
                  19.0359,     -0.170957,    4.79808e-4, -4.14699e-7,   &
                   1.07605,     1.48932,    -7.96052e-3,  7.61229e-6,   &
                   6.0916,     -1.25378,     9.39836e-3, -1.74927e-5,   &
                   0.31176,     1.64009,    -3.43852e-3, -1.09753e-5,   &
                  -2.00738e-2, -0.752115,    5.25813e-3, -8.98038e-6,   &
                   0.165536,    3.26623,    -4.89703e-2,  1.46967e-4,   &
                   6.52645,    -0.258002,    1.43456e-3, -2.02036e-6,   &
                   3.68024,    -0.204098,    1.06259e-3, -1.2656e-6 ,   &
                  -6.6514e-2,  -7.82382,     1.22938e-2,  6.18554e-5,   &
                   0.65874,     0.190542,   -1.65718e-3,  3.41744e-6,   &
                   5.99321e-2,  5.96475,    -3.62432e-2,  4.93337e-5,   &
                  -0.732731,   -1.84179e-2,  1.47186e-4, -2.37711e-7,   &
                   0.728429,    3.64736,    -2.7422e-2,   4.93478e-5,   &
                  41.3016,     -0.35752,     9.04383e-4, -5.73788e-7,   &
                  -0.160336,    8.89881e-3, -5.39514e-5,  8.39522e-8,   &
                   8.57868,    -0.112358,    4.72626e-4, -6.48365e-7,   &
                   5.30167e-2, -1.98815,     1.57827e-2, -2.93564e-5,   &
                  -2.32736,     2.34646e-2, -7.6519e-5,   8.0459e-8 /


!
! apply limits to input parameters
!
        temp    = max( 240.15, min (300.15, temp_in ) )
        rh      = max( 0.05,   min (0.95,   rh_in ) )
        so4vol  = max( 1.0e4,  min (1.0e9,  so4vol_in ) )
        nh3conc = max( 0.1,    min (100.0,  nh3conc_in ) )

!
! calculate the functions of polynomials, fi(t), i = 1,  npoly
! at temperature j
!
        do n = 1, npoly
            f ( n )   = a ( 1, n ) + a ( 2, n ) * temp   &
                      + a ( 3, n ) * ( temp ) ** 2.0   &
                      + a ( 4, n ) * ( temp ) ** 3.0
        end do

!
! calculate ln (j)
!
        log_rh = log ( rh )
        log_nh3conc = log ( nh3conc )
        log_so4vol = log ( so4vol )
        rateloge = -84.7551   &
                 + f ( 1 ) / log_so4vol   &
                 + f ( 2 )  * ( log_so4vol )   &
                 + f ( 3 )  * ( log_so4vol ) **2.0   &
                 + f ( 4 )  * ( log_nh3conc )   &
                 + f ( 5 )  * ( log_nh3conc ) **2.0   &
                 + f ( 6 )  * rh   &
                 + f ( 7 )  * ( log_rh )   &
                 + f ( 8 )  * ( log_nh3conc /   &
                              log_so4vol )   &
                 + f ( 9 )  * ( log_nh3conc *   &
                              log_so4vol )   &
                 + f ( 10 ) * rh  *   &
                              ( log_so4vol )   &
                 + f ( 11 ) * ( rh /   &
                              log_so4vol )   &
                 + f ( 12 ) * ( rh *   &
                              log_nh3conc )   &
                 + f ( 13 ) * ( log_rh /   &
                              log_so4vol )   &
                 + f ( 14 ) * ( log_rh *   &
                              log_nh3conc )   &
                 + f ( 15 ) * (( log_nh3conc ) ** 2.0   &
                              / log_so4vol )   &
                 + f ( 16 ) * ( log_so4vol *   &
                              ( log_nh3conc ) ** 2.0 )   &
                 + f ( 17 ) * (( log_so4vol ) ** 2.0 *   &
                              log_nh3conc )   &
                 + f ( 18 ) * ( rh  *   &
                              ( log_nh3conc ) ** 2.0 )   &
                 + f ( 19 ) * ( rh  *  log_nh3conc   &
                              / log_so4vol )   &
                 + f ( 20 ) * (( log_so4vol ) ** 2.0 *   &
                              ( log_nh3conc ) ** 2.0 )

        ratenuclt = exp ( rateloge )
!
! calculate number of molecules and critical radius of the cluster
!
        cnum_h2so4 = 38.1645 + 0.774106 * rateloge   &
                   + 0.00298879 * ( rateloge ) ** 2.0   &
                   - 0.357605 * temp   &
                   - 0.00366358 * temp * rateloge   &
                   + 0.0008553 * ( temp ) ** 2.0

        cnum_nh3   = 26.8982 + 0.682905 * rateloge   &
                   + 0.00357521 * ( rateloge ) ** 2.0   &
                   - 0.265748 * temp   &
                   - 0.00341895 * temp * rateloge   &
                   + 0.000673454 * ( temp ) ** 2.0

        cnum_tot   = 79.3484 + 1.7384 * rateloge   &
                   + 0.00711403 * ( rateloge ) ** 2.0   &
                   - 0.744993 * temp   &
                   - 0.00820608 * temp * rateloge   &
                   + 0.0017855 * ( temp ) ** 2.0

        radius_cluster = 0.141027 - 0.00122625 * rateloge   &
                   - 7.82211e-6 * ( rateloge ) ** 2.0   &
                   - 0.00156727 * temp   &
                   - 0.00003076 * temp * rateloge   &
                   + 0.0000108375 * ( temp ) ** 2.0

        return
        end subroutine ternary_nuc_napari



!-----------------------------------------------------------------------
        subroutine wexler_nuc_mosaic_1box(   &
           dtnuc, temp_in, rh_in, cair,   &
           qh2so4_avg, qh2so4_cur, qnh3_avg, qnh3_cur,   &
           nsize, maxd_asize, volumlo_sect, volumhi_sect,   &
           isize_nuc, qnuma_del, qso4a_del, qnh4a_del,   &
           qh2so4_del, qnh3_del, dens_nh4so4a )
!.......................................................................
!
! calculates new particle production from h2so4-h2o binary nucleation
!    over timestep dtnuc, using the wexler et al. (1994) parameterization
!
! the size of new particles is the lower-bound size of the host code's 
!    smallest size bin.  their composition is so4 and nh4, since the nuclei
!    would incorporate nh3 as they grow from ~1 nm to the lower-bound size.
!    (the new particle composition can be forced to pure so4 by setting
!    the qnh3_avg & qnh3_cur input arguments to 0.0).
!
! revision history
!    coded by rc easter, pnnl, 20-mar-2006
!
! key routines called: none
!
! references:
!    wexler, a. s., f. w. lurmann, and j. h. seinfeld,
!       modelling urban and regional aerosols -- i.  model development,
!       atmos. environ., 28, 531-546, 1994.
!
!.......................................................................
        implicit none

! subr arguments (in)
        real, intent(in) :: dtnuc             ! nucleation time step (s)
        real, intent(in) :: temp_in           ! temperature, in k
        real, intent(in) :: rh_in             ! relative humidity, as fraction
        real, intent(in) :: cair              ! dry-air molar density (mole-air/cm3)

        real, intent(in) :: qh2so4_avg, qh2so4_cur   ! gas h2so4 mixing ratios (ppm)
        real, intent(in) :: qnh3_avg, qnh3_cur       ! gas nh3 mixing ratios (ppm)
             ! qxxx_cur = current value (at end of condensation)
             ! qxxx_avg = average value (from start to end of condensation)

        integer, intent(in) :: nsize                    ! number of aerosol size bins
        integer, intent(in) :: maxd_asize               ! dimension for volumlo_sect, ...
        real, intent(in) :: volumlo_sect(maxd_asize)    ! dry volume at lower bnd of bin (cm3)
        real, intent(in) :: volumhi_sect(maxd_asize)    ! dry volume at upper bnd of bin (cm3)

! subr arguments (out)
        integer, intent(out) :: isize_nuc     ! size bin into which new particles go
        real, intent(out) :: qnuma_del        ! change to aerosol number mixing ratio (#/kg)
        real, intent(out) :: qso4a_del        ! change to aerosol so4 mixing ratio (ug/kg)
        real, intent(out) :: qnh4a_del        ! change to aerosol nh4 mixing ratio (ug/kg)
        real, intent(out) :: qh2so4_del       ! change to gas h2so4 mixing ratio (ppm)
        real, intent(out) :: qnh3_del         ! change to gas nh3 mixing ratio (ppm)
                                              ! aerosol changes are > 0; gas changes are < 0

! subr arguments (inout)
        real, intent(inout) :: dens_nh4so4a   ! dry-density of the new nh4-so4 aerosol mass (g/cm3)
                                              ! use 'in' value only if it is between 1.6-2.0 g/cm3

! local variables
        integer i
        integer, save :: icase = 0, icase_reldiffmax = 0

        real, parameter :: pi = 3.1415926536
        real, parameter :: avogad = 6.022e23   ! avogadro number (molecules/mole)
        real, parameter :: mw_air = 28.966     ! dry-air mean molecular weight (g/mole)

! dry densities (g/cm3) molecular weights of aerosol 
! ammsulf, ammbisulf, and sulfacid (from mosaic  dens_electrolyte values)
        real, parameter :: dens_ammsulf   = 1.769
        real, parameter :: dens_ammbisulf = 1.78
        real, parameter :: dens_sulfacid  = 1.841

! molecular weights (g/mole) of aerosol ammsulf, ammbisulf, and sulfacid
!    for ammbisulf and sulfacid, use 114 & 96 here rather than 115 & 98
!    because we don't keep track of aerosol hion mass
        real, parameter :: mw_ammsulf   = 132.0
        real, parameter :: mw_ammbisulf = 114.0
        real, parameter :: mw_sulfacid  =  96.0
! molecular weights of aerosol sulfate and ammonium
        real, parameter :: mw_so4a      =  96.0
        real, parameter :: mw_nh4a      =  18.0

        real, save :: reldiffmax = 0.0

        real ch2so4_crit              ! critical h2so4 conc (ug/m3)
        real dens_part                ! "grown" single-particle dry density (g/cm3)
        real duma, dumb, dumc, dume
        real dum_m1, dum_m2, dum_m3, dum_n1, dum_n2, dum_n3
        real fogas, foso4a, fonh4a, fonuma
        real mass_part                ! "grown" single-particle mass (g)
        real molenh4a_per_moleso4a    ! (mole aerosol nh4)/(mole aerosol so4)
        real qh2so4_crit              ! critical h2so4 mixrat (ppm)
        real qh2so4_avail             ! amount of h2so4 available for new particles (ppm)
        real vol_part                 ! "grown" single-particle volume (cm3)


!
! initialization output arguments with "zero nucleation" values
!
        isize_nuc = 1
        qnuma_del = 0.0
        qso4a_del = 0.0
        qnh4a_del = 0.0
        qh2so4_del = 0.0
        qnh3_del = 0.0

!
! calculate critical h2so4 concentration (ug/m3) and mixing ratio (mole/mole-air)
!
        ch2so4_crit = 0.16 * exp( 0.1*temp_in - 3.5*rh_in - 27.7 )
! ch2so4 = (ug-h2so4/m3-air)
! ch2so4*1.0e-12/mwh2so4 = (mole-h2so4/cm3-air)
        qh2so4_crit = (ch2so4_crit*1.0e-12/98.0)/cair
        qh2so4_avail = (qh2so4_cur - qh2so4_crit)*dtnuc

! if "available" h2so4 vapor < 4.0e-18 mole/mole-air ~= 1.0e2 molecules/cm3, 
! exit with new particle formation = 0
        if (qh2so4_avail .le. 4.0e-18) then
           return
        end if

! determine size bin into which the new particles go
        isize_nuc = 1
        vol_part = volumlo_sect(1)

!
! determine composition and density of the "grown particles"
! the grown particles are assumed to be liquid 
!    (since critical clusters contain water)
!    so any (nh4/so4) molar ratio between 0 and 2 is allowed
! assume that the grown particles will have 
!    (nh4/so4 molar ratio) = min( 2, (nh3/h2so4 gas molar ratio) )
!
        if (qnh3_cur .ge. qh2so4_avail) then
! combination of ammonium sulfate and ammonium bisulfate
! dum_n1 & dum_n2 = mole fractions of the ammsulf & ammbisulf
           dum_n1 = (qnh3_cur/qh2so4_avail) - 1.0
           dum_n1 = max( 0.0, min( 1.0, dum_n1 ) )
           dum_n2 = 1.0 - dum_n1
           dum_n3 = 0.0
        else
! combination of ammonium bisulfate and sulfuric acid
! dum_n2 & dum_n3 = mole fractions of the ammbisulf & sulfacid
           dum_n1 = 0.0
           dum_n2 = (qnh3_cur/qh2so4_avail)
           dum_n2 = max( 0.0, min( 1.0, dum_n2 ) )
           dum_n3 = 1.0 - dum_n2
	end if

        dum_m1 = dum_n1*mw_ammsulf
        dum_m2 = dum_n2*mw_ammbisulf
        dum_m3 = dum_n3*mw_sulfacid
        dens_part = (dum_m1 + dum_m2 + dum_m3)/   &
           ((dum_m1/dens_ammsulf) + (dum_m2/dens_ammbisulf)   &
                                  + (dum_m3/dens_sulfacid))
! 25-jul-2006 - use 'in' value only if it is between 1.6-2.0 g/cm3
	if (abs(dens_nh4so4a-1.8) .le. 0.2) then
	    dens_part = dens_nh4so4a
	else
            dens_nh4so4a = dens_part
	end if
        mass_part  = vol_part*dens_part 
        molenh4a_per_moleso4a = 2.0*dum_n1 + dum_n2


! changes to h2so4 & nh3 gas (in mole/mole-air), limited by amounts available
        duma = 0.9999
        qh2so4_del = min( duma*qh2so4_cur, qh2so4_avail )
        qnh3_del   = min( duma*qnh3_cur, qh2so4_del*molenh4a_per_moleso4a )
        qh2so4_del = -qh2so4_del
        qnh3_del   = -qnh3_del

! changes to so4 & nh4 aerosol (in mole/mole-air)
        qso4a_del = -qh2so4_del
        qnh4a_del =   -qnh3_del
! change to aerosol number (in #/mole-air)
        qnuma_del = (qso4a_del*mw_so4a + qnh4a_del*mw_nh4a)/mass_part


        return
        end subroutine wexler_nuc_mosaic_1box


!-----------------------------------------------------------------------
        subroutine riccobono_nuc_msa(   &
           dtnuc, temp_in, rh_in, cair, qh2so4_avg, qh2so4_cur, qnh3_avg, &
           maxd_asize, volumlo_sect, volumhi_sect,   &
           qnuma_del, qh2so4_del, qmsa_avg, qmsa_cur, qmsa_del, qmsaa_del, j2_msa, nuc_msa_fac )

           ! PURPOSE: Compute new particle formation rates for the nucleation of MSA with sulfuric acid, following Riccobono et al. 2014 (2)
           ! the size of new particles is the lower-bound size of the smallest size bin, i.e. all newly nucleated particles coagulate into the first bin
           
           ! REFERENCES
           ! (1) Riccobono et al. (2014) - doi: 10.1126/science.1243527
           ! (2) Hodshire et al. (2019) - doi: 10.5194/acp-19-3137-2019

           implicit none

           real, intent(in) :: dtnuc ! nucleation time step (s)
           real, intent(in) :: temp_in ! temperature, in k
           real, intent(in) :: rh_in ! relative humidity, as fraction
           real, intent(in) :: cair ! dry-air molar density (mole-air/cm3)

           real, intent(in) :: qh2so4_avg, qh2so4_cur, qnh3_avg
           real, intent(in) :: qmsa_avg, qmsa_cur

           integer, intent(in) :: maxd_asize ! dimension for volumlo_sect
           real, intent(in) :: volumlo_sect(maxd_asize) ! dry volume at lower bnd of bin (cm3)
           real, intent(in) :: volumhi_sect(maxd_asize) ! dry volume at upper bnd of bin (m3)

           real, intent(inout) :: qnuma_del ! change to aerosol number mixing ratio (#/mole-air)
           real, intent(out) :: qh2so4_del ! change to gas h2so4 mixing ratio (mole/mole-air)
           real, intent(out) :: qmsa_del ! change to gas msa mixing ratio (mole/mole-air)
           real, intent(out) :: qmsaa_del ! change to msa aerosol mixing ratio (mole/mole-air) 
           real, intent(in) :: nuc_msa_fac ! correction factor to nucleation rate - unitless

           real, parameter :: dens_msaa = 1.48 ! g/cm3

           real, parameter :: pi = 3.1415926536
           real, parameter :: avogad = 6.022e23 ! avogadro number (molecules/mole)
           real, parameter :: mw_air = 28.966 ! dry-air mean molecular weight (g/mole)
           real, parameter :: mw_msaa = 96.0 ! msa aerosol molecular weight (g/mole)

           real, parameter :: k_sa1 = 1.1e-14 ! cm3/s
           real, parameter :: k_sa2 = 3.2e-14 ! cm3/s
           real, parameter :: k_m = 3.27e-21 ! cm6/s

           real, parameter :: qnh3_crit_h2019 = 10.0 ! threshold nh3 concentration from (2)
           ! parameters for critical temperature formula from (2) - K
           real, parameter :: a_h = 252.0
           real, parameter :: b_h = 0.619
           real, parameter :: c_h = 0.0349
           real, parameter :: d_h = 0.00056
           real, parameter :: e_h = 0.00000332
	   
           real vol_part  ! "grown" single-particle volume (cm3)

           real, intent(out) :: j2_msa ! nucleation rate - #/cm3
           real ppm_to_num ! conversion factor mole/mole_air to #/cm3
           real qh2so4_avg_num ! number mixing ratio of h2so4 - #/cm3
           real qmsa_avg_num ! number mixing ratio of msa gas - #/cm3
           real t_trans ! transition temperature - K
           real rh100 ! relative humidity as percentage
           real qnh3_avg_ppt ! mixing ratio of nh3 - ppt

           qmsa_del = 0.0
           qmsaa_del = 0.0
           rh100 = rh_in*100.0
           qnh3_avg_ppt = qnh3_avg*1e9
	   
           j2_msa = 0.0

           ! input gas mixing ratios are in ppm
           ppm_to_num = cair*avogad ! conversion from mole/mole-air to #/cm3
           qh2so4_avg_num = qh2so4_avg*ppm_to_num
           qmsa_avg_num = qmsa_avg*ppm_to_num

           ! nucleation rate from Riccobono et al. (2014) - ref (1)
           j2_msa = k_m*qh2so4_avg_num*qh2so4_avg_num*qmsa_avg_num
           ! apply tuning factor
           j2_msa = j2_msa*nuc_msa_fac

           ! convert nucleation rate to wrf-chem expected units
           j2_msa = max(0.0, j2_msa*dtnuc) ! # particles per cm3 per nucleation time step
           j2_msa = j2_msa/cair ! # per mole-air per nucleation time step


           ! temperature threshold from Hodshire et al. (2019) - ref (2)
           t_trans = a_h-b_h*rh100+c_h*rh100*rh100-d_h*rh100*rh100*rh100+e_h*rh100*rh100*rh100*rh100


           ! New particle formation is only possible in specific temperature, RH conditions,
           ! and depends on the availability of bases, represented here by NH3 (Hodshire et al. (2019) - ref (2))

           ! In low base conditions, MSA acts as a Volatile or semivolatile organic compound, and does not contribute to new particle formation
           if (qnh3_avg_ppt .le. qnh3_crit_h2019) then
              j2_msa = 0.0
           else
           ! In excess-base conditions, MSA volatility depends on T and RH
           ! If t>t_trans, then MSA is treated as a semivolatile species that undergoes condensation but not nucleation 
              if (temp_in .ge. t_trans) then
                 j2_msa = 0.0
              endif
           endif

           vol_part = volumlo_sect(1)

           qmsaa_del = j2_msa*vol_part*dens_msaa/mw_msaa/(mw_air) ! # particles per kg to mole msa per mole-air
           qmsa_del = -qmsaa_del
           qnuma_del = qnuma_del+j2_msa

           j2_msa = j2_msa/(mw_air)*1.0e3 ! # per kg air per nucleation time step for output

        return
        end subroutine riccobono_nuc_msa
!.......................................................................


!-----------------------------------------------------------------------



	
!-----------------------------------------------------------------------
        subroutine kulmala1998(&
           dtnuc, temp_in, rh_in, cair, qv_in,  &
           qh2so4_avg, qh2so4_cur, &
           isize_nuc, maxd_asize, volumlo_sect, volumhi_sect,   &
           qnuma_del, qh2so4_del, qso4a_del, dens_so4_aer )
!-----------------------------------------------------------------------
! from Kulmala et al., 1998 - Parameterizations for sulfuric acid/water nucleation rates

	implicit none

	real, intent(in) :: dtnuc ! nucleation time step (s)
    real, intent(in) :: temp_in ! temperature, in k
    real, intent(in) :: rh_in ! relative humidity, as fraction
    real, intent(in) :: cair ! dry-air molar density (mole-air/cm3)
	real, intent(in) :: qv_in ! specific humidity
	real, intent(in) :: dens_so4_aer
	
    real, intent(in) :: qh2so4_avg, qh2so4_cur   ! gas h2so4 mixing ratios
    real, intent(out) :: qh2so4_del, qso4a_del ! change to gas h2so4 and so4 aerosol
	real, intent(inout) :: qnuma_del ! change to aerosol number in #/cm3

	integer, intent(in) :: maxd_asize ! dimension for volumlo_sect
    real, intent(in) :: volumlo_sect(maxd_asize) ! dry volume at lower bnd of bin (cm3)        
	real, intent(in) :: volumhi_sect(maxd_asize) ! dry volume at upper bnd of bin (m3)         
	integer, intent(out) :: isize_nuc ! size bin into which new particles go

	real vol_part, ppm_to_num
	
	real xal, flux
	real zp, ksi, jsulf, ra, jloc
	real nac, delta, psat, nsulf
	real te, rh, nav, nwv
	real, parameter :: avogad = 6.022e23 ! avogadro number (molecules/mole)
	real, parameter :: R = 8.31446261815324   ! molar gas constant R in [J mol-1 K-1]
	real, parameter :: t0k = 273.15
	real, parameter :: mh2o=18.0
	real, parameter :: mh2so4=98.0
	real, parameter :: mw_air=28.966 ! dry-air mean molecular weight (g/mole)
	real, parameter :: mw_so4a=96.0 
	
	te = temp_in
	rh = rh_in

	jloc = 0.0
	qh2so4_del = 0.0
	qso4a_del = 0.0
	
	ppm_to_num = cair*avogad ! conversion factor from mole/mole-air to #/cm3
	nav = qh2so4_avg*ppm_to_num
	nwv = qv_in*ppm_to_num 

	te = max(min(te,298.0),233.0)
	rh = max(min(rh,1.0),0.1)

	! relative acidity calculation
	psat=exp(27.7849-10156.0/te)
	zp=9.869*nav*R*te/avogad
	ra=zp/psat
 
	xal = 1.2233-(0.0154*ra)/(ra+rh)+0.0102*log(nav)-0.0415*log(nwv)+0.0016*te

	nac = exp(-14.5125+0.1335*te-10.5462*rh+1958.4*rh/te)

	nsulf = log(nav/nac)
	delta = 1.0+(te-t0k)/t0k
	ksi = 25.1289*nsulf-4890.8*nsulf/te-1743.3/te-2.2479*delta*nsulf*rh+7643.4*xal/te-1.9712*xal*delta/rh

	jsulf = exp(ksi) ! #/cm3/s
	jsulf = max(min(jsulf,1.0e5),1.0e-5)

	jloc = jsulf*dtnuc/cair ! #/mole-air per nucleation time step
	
	isize_nuc = 1
	vol_part = volumlo_sect(1)
	
	qso4a_del = jloc*vol_part*dens_so4_aer/mw_so4a/(mw_air) ! #/mole-air  to mole_so4/mole-air
	qh2so4_del = -qso4a_del
	qnuma_del = qnuma_del+jloc

		
	return
	end subroutine kulmala1998

!----------------------------------------



!-----------------------------------------------------------------------
subroutine vehkamaki2002(&
           dtnuc, temp_in, rh_in, cair,  &
           qh2so4_avg, qh2so4_cur, &
           isize_nuc, maxd_asize, volumlo_sect, volumhi_sect,   &
           qnuma_del, qh2so4_del, qso4a_del, dens_so4_aer )
!----------------------------------------------------------------------- 
! from Vehkamaki et al., 2002 - An improved parameterization for sulfuric acid–water
! nucleation rates for tropospheric and stratospheric conditions, 10.1029/2002JD002184
	
    implicit none

    real, intent(in) :: dtnuc ! nucleation time step (s)
    real, intent(in) :: temp_in ! temperature, in k
    real, intent(in) :: rh_in ! relative humidity, as fraction
    real, intent(in) :: cair ! dry-air molar density (mole-air/cm3)
    real, intent(in) :: dens_so4_aer

    real, intent(in) :: qh2so4_avg, qh2so4_cur ! gas h2so4 mixing ratios
    real, intent(out) :: qh2so4_del, qso4a_del ! change to gas h2so4 and so4 aerosol
    real, intent(inout) :: qnuma_del ! change to aerosol number in #/cm3
	
    integer, intent(in) :: maxd_asize ! dimension for volumlo_sect
    real, intent(in) :: volumlo_sect(maxd_asize) ! dry volume at lower bnd of bin (cm3)
    real, intent(in) :: volumhi_sect(maxd_asize) ! dry volume at upper bnd of bin (m3)
    integer, intent(out) :: isize_nuc ! size bin into which new particles go

    real vol_part, ppm_to_num, jloc

	real :: acoe, bcoe, ccoe, dcoe, ecoe, fcoe, gcoe, hcoe, icoe, jcoe
	real :: te,rh 
	real :: so4vol ! concentration of h2so4 (molecules cm-3
	real :: crit_x
	real :: tmpa
	real :: ratenucl ! binary nucleation rate, j (# cm-3 s-1)
	real :: cnum_tot ! total number of molecules
	real :: fj

	real, parameter :: avogad = 6.022e23

    real, parameter :: mw_air=28.966 ! dry-air mean molecular weight (g/mole)
    real, parameter :: mw_so4a=96.0
	
!       initial conditions
	te=temp_in
	rh=rh_in

	ppm_to_num = cair*avogad ! conversion factor from mole/mole-air to #/cm3
	so4vol=qh2so4_avg*ppm_to_num

	te = max(min(te,300.15),230.15)
    rh = max(min(rh,1.0),0.1)
    so4vol= max(min(so4vol,1.0e11),1.0e4)
		
!       calc sulfuric acid mole fraction in critical cluster
    crit_x = 0.740997d0 - 0.00266379d0 * te              &
         	- 0.00349998d0 * log (so4vol)               &
         	+ 0.0000504022d0 * te * log (so4vol)        &
         	+ 0.00201048d0  * log (rh)                  &
         	- 0.000183289d0 * te * log (rh)             &
         	+ 0.00157407d0  * (log (rh)) ** 2.0d0       &
         	- 0.0000179059d0 * te * (log (rh)) ** 2.0d0 &
         	+ 0.000184403d0  * (log (rh)) ** 3.0d0      &
         	- 1.50345d-6 * te * (log (rh)) ** 3.0d0

!       calc nucleation rate
	acoe    = 0.14309d0+2.21956d0*te    &
	 - 0.0273911d0 * te**2.0d0 &
	 + 0.0000722811d0 * te**3.0d0 + 5.91822d0/crit_x
	 
	bcoe    = 0.117489d0 + 0.462532d0 *te &
	 - 0.0118059d0 * te**2.0d0   &
	 + 0.0000404196d0 * te**3.0d0 + 15.7963d0/crit_x
	 
    ccoe    = -0.215554d0-0.0810269d0 * te &
	 + 0.00143581d0 * te**2.0d0   &
	 - 4.7758d-6 * te**3.0d0      &
	 - 2.91297d0/crit_x
	 
    dcoe    = -3.58856d0+0.049508d0 * te &
	 - 0.00021382d0 * te**2.0d0 &
	 + 3.10801d-7 * te**3.0d0   &
	 - 0.0293333d0/crit_x
	 
	ecoe    = 1.14598d0 - 0.600796d0 * te  &
	 + 0.00864245d0 * te**2.0d0   &
	 - 0.0000228947d0 * te**3.0d0 &
	 - 8.44985d0/crit_x
	 
	fcoe    = 2.15855d0 + 0.0808121d0 * te &
	 -0.000407382d0 * te**2.0d0   &
	 -4.01957d-7 * te**3.0d0      &
	 + 0.721326d0/crit_x
	 
	gcoe    = 1.6241d0 - 0.0160106d0 * te  &
	 + 0.0000377124d0 * te**2.0d0 &
	 + 3.21794d-8 * te**3.0d0     &
	 - 0.0113255d0/crit_x
	 
	hcoe    = 9.71682d0 - 0.115048d0 * te &
	 + 0.000157098d0 * te**2.0d0 &
	 + 4.00914d-7 * te**3.0d0    &
	 + 0.71186d0/crit_x
	 
	icoe    = -1.05611d0 + 0.00903378d0 * te &
	 - 0.0000198417 * te**2.0d0     &
	 + 2.46048d-8  * te**3.0d0      &
	 - 0.0579087d0/crit_x
	 
	jcoe    = -0.148712d0 + 0.00283508d0 * te &
	 - 9.24619d-6  * te**2.0d0       &
	 + 5.00427d-9 * te**3.0d0        &
	 - 0.0127081d0/crit_x

	tmpa    = (                                    &
	 acoe                                 &
	 + bcoe * log (rh)                    &
	 + ccoe * ( log (rh))**2.0d0          &
	 + dcoe * ( log (rh))**3.0d0          &
	 + ecoe * log (so4vol)                &
	 + fcoe * (log (rh)) * (log (so4vol)) &
	 + gcoe * ((log (rh) ) **2.0d0)       &
	 * (log (so4vol))                     &
	 + hcoe * (log (so4vol)) **2.0d0      &
	 + icoe * log (rh)                    &
	 * ((log (so4vol)) **2.0d0)           &
	 + jcoe * (log (so4vol)) **3.0d0      &
	 )

	ratenucl = exp ( tmpa )

!       calc number of molecules in critical cluster
	acoe    = -0.00295413d0 - 0.0976834d0 * te &
	 + 0.00102485d0 * te**2.0d0       &
	 - 2.18646d0-6 * te**3.0d0 - 0.101717d0 / crit_x
	 
	bcoe    = -0.00205064d0 - 0.00758504d0 * te &
	 + 0.000192654d0 * te**2.0         &
	 - 6.7043d-7 * te**3.0d0 - 0.255774d0 / crit_x
	 
    ccoe    = +0.00322308d0 + 0.000852637d0 * te &
	 - 0.0000154757d0 * te**2.0d0       &
	 + 5.66661d-8 * te**3.0d0           &
	 + 0.0338444d0 / crit_x
	 
    dcoe    = +0.0474323d0 - 0.000625104d0 * te &
	 + 2.65066d-6 * te**2.0d0          &
	 - 3.67471d-9 * te**3.0d0          &
	 - 0.000267251d0 / crit_x

	ecoe    = -0.0125211d0 + 0.00580655d0 * te &
	 - 0.000101674 * te**2.0d0        &
	 + 2.88195d-7 * te**3.0d0         &
	 + 0.0942243d0/crit_x
	 
	fcoe    = -0.038546d0 - 0.000672316 * te &
	 + 2.60288d-6 * te**2.0d0       &
	 + 1.19416d-8 * te**3.0d0       &
	 - 0.00851515d0 / crit_x
	 
	gcoe    = -0.0183749d0 + 0.000172072d0 * te &
	 - 3.71766d-7 * te**2.0d0          &
	 - 5.14875d-10 * te**3.0d0         &
	 + 0.00026866d0 / crit_x
	 
	hcoe    = -0.0619974d0 + 0.000906958d0 * te &
	 - 9.11728d-7 * te**2.0d0          &
	 - 5.36796d-9 * te**3.0d0          &
	 - 0.00774234d0 / crit_x
	 
	icoe    = +0.0121827d0 - 0.00010665d0 * te &
	 + 2.5346d-7 * te**2.0d0          &
	 - 3.63519d-10 * te**3.0d0        &
	 + 0.000610065d0 / crit_x
	 
	jcoe    = +0.000320184d0 - 0.0000174762d0 * te &
	 + 6.06504d-8 * te**2.0d0             &
	 - 1.4177d-11 * te**3.0d0             &
	 + 0.000135751d0 / crit_x
	 
	cnum_tot = exp (   &
	 acoe     &
	 + bcoe * log (rh)                    &
	 + ccoe * ( log (rh))**2.0d0          &
	 + dcoe * ( log (rh))**3.0d0          &
	 + ecoe * log (so4vol)                &
	 + fcoe * (log (rh)) * (log (so4vol)) &
	 + gcoe * ((log (rh) ) **2.0d0)       &
	 * (log (so4vol))                     &
	 + hcoe * (log (so4vol)) **2.0d0      &
	 + icoe * log (rh)                    &
	 * ((log (so4vol)) **2.0d0)    &
	 + jcoe * (log (so4vol)) **3.0d0      &
	 )

	if (ratenucl < 1.0d-7 .or. cnum_tot < 4.d0) return
	if (ratenucl > 1.0d11) ratenucl=1.0d11

	fj=ratenucl

	jloc = fj*dtnuc/cair

	isize_nuc = 1
	vol_part = volumlo_sect(1)

	qso4a_del = jloc*vol_part*dens_so4_aer/mw_so4a/(mw_air) ! #/mole-air  to mole_so4/mole-air
	qh2so4_del = -qso4a_del
	qnuma_del = qnuma_del+jloc

end subroutine vehkamaki2002

!-------------------------
	

	
	end module module_mosaic_newnuc
