!=======================================================================
! module_mp_microwcp
!=======================================================================
! HOW TO READ THE VARIABLE NAMES
!-----------------------------------------------------------------------
! Prognostic fields:
!   qv,qc,qr,qi,qs,qg   mixing ratios   [kg kg-1]  vapour / cloud /
!                       rain / cloud-ice / snow / graupel
!   nc,nr,ni            number concentrations [kg-1] for the 2-moment
!                       categories (cloud,rain,ice). snow & graupel are
!                       1-moment: their number is diagnosed, not carried.
!   T                   temperature [K]   (th*pii inside the column code)
!
! Environment / grid (read-only inside the physics):
!   P    pressure [Pa]      rho  air density [kg m-3]
!   w    vertical velocity [m s-1]   dz  layer thickness [m]
!   pii  Exner function (p/P00)**(R/cp)   th  potential temperature [K]
!
! The "q"/"n" prefixes:  q = mass mixing ratio,  n = number.
! The category letters:  c=cloud, r=rain, i=cloud-ice, s=snow, g=graupel.
!   (so qg = graupel mass, ni = ice number, etc.)
!
! Category INDEX constants (for the per-category parameter arrays):
!   IC=cloud  IR=rain  II=ice  IS=snow  IG=graupel  ; NCAT=5
!   e.g. av(IS) is the snow fall-speed coefficient.
!
! Per-category PSD (particle size distribution) parameters / diagnostics:
!   mu,bm,am   shape, mass-diameter exponent/coeff:  mass = am*D^bm
!   av,bv      fall-speed law:  V(D) = av*D^bv
!   Gb1,G1,Gq,GN   precomputed Gamma-function ratios (see mp_microwcp_init)
!   lam(:,c)   PSD slope lambda [m-1]    (n(D) ~ D^mu * exp(-lam*D))
!   Dm (:,c)   mass-weighted mean diameter [m]
!   Nb (:,c)   bulk number concentration [m-3]  (= specific n * rho)
!     "bulk" = per m3 of air;  "specific" = per kg of air (n*rho = bulk).
!
! Local tendency bookkeeping inside each process subroutine:
!   <x>0   a FROZEN copy of field x at entry (e.g. qc0, T0k). Every
!          sub-process rate is computed from these frozen values, so the
!          order of sub-processes does not matter (matches the prototype's
!          "sum all tendencies, then apply once" design; Thompson does the
!          same with its prr_*/prv_* accumulators).
!   a<x>   the ACCUMULATED tendency for field x [unit s-1] (e.g. aqc, aT,
!          anr). After all sub-processes: x = max(x0 + a<x>*dt, 0).
!   d<x>   a single-process rate contribution before it is added into a<x>.
!
! Suffix "1d":  a single vertical column (length nz) carved out of a 3-D
!   tile array by the driver (matches Thompson's qv1d, t1d, ... naming).
!
! Thresholds:  Q*_MIN / N_MIN are floors below which a category is treated
!   as absent (skip the process) or clamped, to avoid divide-by-zeros.
!=======================================================================
module module_mp_microwcp
  implicit none
  private

  public :: mp_microwcp_init, mp_microwcp_aero_init, mp_microwcp_column, mp_microwcp_driver

  integer, parameter :: rk = kind(1.0d0)   ! work in double precision

  ! ---- category indices -------------------------------------------
  integer, parameter :: IC=1, IR=2, II=3, IS=4, IG=5   ! cloud rain ice snow graupel
  integer, parameter :: NCAT=5

  ! ---- physical constants -----------------
  real(rk), parameter :: PI    = 3.14159265358979_rk
  real(rk), parameter :: CP    = 1004.64_rk
  real(rk), parameter :: RD    = 287.04_rk
  real(rk), parameter :: RV    = 461.5_rk
  real(rk), parameter :: P00   = 1.0e5_rk
  real(rk), parameter :: L_VAP = 2.501e6_rk
  real(rk), parameter :: L_SUB = 2.838e6_rk
  real(rk), parameter :: L_FUS = L_SUB - L_VAP
  real(rk), parameter :: E0    = 611.2_rk
  real(rk), parameter :: T0    = 273.15_rk
  real(rk), parameter :: T_FREEZE = T0
  real(rk), parameter :: T_HOM = 235.15_rk
  real(rk), parameter :: RHO_W = 1000.0_rk
  real(rk), parameter :: RHO_I = 917.0_rk
  real(rk), parameter :: RHO_S = 100.0_rk
  real(rk), parameter :: RHO_G = 400.0_rk
  real(rk), parameter :: Q_MIN = 1.0e-12_rk
  real(rk), parameter :: N_MIN = 1.0_rk
  real(rk), parameter :: QC_MIN= 1.0e-8_rk
  real(rk), parameter :: QR_MIN= 1.0e-8_rk
  real(rk), parameter :: QI_MIN= 1.0e-10_rk
  real(rk), parameter :: QS_MIN= 1.0e-10_rk
  real(rk), parameter :: QG_MIN= 1.0e-10_rk

  ! ---- diffusion / vapour constants ---------------------
  real(rk), parameter :: D_V  = 2.26e-5_rk
  real(rk), parameter :: K_A  = 2.43e-2_rk
  real(rk), parameter :: NU_AIR = 1.5e-5_rk
  real(rk), parameter :: M_ICE_MIN = 1.0e-12_rk

  ! Deposition ventilation mass-weighted Gamma ratio, matching the
  ! Only ice/snow/graupel deposit; cloud/rain entries unused.
  !                            cloud   rain   ice     snow      graupel
  real(rk), parameter :: vent_Gq(NCAT) = (/ 0.0_rk, 0.0_rk, 1.0_rk, 0.8873_rk, 0.8912_rk /)

  ! ---- PSD mass-diameter params: m = am*D^bm,  shape mu ----------
  !                          cloud      rain     ice      snow     graupel
  real(rk), parameter :: mu(NCAT) = (/ 3.0_rk, 0.0_rk, 0.0_rk, 0.0_rk, 0.0_rk /)
  real(rk), parameter :: bm(NCAT) = (/ 3.0_rk, 3.0_rk, 3.0_rk, 2.0_rk, 3.0_rk /)
  real(rk), parameter :: am(NCAT) = (/ PI/6.0_rk*RHO_W, PI/6.0_rk*RHO_W, &
                                       PI/6.0_rk*RHO_I, 0.069_rk,         &
                                       PI/6.0_rk*RHO_G /)
  real(rk), parameter :: N0_graupel = 4.0e4_rk

  ! ---- fall-speed params: V = av*D^bv ---------
  real(rk), parameter :: av(NCAT) = (/ 3.0e7_rk, 841.99_rk, 700.0_rk, 11.72_rk, 19.3_rk /)
  real(rk), parameter :: bv(NCAT) = (/ 2.0_rk,   0.8_rk,    1.0_rk,   0.41_rk,  0.37_rk /)

  ! ---- gamma factors, filled by mp_microwcp_init ---------------
  real(rk) :: Gb1(NCAT)   ! gamma(mu+bm+1)   mass-moment norm
  real(rk) :: G1(NCAT)    ! gamma(mu+1)      number-moment norm
  real(rk) :: Gq(NCAT)    ! gamma(mu+bv+1)/gamma(mu+1)   mass-weighted V
  real(rk) :: GN(NCAT)    ! gamma(mu+bv+1)/gamma(mu+2)   number-weighted V
  real(rk) :: G6(NCAT)    ! gamma(mu+7)      6th-moment (reflectivity)
  logical  :: initialized = .false.

  ! ---- diagnostics constants (reflectivity / effective radii) ----
  real(rk), parameter :: Nt_c   = 100.0e6_rk   ! fixed cloud number when no aerosol [m-3]
  real(rk), parameter :: K2_ICE = 0.176_rk     ! ice dielectric factor |K|^2
  real(rk), parameter :: K2_LIQ = 0.93_rk      ! liquid dielectric factor |K|^2
  !                              cloud   rain    ice     snow    graupel
  real(rk), parameter :: rho_cat(NCAT) = (/ RHO_W, RHO_W, RHO_I, RHO_S, RHO_G /)
  ! WRF effective-radius bounds [m] (Thompson-like clamps for radiation)
  real(rk), parameter :: REQC_MIN=2.51e-6_rk, REQC_MAX=50.e-6_rk
  real(rk), parameter :: REQI_MIN=4.99e-6_rk, REQI_MAX=125.e-6_rk
  real(rk), parameter :: REQS_MIN=4.99e-6_rk, REQS_MAX=999.e-6_rk

contains

  !=================================================================
  ! mp_microwcp_init
  !   Precompute the per-category Gamma-function ratios that appear in
  !   the PSD moment relations and fall-speed integrals. Call once before
  !   the first column (the driver and the column routine call it lazily).
  !   In a production port these become compile-time PARAMETER literals,
  !   so no Gamma function is evaluated at run time.
  !     Gb1 = gamma(mu+bm+1)              (mass-moment normalisation)
  !     G1  = gamma(mu+1)                 (number-moment normalisation)
  !     Gq  = gamma(mu+bv+1)/gamma(mu+1)  (mass-weighted bulk fall speed)
  !     GN  = gamma(mu+bv+1)/gamma(mu+2)  (number-weighted bulk fall speed)
  !=================================================================
  subroutine mp_microwcp_init()
    integer :: c
    do c = 1, NCAT
      Gb1(c) = gamma(mu(c) + bm(c) + 1.0_rk)
      G1(c)  = gamma(mu(c) + 1.0_rk)
      ! mass-weighted bulk fall speed: weight V(D)=av*D^bv by mass m(D)=am*D^bm,
      ! so the gamma ratio must include the mass-diameter exponent bm.
      ! (Previously used gamma(mu+bv+1)/gamma(mu+1), the NUMBER-weighted factor,
      !  which made q sediment ~3x too slow for rain -> surface accumulation.)
      Gq(c)  = gamma(mu(c) + bm(c) + bv(c) + 1.0_rk) / gamma(mu(c) + bm(c) + 1.0_rk)
      GN(c)  = gamma(mu(c) + bv(c) + 1.0_rk) / gamma(mu(c) + 1.0_rk)
      G6(c)  = gamma(mu(c) + 7.0_rk)
    end do
    initialized = .true.
  end subroutine mp_microwcp_init

  !=================================================================
  ! mp_microwcp_aero_init
  !   Fill a default background aerosol profile (water-friendly nwfa,
  !   ice-friendly nifa) when the input fields are essentially zero, so
  !   the scheme has CCN/INP to nucleate from even when real/met_em
  !   provides none. Mirrors Thompson's thompson_init exactly (same
  !   exponential, boundary-layer-following profile and constants), so
  !   MICROWCP and Thompson start from the same aerosol field.
  !   Default values are used per gridpoint; existing (nonzero) input is
  !   left untouched. Arrays are WRF default REAL (single precision).
  !=================================================================
  subroutine mp_microwcp_aero_init(hgt, nwfa, nifa,                    &
       ids,ide,jds,jde, kds,kde, ims,ime,jms,jme, kms,kme,             &
       its,ite,jts,jte, kts,kte)
    integer, intent(in) :: ids,ide,jds,jde,kds,kde
    integer, intent(in) :: ims,ime,jms,jme,kms,kme
    integer, intent(in) :: its,ite,jts,jte,kts,kte
    real, intent(in)    :: hgt(ims:ime,kms:kme,jms:jme)
    real, intent(inout) :: nwfa(ims:ime,kms:kme,jms:jme)
    real, intent(inout) :: nifa(ims:ime,kms:kme,jms:jme)
    ! background aerosol constants (== module_mp_thompson.F)
    real, parameter :: naCCN0=300.0e6, naCCN1=50.0e6
    real, parameter :: naIN0 =  1.5e6, naIN1 = 0.5e6
    real, parameter :: aero_eps = 1.0e-15
    real    :: h_01, niCCN3, niIN3
    integer :: i,j,k

    ! ---- water-friendly (CCN) ----
    if (maxval(nwfa(its:min(ide-1,ite),:,jts:min(jde-1,jte))) < aero_eps) then
      do j = jts, min(jde-1,jte)
      do i = its, min(ide-1,ite)
        if      (hgt(i,1,j) <= 1000.0) then ; h_01 = 0.8
        else if (hgt(i,1,j) >= 2500.0) then ; h_01 = 0.01
        else ; h_01 = 0.8*cos(hgt(i,1,j)*0.001 - 1.0) ; end if
        niCCN3 = -1.0*log(naCCN1/naCCN0)/h_01
        do k = kts, kte
          nwfa(i,k,j) = naCCN1 + naCCN0*exp(-((hgt(i,k,j)-hgt(i,1,j))/1000.)*niCCN3)
        end do
      end do
      end do
    end if

    ! ---- ice-friendly (INP) ----
    if (maxval(nifa(its:min(ide-1,ite),:,jts:min(jde-1,jte))) < aero_eps) then
      do j = jts, min(jde-1,jte)
      do i = its, min(ide-1,ite)
        if      (hgt(i,1,j) <= 1000.0) then ; h_01 = 0.8
        else if (hgt(i,1,j) >= 2500.0) then ; h_01 = 0.01
        else ; h_01 = 0.8*cos(hgt(i,1,j)*0.001 - 1.0) ; end if
        niIN3 = -1.0*log(naIN1/naIN0)/h_01
        do k = kts, kte
          nifa(i,k,j) = naIN1 + naIN0*exp(-((hgt(i,k,j)-hgt(i,1,j))/1000.)*niIN3)
        end do
      end do
      end do
    end if
  end subroutine mp_microwcp_aero_init

  !=================================================================
  ! Thermodynamics (== micro/thermo.py).  All elemental: act on scalars
  ! or whole arrays.
  !=================================================================

  ! Saturation vapour pressure over liquid water [Pa]  (Bolton 1980).
  elemental real(rk) function es_liquid(T)
    real(rk), intent(in) :: T
    es_liquid = E0 * exp(17.67_rk*(T-T0)/(T-29.65_rk))
  end function

  ! Saturation vapour pressure over ice [Pa]  (Magnus form over ice).
  elemental real(rk) function es_ice(T)
    real(rk), intent(in) :: T
    es_ice = E0 * exp(22.46_rk*(T-T0)/(T-0.53_rk))
  end function

  ! Saturation mixing ratio over liquid [kg kg-1] from T,P.
  elemental real(rk) function qsat_liquid(T,P)
    real(rk), intent(in) :: T, P
    real(rk) :: es
    es = es_liquid(T)
    qsat_liquid = 0.622_rk*es/max(P-es, 1.0_rk)
  end function

  ! Saturation mixing ratio over ice [kg kg-1] from T,P.
  elemental real(rk) function qsat_ice(T,P)
    real(rk), intent(in) :: T, P
    real(rk) :: es
    es = es_ice(T)
    qsat_ice = 0.622_rk*es/max(P-es, 1.0_rk)
  end function

  !=================================================================
  ! PSD diagnostics 
  !=================================================================

  ! Temperature-dependent snow intercept N0 [m-4] (Houze et al. 1979),
  ! used because snow is 1-moment (its N0 is prescribed, not predicted).
  elemental real(rk) function snow_N0(T)
    real(rk), intent(in) :: T
    real(rk) :: Tc
    Tc = min(T-T0, -0.1_rk)        ! clamp to below freezing
    snow_N0 = 2.0e6_rk * exp(-0.12_rk*Tc)
  end function

  ! Mass-weighted bulk fall speed [m s-1] for category c:  V_q = av*Gq/lam^bv.
  elemental real(rk) function fall_speed_q(lam, c)
    real(rk), intent(in) :: lam
    integer,  intent(in) :: c
    fall_speed_q = av(c)*Gq(c) / lam**bv(c)
  end function

  ! Number-weighted bulk fall speed [m s-1] for category c:  V_N = av*GN/lam^bv.
  elemental real(rk) function fall_speed_N(lam, c)
    real(rk), intent(in) :: lam
    integer,  intent(in) :: c
    fall_speed_N = av(c)*GN(c) / lam**bv(c)
  end function

  !-----------------------------------------------------------------
  ! diagnose
  !   Given the prognostic fields of one column, compute the PSD
  !   descriptors every process needs: slope lam, mean diameter Dm, and
  !   bulk number Nbulk for all 5 categories (plus the snow/graupel
  !   intercepts N0s,N0g). Cloud/rain/ice are 2-moment (lam from number &
  !   mass); snow/graupel are 1-moment (lam from prescribed N0 & mass,
  !   then number diagnosed). Called fresh before each process group so
  !   the sizes stay consistent with the updated state.
  !-----------------------------------------------------------------
  subroutine diagnose(nz, T, rho, qc,qr,qi,qs,qg, nc,nr,ni, &
                      lam, Dm, Nbulk, N0s, N0g)
    integer,  intent(in)  :: nz
    real(rk), intent(in)  :: T(nz), rho(nz)
    real(rk), intent(in)  :: qc(nz),qr(nz),qi(nz),qs(nz),qg(nz)
    real(rk), intent(in)  :: nc(nz),nr(nz),ni(nz)
    real(rk), intent(out) :: lam(nz,NCAT), Dm(nz,NCAT), Nbulk(nz,NCAT)
    real(rk), intent(out) :: N0s(nz), N0g(nz)
    integer  :: k
    real(rk) :: q, Nb, ratio, N0

    do k = 1, nz
      ! ---- 2-moment: cloud, rain, ice  (lam from N_bulk, q) ----
      call lam_from_Nq(qc(k), nc(k)*rho(k), rho(k), IC, lam(k,IC), Nbulk(k,IC))
      call lam_from_Nq(qr(k), nr(k)*rho(k), rho(k), IR, lam(k,IR), Nbulk(k,IR))
      call lam_from_Nq(qi(k), ni(k)*rho(k), rho(k), II, lam(k,II), Nbulk(k,II))
      ! ---- 1-moment: snow, graupel  (N0 first, then lam, then N) ----
      N0s(k) = snow_N0(T(k))
      N0g(k) = N0_graupel
      lam(k,IS) = lam_from_N0q(qs(k), N0s(k), rho(k), IS)
      lam(k,IG) = lam_from_N0q(qg(k), N0g(k), rho(k), IG)
      Nbulk(k,IS) = N0s(k)*G1(IS)/lam(k,IS)**(mu(IS)+1.0_rk)
      Nbulk(k,IG) = N0g(k)*G1(IG)/lam(k,IG)**(mu(IG)+1.0_rk)
      ! ---- mean diameters ----
      Dm(k,IC) = (mu(IC)+bm(IC)+1.0_rk)/lam(k,IC)
      Dm(k,IR) = (mu(IR)+bm(IR)+1.0_rk)/lam(k,IR)
      Dm(k,II) = (mu(II)+bm(II)+1.0_rk)/lam(k,II)
      Dm(k,IS) = (mu(IS)+bm(IS)+1.0_rk)/lam(k,IS)
      Dm(k,IG) = (mu(IG)+bm(IG)+1.0_rk)/lam(k,IG)
    end do
  contains
    ! 2-moment slope: solve for lam from bulk number Nbin and mass qin.
    ! lam = [ am*N*gamma(mu+bm+1) / (rho*q*gamma(mu+1)) ] ** (1/bm)
    subroutine lam_from_Nq(qin, Nbin, rhoin, c, lamout, Nbout)
      real(rk), intent(in)  :: qin, Nbin, rhoin
      integer,  intent(in)  :: c
      real(rk), intent(out) :: lamout, Nbout
      real(rk) :: qs_, Ns_
      qs_ = max(qin, Q_MIN)
      Ns_ = max(Nbin, N_MIN)
      Nbout = Nbin              ! pass the bulk number straight through
      ratio = (am(c)*Ns_*Gb1(c)) / (rhoin*qs_*G1(c))
      lamout = ratio**(1.0_rk/bm(c))
    end subroutine
    ! 1-moment slope: solve for lam from prescribed intercept N0in and mass.
    ! lam = [ N0*am*gamma(mu+bm+1) / (rho*q) ] ** (1/(mu+bm+1))
    real(rk) function lam_from_N0q(qin, N0in, rhoin, c)
      real(rk), intent(in) :: qin, N0in, rhoin
      integer,  intent(in) :: c
      real(rk) :: qs_
      qs_ = max(qin, Q_MIN)
      lam_from_N0q = (N0in*am(c)*Gb1(c)/(rhoin*qs_))**(1.0_rk/(mu(c)+bm(c)+1.0_rk))
    end function
  end subroutine diagnose

  !=================================================================
  ! upstream_flux
  !   First-order upwind vertical settling of one quantity q falling at
  !   speed V, returned as a time-averaged tendency dq [unit s-1].
  !   Index 1 = lowest model level; flux is downward, so the bottom layer
  !   loses mass to the surface (that loss is what becomes precip). The
  !   explicit scheme is only CFL-stable for Courant<=1, so the routine
  !   sub-steps internally (nsub) using this column's own V_max.
  !   == micro/sedimentation.py::_upstream_flux
  !=================================================================
  subroutine upstream_flux(nz, q, V, rho, dz, dt, dq)
    integer,  intent(in)  :: nz
    real(rk), intent(in)  :: q(nz), V(nz), rho(nz), dz(nz), dt
    real(rk), intent(out) :: dq(nz)
    real(rk) :: qw(nz), q0(nz), flux(nz), dqs(nz)
    real(rk) :: Vmax, dzmin, dt_sub, fin, fout
    integer  :: nsub, s, k
    qw = q; q0 = q
    Vmax = maxval(V); dzmin = minval(dz)
    if (Vmax > 0.0_rk) then
      dt_sub = dzmin/Vmax
      nsub = max(1, ceiling(dt/dt_sub))
      dt_sub = dt/real(nsub, rk)
    else
      nsub = 1; dt_sub = dt
    end if
    do s = 1, nsub
      do k = 1, nz
        flux(k) = rho(k)*V(k)*qw(k)
      end do
      do k = 1, nz
        fout = flux(k)
        if (k < nz) then
          fin = flux(k+1)
        else
          fin = 0.0_rk
        end if
        dqs(k) = (fin - fout)/max(rho(k)*dz(k), 1.0e-10_rk)
      end do
      do k = 1, nz
        qw(k) = max(qw(k) + dqs(k)*dt_sub, 0.0_rk)
      end do
    end do
    dq = (qw - q0)/dt
  end subroutine upstream_flux

  !=================================================================
  ! mp_microwcp_column   (<-> Thompson's mp_thompson)
  !   Advance ONE atmospheric column by one time step dt. This is the
  !   heart of the scheme: it runs the operator-split sequence
  !     1 nucleation -> 2 warm -> 3 cold -> 4 collection
  !     -> 5 sedimentation -> 6 saturation adjustment
  !   re-diagnosing the PSD (diagnose) before each group so particle sizes
  !   stay consistent with the updated state. Prognostic fields are 1-D
  !   (length nz, index 1 = lowest level) and updated in place. Surface
  !   precip is returned per category as depth [mm].
  !   Temperature is carried as FULL T here; the th<->T conversion lives
  !   in the driver (matching Thompson).
  !   == micro/column.py::step_column + micro/stepper.py::step
  !=================================================================
  subroutine mp_microwcp_column(nz, T, P, rho, w, dz, &
       qv,qc,qr,qi,qs,qg, nc,nr,ni, nwfa,nifa, has_aero, dt, &
       pptrain, pptsnow, pptgraul, pptice)
    integer,  intent(in)    :: nz
    real(rk), intent(inout) :: T(nz)
    real(rk), intent(in)    :: P(nz), rho(nz), w(nz), dz(nz)
    real(rk), intent(inout) :: qv(nz),qc(nz),qr(nz),qi(nz),qs(nz),qg(nz)
    real(rk), intent(inout) :: nc(nz),nr(nz),ni(nz),nwfa(nz),nifa(nz)
    logical,  intent(in)    :: has_aero
    real(rk), intent(in)    :: dt
    real(rk), intent(out)   :: pptrain, pptsnow, pptgraul, pptice

    real(rk) :: lam(nz,NCAT), Dm(nz,NCAT), Nb(nz,NCAT), N0s(nz), N0g(nz)
    real(rk) :: dqr_sed(nz), dqi_sed(nz), dqs_sed(nz), dqg_sed(nz)

    if (.not. initialized) call mp_microwcp_init()

    ! ---- 1. nucleation ----
    if (has_aero) call nucleation(nz, T,P,rho, qv,qc,qi, nc,ni, nwfa,nifa, dt)

    ! ---- 2. warm phase ----
    call diagnose(nz,T,rho, qc,qr,qi,qs,qg, nc,nr,ni, lam,Dm,Nb,N0s,N0g)
    call warm(nz, T,P,rho, qv,qc,qr, nc,nr, Dm, dt)

    ! ---- 3. cold phase ----
    call diagnose(nz,T,rho, qc,qr,qi,qs,qg, nc,nr,ni, lam,Dm,Nb,N0s,N0g)
    call cold(nz, T,P,rho, qv,qc,qr,qi,qs,qg, nr,ni, lam,Dm,Nb, dt)

    ! ---- 4. collection ----
    call diagnose(nz,T,rho, qc,qr,qi,qs,qg, nc,nr,ni, lam,Dm,Nb,N0s,N0g)
    call collection(nz, T,rho, qc,qr,qi,qs,qg, nc,nr,ni, lam,Dm,Nb, dt)

    ! ---- 5. sedimentation ----
    call diagnose(nz,T,rho, qc,qr,qi,qs,qg, nc,nr,ni, lam,Dm,Nb,N0s,N0g)
    call sedimentation(nz, rho,dz,dt, qc,qr,qi,qs,qg, nc,nr,ni, lam, &
                       dqr_sed,dqi_sed,dqs_sed,dqg_sed)

    ! ---- 6. saturation adjustment ----
    call saturation_adjustment(nz, T,P, qv,qc,qi,qs,qg)

    ! ---- surface precip [kg m-2 == mm] = -integral(dq_sed*rho*dz)*dt ----
    pptrain  = max(-sum(dqr_sed*rho*dz)*dt, 0.0_rk)
    pptsnow  = max(-sum(dqs_sed*rho*dz)*dt, 0.0_rk)
    pptgraul = max(-sum(dqg_sed*rho*dz)*dt, 0.0_rk)
    pptice   = max(-sum(dqi_sed*rho*dz)*dt, 0.0_rk)
  end subroutine mp_microwcp_column

  !-----------------------------------------------------------------
  !-----------------------------------------------------------------
  ! nucleation  
  !   Birth of new hydrometeors from vapour/aerosol:
  !     - CCN activation (Twomey): vapour supersaturation + water-friendly
  !       aerosol nwfa -> new cloud droplets (nc,qc).
  !     - Heterogeneous ice nucleation (DeMott 2010): ice-friendly aerosol
  !       nifa, between -38 and 0 C -> new ice crystals (ni,qi).
  !     - Homogeneous freezing: below -38 C all cloud instantly -> ice.
  !   Seed condensate is drawn from vapour (qv) so water is conserved.
  !   Only runs when aerosols are present (has_aero in the column routine).
  !-----------------------------------------------------------------
  subroutine nucleation(nz, T,P,rho, qv,qc,qi, nc,ni, nwfa,nifa, dt)
    integer,  intent(in)    :: nz
    real(rk), intent(in)    :: P(nz),rho(nz),nwfa(nz),nifa(nz),dt
    real(rk), intent(inout) :: T(nz),qv(nz),qc(nz),qi(nz),nc(nz),ni(nz)
    ! All sub-process rates are evaluated from the FROZEN input state and
    ! accumulated, then applied once (matches stepper.step's tendency sum).
    integer  :: k
    real(rk) :: T0k,qv0,qc0,qi0,nc0,ni0
    real(rk) :: aqv,aqc,aqi,anc,ani,aT
    real(rk) :: qsat,S,Nact,dNc,m_drop,dqc
    real(rk) :: na_cm3,Tsc,Ninp_L,Ninp,dNi,dqi
    real(rk) :: a_dm,b_dm,c_dm,d_dm
    a_dm=0.0000594_rk; b_dm=3.33_rk; c_dm=0.0264_rk; d_dm=0.0827_rk
    m_drop = (4.0_rk/3.0_rk)*PI*(1.0e-6_rk)**3*1000.0_rk
    do k = 1, nz
      T0k=T(k); qv0=qv(k); qc0=qc(k); qi0=qi(k); nc0=nc(k); ni0=ni(k)
      aqv=0.0_rk; aqc=0.0_rk; aqi=0.0_rk; anc=0.0_rk; ani=0.0_rk; aT=0.0_rk
      ! --- CCN activation (Twomey) ---
      qsat = qsat_liquid(T0k,P(k))
      S = max((qv0-qsat)/max(qsat,Q_MIN), 0.0_rk)
      if (S > 0.0_rk) then
        Nact = min(nwfa(k)*S**0.462_rk, nwfa(k))
        dNc  = max(Nact-nc0,0.0_rk)/dt
        dqc  = dNc*m_drop
        anc=anc+dNc; aqc=aqc+dqc; aqv=aqv-dqc; aT=aT+L_VAP/CP*dqc
      end if
      ! --- heterogeneous ice nucleation (DeMott 2010) ---
      Tsc = max(T0+0.01_rk - T0k, 0.0_rk)            ! 273.16 - T
      if (T0k > T_HOM .and. T0k < T_FREEZE .and. Tsc > 0.0_rk) then
        na_cm3 = max(nifa(k)*rho(k)*1.0e-6_rk, 1.0e-6_rk)
        Ninp_L = a_dm*(Tsc**b_dm)*(na_cm3**(c_dm*Tsc+d_dm))     ! the actual parameterization
        Ninp   = Ninp_L*1.0e3_rk/max(rho(k),1.0e-3_rk)
        dNi    = max(Ninp-ni0,0.0_rk)/dt
        dqi    = dNi*M_ICE_MIN
        ani=ani+dNi; aqi=aqi+dqi; aqv=aqv-dqi; aT=aT+L_SUB/CP*dqi
      end if
      ! --- homogeneous freezing (qc->qi, nc->ni) ---
      if (T0k < T_HOM) then
        dqi = qc0/dt; dNi = nc0/dt
        aqi=aqi+dqi; ani=ani+dNi; aqc=aqc-qc0/dt; anc=anc-nc0/dt
        aT=aT+L_FUS/CP*dqi
      end if
      ! --- apply once ---
      qc(k)=max(qc0+aqc*dt,0.0_rk); qv(k)=max(qv0+aqv*dt,0.0_rk)
      qi(k)=max(qi0+aqi*dt,0.0_rk)
      nc(k)=max(nc0+anc*dt,0.0_rk); ni(k)=max(ni0+ani*dt,0.0_rk)
      T(k) =T0k+aT*dt
    end do
  end subroutine nucleation

  !-----------------------------------------------------------------
  !-----------------------------------------------------------------
  ! warm  
  !   Liquid-only processes that make and shape rain:
  !     - autoconversion (KK2000): cloud droplets coalesce into rain.
  !     - accretion (KK2000): rain sweeps up cloud droplets.
  !     - rain evaporation: rain -> vapour in subsaturated air (cools).
  !     - self-collection / breakup: adjust rain DROP NUMBER nr only
  !       (mass conserved), driving the spectrum toward equilibrium size.
  !-----------------------------------------------------------------
  subroutine warm(nz, T,P,rho, qv,qc,qr, nc,nr, Dm, dt)
    integer,  intent(in)    :: nz
    real(rk), intent(in)    :: P(nz),rho(nz),Dm(nz,NCAT),dt
    real(rk), intent(inout) :: T(nz),qv(nz),qc(nz),qr(nz),nc(nz),nr(nz)
    integer  :: k
    real(rk) :: qv0,qc0,qr0,nc0,nr0
    real(rk) :: aqv,aqc,aqr,anc,anr,aT
    real(rk) :: Nc_cm3,dqr,D0,m_drop,dNr,Dc,n_coal
    real(rk) :: qsat,sub,Cev,dqr_e
    D0=50.0e-6_rk; Dc=15.0e-6_rk
    m_drop=(PI/6.0_rk)*RHO_W*D0**3; n_coal=(D0/Dc)**3
    do k = 1, nz
      qv0=qv(k); qc0=qc(k); qr0=qr(k); nc0=nc(k); nr0=nr(k)
      aqv=0.0_rk; aqc=0.0_rk; aqr=0.0_rk; anc=0.0_rk; anr=0.0_rk; aT=0.0_rk
      ! --- autoconversion ---
      if (qc0 > QC_MIN) then
        Nc_cm3 = max(nc0*rho(k)*1.0e-6_rk, 1.0e-6_rk)
        dqr = min(1350.0_rk*qc0**2.47_rk*Nc_cm3**(-1.79_rk), qc0/dt)
        dNr = dqr/m_drop
        aqr=aqr+dqr; aqc=aqc-dqr; anr=anr+dNr; anc=anc-dNr*n_coal
      end if
      ! --- accretion ---
      if (qc0 > QC_MIN .and. qr0 > QR_MIN) then
        dqr = min(67.0_rk*(qc0*qr0)**1.15_rk, qc0/dt)
        aqr=aqr+dqr; aqc=aqc-dqr
      end if
      ! --- rain evaporation ---
      if (qr0 > QR_MIN) then
        qsat = qsat_liquid(T(k),P(k))
        sub  = max(1.0_rk - qv0/max(qsat,Q_MIN), 0.0_rk)
        Cev  = 1.6_rk + 124.9_rk*(rho(k)*qr0)**0.2046_rk
        dqr_e= min(Cev*sub*qr0**0.525_rk, qr0/dt)
        aqr=aqr-dqr_e; aqv=aqv+dqr_e
        anr=anr - nr0*dqr_e/max(qr0,Q_MIN); aT=aT - L_VAP/CP*dqr_e
      end if
      ! --- rain self-collection ---
      if (qr0 > QR_MIN) then
        anr = anr + max(-5.78e-4_rk*(nr0*rho(k))**2/rho(k), -nr0/dt)
      end if
      ! --- rain breakup ---
      if (qr0 > QR_MIN) then
        anr = anr + min(1000.0_rk*max(Dm(k,IR)-0.9e-3_rk,0.0_rk)*nr0, nr0/dt)
      end if
      ! --- apply once ---
      qc(k)=max(qc0+aqc*dt,0.0_rk); qr(k)=max(qr0+aqr*dt,0.0_rk)
      qv(k)=max(qv0+aqv*dt,0.0_rk)
      nc(k)=max(nc0+anc*dt,0.0_rk); nr(k)=max(nr0+anr*dt,0.0_rk)
      T(k) =T(k)+aT*dt
    end do
  end subroutine warm

  !-----------------------------------------------------------------
  !-----------------------------------------------------------------
  ! cold  
  !   Ice-phase processes that exchange mass with vapour and between ice
  !   species:
  !     - vapour deposition/sublimation on ice, snow, graupel (qv<->qx),
  !       diffusion-limited with a ventilation correction (dep_rate).
  !     - rain freezing (Bigg): supercooled rain -> graupel.
  !     - snow aggregation: cloud-ice -> snow.
  !     - Hallett-Mossop: rime-splintering source of secondary ice
  !       (active only -8..-3 C with cloud present).
  !     - melting: snow/graupel -> rain above 0 C (cools).
  !-----------------------------------------------------------------
  subroutine cold(nz, T,P,rho, qv,qc,qr,qi,qs,qg, nr,ni, lam,Dm,Nb, dt)
    integer,  intent(in)    :: nz
    real(rk), intent(in)    :: P(nz),rho(nz),lam(nz,NCAT),Dm(nz,NCAT),Nb(nz,NCAT),dt
    real(rk), intent(inout) :: T(nz),qv(nz),qc(nz),qr(nz),qi(nz),qs(nz),qg(nz),nr(nz),ni(nz)
    integer  :: k
    real(rk) :: T0k,qv0,qc0,qr0,qi0,qs0,qg0,nr0,ni0
    real(rk) :: aqv,aqc,aqr,aqi,aqs,aqg,anr,ani,aT
    real(rk) :: dqi,dqs,dqg,dqr_f,dqs_a,dNi_h,dqi_h
    real(rk) :: fT,rim,coeff,rate_s,rate_g
    real(rk) :: K_bigg, dTf
    K_bigg=(PI/6.0_rk)*RHO_W*100.0_rk
    do k = 1, nz
      T0k=T(k); qv0=qv(k); qc0=qc(k); qr0=qr(k)
      qi0=qi(k); qs0=qs(k); qg0=qg(k); nr0=nr(k); ni0=ni(k)
      aqv=0.0_rk; aqc=0.0_rk; aqr=0.0_rk; aqi=0.0_rk; aqs=0.0_rk
      aqg=0.0_rk; anr=0.0_rk; ani=0.0_rk; aT=0.0_rk
      ! --- deposition/sublimation on ice/snow/graupel ---
      dqi = dep_rate(qi0,Nb(k,II),Dm(k,II),lam(k,II),qv0,T0k,P(k),rho(k),QI_MIN,dt,II)
      if (dqi < -1.0e-8_rk) ani = ani + ni0*dqi/max(qi0,Q_MIN)
      aqi=aqi+dqi; aqv=aqv-dqi; aT=aT+L_SUB/CP*dqi
      dqs = dep_rate(qs0,Nb(k,IS),Dm(k,IS),lam(k,IS),qv0,T0k,P(k),rho(k),QS_MIN,dt,IS)
      aqs=aqs+dqs; aqv=aqv-dqs; aT=aT+L_SUB/CP*dqs
      dqg = dep_rate(qg0,Nb(k,IG),Dm(k,IG),lam(k,IG),qv0,T0k,P(k),rho(k),QG_MIN,dt,IG)
      aqg=aqg+dqg; aqv=aqv-dqg; aT=aT+L_SUB/CP*dqg
      ! --- rain freezing (Bigg) -> graupel ---
      if (T0k < T_FREEZE .and. qr0 > QR_MIN) then
        dTf = max(T_FREEZE-T0k,0.0_rk)
        dqr_f = min(K_bigg*exp(0.66_rk*dTf)*nr0*Dm(k,IR)**3, qr0/dt)
        aqr=aqr-dqr_f; aqg=aqg+dqr_f
        anr=anr - nr0*dqr_f/max(qr0,Q_MIN); aT=aT+L_FUS/CP*dqr_f
      end if
      ! --- snow aggregation: ice -> snow ---
      if (T0k < T_FREEZE .and. qi0 > QI_MIN .and. qs0 > QS_MIN) then
        dqs_a = min(0.1_rk*0.25_rk*Nb(k,II)*qi0, qi0/dt)
        aqs=aqs+dqs_a; aqi=aqi-dqs_a
        ani=ani - ni0*dqs_a/max(qi0,Q_MIN)
      end if
      ! --- Hallett-Mossop rime splintering ---
      if (T0k > 265.15_rk .and. T0k < 270.15_rk .and. qc0 > QC_MIN) then
        fT = 1.0_rk - abs(T0k-268.15_rk)/(270.15_rk-268.15_rk)
        fT = max(min(fT,1.0_rk),0.0_rk)
        rim = (PI/4.0_rk)*1.0_rk*qc0*(Dm(k,IS)**2*Nb(k,IS)+Dm(k,IG)**2*Nb(k,IG))
        dNi_h = 3.5e8_rk*rim*fT
        dqi_h = dNi_h*M_ICE_MIN
        ani=ani+dNi_h; aqi=aqi+dqi_h
      end if
      ! --- melting: snow/graupel -> rain (T>0) ---
      if (T0k > T_FREEZE) then
        coeff = 2.0_rk*PI*K_A/L_FUS
        rate_s=0.0_rk; rate_g=0.0_rk
        if (qs0 > QS_MIN) rate_s = min(coeff*(T0k-T_FREEZE)*Nb(k,IS)*Dm(k,IS)/(3.0_rk*rho(k)), qs0/dt)
        if (qg0 > QG_MIN) rate_g = min(coeff*(T0k-T_FREEZE)*Nb(k,IG)*Dm(k,IG)/(3.0_rk*rho(k)), qg0/dt)
        aqs=aqs-rate_s; aqg=aqg-rate_g; aqr=aqr+(rate_s+rate_g)
        aT=aT - L_FUS/CP*(rate_s+rate_g)
      end if
      ! --- apply once ---
      qv(k)=max(qv0+aqv*dt,0.0_rk); qc(k)=max(qc0+aqc*dt,0.0_rk)
      qr(k)=max(qr0+aqr*dt,0.0_rk); qi(k)=max(qi0+aqi*dt,0.0_rk)
      qs(k)=max(qs0+aqs*dt,0.0_rk); qg(k)=max(qg0+aqg*dt,0.0_rk)
      nr(k)=max(nr0+anr*dt,0.0_rk); ni(k)=max(ni0+ani*dt,0.0_rk)
      T(k) =T0k+aT*dt
    end do
  contains
    ! dep_rate: diffusion-limited deposition(+)/sublimation(-) rate
    ! [kg kg-1 s-1] for one ice category c. Combines supersaturation over
    ! ice (Si-1), the particle population (Nbx,Dmx), and a ventilation
    ! enhancement fv from the fall-speed Reynolds number. Capped so it
    ! cannot remove more ice than present (sublimation) or more vapour
    ! than present (deposition) in one step.
    real(rk) function dep_rate(qx,Nbx,Dmx,lamx,qvx,Tx,Px,rhx,qmin,dtl,c)
      real(rk), intent(in) :: qx,Nbx,Dmx,lamx,qvx,Tx,Px,rhx,qmin,dtl
      integer,  intent(in) :: c
      real(rk) :: Vf,Re,xv,fv,Si,esi,Ad,Bd,den,r
      Vf = av(c)*vent_Gq(c)/max(lamx,1.0_rk)**bv(c)
      Re = max(Vf*Dmx,0.0_rk)/NU_AIR
      xv = sqrt(Re)*(NU_AIR/D_V)**(1.0_rk/3.0_rk)
      if (xv <= 1.0_rk) then
        fv = 1.0_rk + 0.14_rk*xv**2
      else
        fv = 0.78_rk + 0.308_rk*xv
      end if
      Si  = qvx/max(qsat_ice(Tx,Px),Q_MIN)
      esi = es_ice(Tx)
      ! psychrometric denominator (Pruppacher & Klett 1997, Eq. 13-28):
      !   Ad = latent-heat conduction term ; Bd = vapour diffusion term
      Ad  = L_SUB**2/max(K_A*RV*Tx**2, 1.0e-30_rk)
      Bd  = RV*Tx/max(D_V*esi, 1.0e-30_rk)
      den = rhx*(Ad+Bd)
      ! capacitance integral int((D/2) n(D) dD) = Nb*(mu+1)/(2*lam)
      r   = 4.0_rk*PI*(Si-1.0_rk)*fv &
            *(0.5_rk*Nbx*(mu(c)+1.0_rk)/max(lamx,1.0_rk)) &
            /max(den,1.0e-30_rk)
      if (qx <= qmin) r = 0.0_rk
      if (r < -1.0e-8_rk) r = max(r, -qx/dtl)
      if (r > 0.0_rk) r = min(r,  qvx/dtl)
      dep_rate = r
    end function
  end subroutine cold

  !-----------------------------------------------------------------
  !-----------------------------------------------------------------
  ! collection  
  !   Mixed-phase collision-collection (riming and wet growth), active
  !   only below 0 C. A faster-falling "collector" sweeps up a slower
  !   "collected" species; the differential fall speed sets the rate:
  !     1 cloud+ice  -> snow      2 cloud+snow -> graupel
  !     3 cloud+graupel -> graupel 4 rain+ice  -> graupel
  !     5 snow+rain  -> graupel
  !   Each freezes liquid onto ice, so it warms the air (L_fus).
  !-----------------------------------------------------------------
  subroutine collection(nz, T,rho, qc,qr,qi,qs,qg, nc,nr,ni, lam,Dm,Nb, dt)
    integer,  intent(in)    :: nz
    real(rk), intent(in)    :: rho(nz),lam(nz,NCAT),Dm(nz,NCAT),Nb(nz,NCAT),dt
    real(rk), intent(inout) :: T(nz),qc(nz),qr(nz),qi(nz),qs(nz),qg(nz),nc(nz),nr(nz),ni(nz)
    integer  :: k
    real(rk) :: T0k,qc0,qr0,qi0,qs0,qg0,nc0,nr0,ni0
    real(rk) :: aqc,aqr,aqi,aqs,aqg,anc,anr,ani,aT
    real(rk) :: Vrel,rate,rate_qi,rate_qr
    do k = 1, nz
      if (T(k) >= T_FREEZE) cycle
      T0k=T(k); qc0=qc(k); qr0=qr(k); qi0=qi(k); qs0=qs(k); qg0=qg(k)
      nc0=nc(k); nr0=nr(k); ni0=ni(k)
      aqc=0.0_rk; aqr=0.0_rk; aqi=0.0_rk; aqs=0.0_rk; aqg=0.0_rk
      anc=0.0_rk; anr=0.0_rk; ani=0.0_rk; aT=0.0_rk
      ! --- 1. cloud-ice riming: qc+qi -> qs ---
      if (qc0 > QC_MIN .and. qi0 > QI_MIN) then
        Vrel = abs(fall_speed_q(lam(k,II),II)-fall_speed_q(lam(k,IC),IC))
        rate = min(max(0.7_rk*(PI/4.0_rk)*Dm(k,II)**2*Nb(k,II)*Vrel*qc0,0.0_rk), qc0/dt)
        anc=anc + nc0*(-rate)/max(qc0,Q_MIN)
        ani=ani + ni0*(-rate)/max(qi0,Q_MIN)
        aqc=aqc-rate; aqs=aqs+rate; aT=aT+L_FUS/CP*rate
      end if
      ! --- 2. cloud-snow riming: qc+qs -> qs (riming grows snow) ---
      if (qc0 > QC_MIN .and. qs0 > QS_MIN) then
        Vrel = abs(fall_speed_q(lam(k,IS),IS)-fall_speed_q(lam(k,IC),IC))
        rate = min(max(0.8_rk*(PI/4.0_rk)*Dm(k,IS)**2*Nb(k,IS)*Vrel*qc0,0.0_rk), qc0/dt)
        anc=anc + nc0*(-rate)/max(qc0,Q_MIN)
        aqc=aqc-rate; aqs=aqs+rate; aT=aT+L_FUS/CP*rate
      end if
      ! --- 3. cloud-graupel riming: qc+qg -> qg ---
      if (qc0 > QC_MIN .and. qg0 > QG_MIN) then
        Vrel = abs(fall_speed_q(lam(k,IG),IG)-fall_speed_q(lam(k,IC),IC))
        rate = min(max(1.0_rk*(PI/4.0_rk)*Dm(k,IG)**2*Nb(k,IG)*Vrel*qc0,0.0_rk), qc0/dt)
        anc=anc + nc0*(-rate)/max(qc0,Q_MIN)
        aqc=aqc-rate; aqg=aqg+rate; aT=aT+L_FUS/CP*rate
      end if
      ! --- 4. rain-ice collection: qr+qi -> qg ---
      if (qi0 > QI_MIN .and. qr0 > QR_MIN) then
        Vrel = abs(fall_speed_q(lam(k,IR),IR)-fall_speed_q(lam(k,II),II))
        rate_qi = min(max(1.0_rk*(PI/4.0_rk)*Dm(k,IR)**2*Nb(k,IR)*Vrel*qi0,0.0_rk), qi0/dt)
        rate_qr = min(qr0*rate_qi/max(qi0,Q_MIN), qr0/dt)
        ani=ani + ni0*(-rate_qi)/max(qi0,Q_MIN)
        anr=anr + nr0*(-rate_qr)/max(qr0,Q_MIN)
        aqi=aqi-rate_qi; aqr=aqr-rate_qr; aqg=aqg+(rate_qi+rate_qr)
        aT=aT+L_FUS/CP*rate_qr
      end if
      ! --- 5. snow collects rain: qr -> qg (rain freezes; snow is collector) ---
      if (qr0 > QR_MIN .and. qs0 > QS_MIN) then
        Vrel = abs(fall_speed_q(lam(k,IR),IR)-fall_speed_q(lam(k,IS),IS))
        rate_qr = min(max(1.0_rk*(PI/4.0_rk)*Dm(k,IS)**2*Nb(k,IS)*Vrel*qr0,0.0_rk), qr0/dt)
        anr=anr + nr0*(-rate_qr)/max(qr0,Q_MIN)
        aqr=aqr-rate_qr; aqg=aqg+rate_qr
        aT=aT+L_FUS/CP*rate_qr
      end if
      ! --- apply once ---
      qc(k)=max(qc0+aqc*dt,0.0_rk); qr(k)=max(qr0+aqr*dt,0.0_rk)
      qi(k)=max(qi0+aqi*dt,0.0_rk); qs(k)=max(qs0+aqs*dt,0.0_rk)
      qg(k)=max(qg0+aqg*dt,0.0_rk)
      nc(k)=max(nc0+anc*dt,0.0_rk); nr(k)=max(nr0+anr*dt,0.0_rk)
      ni(k)=max(ni0+ani*dt,0.0_rk); T(k)=T0k+aT*dt
    end do
  end subroutine collection

  !-----------------------------------------------------------------
  !-----------------------------------------------------------------
  ! sedimentation  
  !   Gravitational settling of every category through the column via the
  !   upwind scheme above. 2-moment species (cloud,rain,ice) settle mass
  !   and number separately (mass- vs number-weighted fall speed); snow
  !   and graupel settle mass only. Returns the per-step settling
  !   tendencies for rain/ice/snow/graupel (dq*_sed), which the column
  !   routine integrates over the column to get surface precip.
  !-----------------------------------------------------------------
  subroutine sedimentation(nz, rho,dz,dt, qc,qr,qi,qs,qg, nc,nr,ni, lam, &
                           dqr_sed,dqi_sed,dqs_sed,dqg_sed)
    integer,  intent(in)    :: nz
    real(rk), intent(in)    :: rho(nz),dz(nz),dt,lam(nz,NCAT)
    real(rk), intent(inout) :: qc(nz),qr(nz),qi(nz),qs(nz),qg(nz),nc(nz),nr(nz),ni(nz)
    real(rk), intent(out)   :: dqr_sed(nz),dqi_sed(nz),dqs_sed(nz),dqg_sed(nz)
    real(rk) :: Vq(nz),Vn(nz),dq(nz),dN(nz),Nbk(nz),dummy(nz)
    integer  :: k
    ! cloud (2-moment)
    do k=1,nz; Vq(k)=fall_speed_q(lam(k,IC),IC); Vn(k)=fall_speed_N(lam(k,IC),IC); end do
    call upstream_flux(nz,qc,Vq,rho,dz,dt,dq); qc=max(qc+dq*dt,0.0_rk)
    Nbk=nc*rho; call upstream_flux(nz,Nbk,Vn,rho,dz,dt,dN); nc=max(nc+dN/rho*dt,0.0_rk)
    ! rain
    do k=1,nz; Vq(k)=fall_speed_q(lam(k,IR),IR); Vn(k)=fall_speed_N(lam(k,IR),IR); end do
    call upstream_flux(nz,qr,Vq,rho,dz,dt,dqr_sed); qr=max(qr+dqr_sed*dt,0.0_rk)
    Nbk=nr*rho; call upstream_flux(nz,Nbk,Vn,rho,dz,dt,dN); nr=max(nr+dN/rho*dt,0.0_rk)
    ! ice
    do k=1,nz; Vq(k)=fall_speed_q(lam(k,II),II); Vn(k)=fall_speed_N(lam(k,II),II); end do
    call upstream_flux(nz,qi,Vq,rho,dz,dt,dqi_sed); qi=max(qi+dqi_sed*dt,0.0_rk)
    Nbk=ni*rho; call upstream_flux(nz,Nbk,Vn,rho,dz,dt,dN); ni=max(ni+dN/rho*dt,0.0_rk)
    ! snow (1-moment)
    do k=1,nz; Vq(k)=fall_speed_q(lam(k,IS),IS); end do
    call upstream_flux(nz,qs,Vq,rho,dz,dt,dqs_sed); qs=max(qs+dqs_sed*dt,0.0_rk)
    ! graupel (1-moment)
    do k=1,nz; Vq(k)=fall_speed_q(lam(k,IG),IG); end do
    call upstream_flux(nz,qg,Vq,rho,dz,dt,dqg_sed); qg=max(qg+dqg_sed*dt,0.0_rk)
  end subroutine sedimentation

  !-----------------------------------------------------------------
  !-----------------------------------------------------------------
  ! saturation_adjustment
  !   Final warm-phase cleanup: where T>=0 C, condense excess vapour to
  !   cloud (or evaporate cloud) so qv relaxes to qsat_liquid, updating T
  !   for the latent heat. Three Newton iterations handle qsat's own
  !   T-dependence. The ice phase is deliberately NOT adjusted here -- it
  !   is left to the diffusion-limited rates in `cold` (see satadj.py).
  !-----------------------------------------------------------------
  subroutine saturation_adjustment(nz, T,P, qv,qc,qi,qs,qg)
    integer,  intent(in)    :: nz
    real(rk), intent(in)    :: P(nz)
    real(rk), intent(inout) :: T(nz),qv(nz),qc(nz),qi(nz),qs(nz),qg(nz)
    integer  :: k, it
    real(rk) :: ql,excess,denom,dq
    do k = 1, nz
      do it = 1, 3
        if (T(k) >= T_FREEZE) then
          ql = qsat_liquid(T(k),P(k))
          excess = qv(k) - ql
          denom  = 1.0_rk + L_VAP**2*ql/(RV*CP*T(k)**2)
          dq = excess/denom
          dq = max(dq, -qc(k))
          qc(k) = qc(k) + dq
          qv(k) = qv(k) - dq
          T(k)  = T(k) + L_VAP/CP*dq
        end if
        qc(k)=max(qc(k),0.0_rk); qi(k)=max(qi(k),0.0_rk)
        qs(k)=max(qs(k),0.0_rk); qg(k)=max(qg(k),0.0_rk); qv(k)=max(qv(k),0.0_rk)
      end do
    end do
  end subroutine saturation_adjustment

  !-----------------------------------------------------------------
  ! column_diag  
  !   Optional radiation/radar diagnostics for one column:
  !     re_qc,re_qi,re_qs  effective radii [m] (clamped to WRF bounds)
  !                        Reff = (mu+3)/(2*lam); re_qs extends the
  !                        prototype (which defines only cloud & ice) with
  !                        the same form so RRTMG has a snow radius.
  !     dbz                equivalent 10-cm reflectivity [dBZ], the linear
  !                        sum of per-category 6th moments (Rayleigh).
  !-----------------------------------------------------------------
  subroutine column_diag(nz, T,rho, qc,qr,qi,qs,qg, nc,nr,ni, &
                         re_qc,re_qi,re_qs, dbz)
    integer,  intent(in)  :: nz
    real(rk), intent(in)  :: T(nz),rho(nz)
    real(rk), intent(in)  :: qc(nz),qr(nz),qi(nz),qs(nz),qg(nz),nc(nz),nr(nz),ni(nz)
    real(rk), intent(out) :: re_qc(nz),re_qi(nz),re_qs(nz),dbz(nz)
    real(rk) :: lam(nz,NCAT),Dm(nz,NCAT),Nb(nz,NCAT),N0s(nz),N0g(nz)
    real(rk) :: ze
    integer  :: k
    call diagnose(nz,T,rho, qc,qr,qi,qs,qg, nc,nr,ni, lam,Dm,Nb,N0s,N0g)
    do k = 1, nz
      re_qc(k) = min(max((mu(IC)+3.0_rk)/(2.0_rk*lam(k,IC)), REQC_MIN), REQC_MAX)
      re_qi(k) = min(max((mu(II)+3.0_rk)/(2.0_rk*lam(k,II)), REQI_MIN), REQI_MAX)
      re_qs(k) = min(max((mu(IS)+3.0_rk)/(2.0_rk*lam(k,IS)), REQS_MIN), REQS_MAX)
      ze = ze_cat(Nb(k,IC),lam(k,IC),IC,.true.)  &
         + ze_cat(Nb(k,IR),lam(k,IR),IR,.true.)  &
         + ze_cat(Nb(k,II),lam(k,II),II,.false.) &
         + ze_cat(Nb(k,IS),lam(k,IS),IS,.false.) &
         + ze_cat(Nb(k,IG),lam(k,IG),IG,.false.)
      dbz(k) = 10.0_rk*log10(max(ze, 1.0e-4_rk))
    end do
  contains
    ! linear equivalent reflectivity [mm6 m-3] for one category.
    real(rk) function ze_cat(Nbx,lamx,c,liquid)
      real(rk), intent(in) :: Nbx,lamx
      integer,  intent(in) :: c
      logical,  intent(in) :: liquid
      real(rk) :: m6
      m6 = Nbx*G6(c)/(lamx**6 * G1(c))
      if (liquid) then
        ze_cat = 1.0e18_rk*m6
      else
        ze_cat = (rho_cat(c)/RHO_W)**2*(K2_LIQ/K2_ICE)*1.0e18_rk*m6
      end if
    end function
  end subroutine column_diag

  !=================================================================
  ! mp_microwcp_driver   (<-> Thompson's mp_gt_driver)
  !   The WRF-facing entry point and intended drop-in replacement for
  !   mp_gt_driver. Interface notes for WRF compatibility:
  !
  !   * Arrays are declared default REAL (single, or whatever WRF is
  !     compiled with) and (ims:ime, kms:kme, jms:jme) -- k is the MIDDLE
  !     index. The per-column work is done in DOUBLE precision (the 1-D
  !     temporaries), exactly as Thompson does, so the validated physics
  !     core is unchanged.
  !   * Aerosol-aware path (nucleation) runs only if nwfa & nifa are
  !     PRESENT. If nc is absent it defaults to Nt_c/rho (1-moment-cloud
  !     fallback), matching Thompson's non-aerosol behaviour.
  !   * Precip: writes the per-step buckets *NCV and ACCUMULATES the
  !     running totals RAINNC/SNOWNC/GRAUPELNC (mm).
  !   * Optional radiation/radar diagnostics re_cloud/re_ice/re_snow
  !     (when has_reqc/i/s==1) and refl_10cm (when do_radar_ref==1).
  !   * Honours tile-vs-domain edges via MIN(ite,ide-1) like Thompson.
  !
  !   Call mp_microwcp_init() ONCE from WRF's serial mp_init before the
  !   threaded tile loop (the lazy guard here is only for standalone use).
  !=================================================================
  subroutine mp_microwcp_driver(qv,qc,qr,qi,qs,qg, ni,nr,nc, nwfa,nifa, &
       th,pii,p,w,dz, dt_in, itimestep,                                   &
       RAINNC,RAINNCV, SNOWNC,SNOWNCV, GRAUPELNC,GRAUPELNCV, SR,          &
       refl_10cm, do_radar_ref,                                          &
       re_cloud, re_ice, re_snow, has_reqc, has_reqi, has_reqs,          &
       ids,ide, jds,jde, kds,kde,                                        &
       ims,ime, jms,jme, kms,kme,                                        &
       its,ite, jts,jte, kts,kte)
    ! ---- dimensions ----
    integer, intent(in) :: ids,ide, jds,jde, kds,kde
    integer, intent(in) :: ims,ime, jms,jme, kms,kme
    integer, intent(in) :: its,ite, jts,jte, kts,kte
    ! ---- prognostic (INOUT), environment (IN) : default REAL = WRF kind
    real, dimension(ims:ime,kms:kme,jms:jme), intent(inout) :: &
                          qv,qc,qr,qi,qs,qg, ni,nr, th
    real, dimension(ims:ime,kms:kme,jms:jme), optional, intent(inout) :: &
                          nc, nwfa, nifa
    real, dimension(ims:ime,kms:kme,jms:jme), intent(in) :: pii,p,w,dz
    real, intent(in) :: dt_in
    integer, intent(in) :: itimestep
    ! ---- precip (INOUT) ----
    real, dimension(ims:ime,jms:jme), intent(inout) :: RAINNC, RAINNCV, SR
    real, dimension(ims:ime,jms:jme), optional, intent(inout) :: &
                          SNOWNC, SNOWNCV, GRAUPELNC, GRAUPELNCV
    ! ---- optional diagnostics ----
    real, dimension(ims:ime,kms:kme,jms:jme), optional, intent(inout) :: &
                          refl_10cm, re_cloud, re_ice, re_snow
    integer, optional, intent(in) :: do_radar_ref, has_reqc, has_reqi, has_reqs

    ! ---- locals : DOUBLE precision column work arrays ----
    integer  :: i,j,k, nz, i_end, j_end
    real(rk) :: dt
    real(rk) :: t1d(kts:kte), p1d(kts:kte), w1d(kts:kte), dz1d(kts:kte), rho(kts:kte)
    real(rk) :: qv1d(kts:kte),qc1d(kts:kte),qr1d(kts:kte),qi1d(kts:kte),qs1d(kts:kte),qg1d(kts:kte)
    real(rk) :: nc1d(kts:kte),nr1d(kts:kte),ni1d(kts:kte),nwfa1d(kts:kte),nifa1d(kts:kte)
    real(rk) :: re_qc(kts:kte), re_qi(kts:kte), re_qs(kts:kte), dbz(kts:kte)
    real(rk) :: pptrain,pptsnow,pptgraul,pptice
    logical  :: has_aero, has_nc, want_refl, want_re

    if (.not. initialized) call mp_microwcp_init()
    dt = real(dt_in, rk)
    nz = kte-kts+1
    has_aero = present(nwfa) .and. present(nifa)
    has_nc   = present(nc)
    want_refl = present(refl_10cm)
    if (present(do_radar_ref)) want_refl = want_refl .and. (do_radar_ref == 1)
    want_re   = present(re_cloud) .and. present(re_ice) .and. present(re_snow)

    ! tile loop, clipped to domain interior (Thompson convention)
    i_end = min(ite, ide-1)
    j_end = min(jte, jde-1)

    do j = jts, j_end
      do i = its, i_end
        ! ---- carve column, th->T, moist rho (single -> double) ----
        do k = kts, kte
          t1d(k)  = real(th(i,k,j),rk)*real(pii(i,k,j),rk)
          p1d(k)  = real(p(i,k,j),rk)
          w1d(k)  = real(w(i,k,j),rk)
          dz1d(k) = real(dz(i,k,j),rk)
          qv1d(k) = real(qv(i,k,j),rk)
          qc1d(k) = real(qc(i,k,j),rk); qr1d(k) = real(qr(i,k,j),rk)
          qi1d(k) = real(qi(i,k,j),rk); qs1d(k) = real(qs(i,k,j),rk); qg1d(k) = real(qg(i,k,j),rk)
          nr1d(k) = real(nr(i,k,j),rk); ni1d(k) = real(ni(i,k,j),rk)
          rho(k)  = 0.622_rk*p1d(k)/(RD*t1d(k)*(qv1d(k)+0.622_rk))
          if (has_nc) then
            nc1d(k) = real(nc(i,k,j),rk)
          else
            nc1d(k) = Nt_c/rho(k)                 ! fixed-Nc fallback
          end if
          if (has_aero) then
            ! clamp >= 0: WRF advection can leave tiny negatives in the aerosol
            ! scalars (roundoff dust relative to the ~1e8 field); negative aerosol
            ! is unphysical and must not enter nucleation.
            nwfa1d(k) = max(real(nwfa(i,k,j),rk),0.0_rk)
            nifa1d(k) = max(real(nifa(i,k,j),rk),0.0_rk)
          else
            nwfa1d(k) = 0.0_rk; nifa1d(k) = 0.0_rk
          end if
        end do

        ! ---- guard: abort only on genuinely broken (NaN/Inf) input from upstream;
        !      tiny negatives / roundoff are normal and handled by the physics ----
        do k = kts, kte
          if ( t1d(k) /= t1d(k) .or. qv1d(k) /= qv1d(k) .or. rho(k) /= rho(k) .or. &
               p1d(k) /= p1d(k) .or. dz1d(k) /= dz1d(k) ) then
            write(0,*) 'MICROWCP NaN input  i,k,j=',i,k,j,' step=',itimestep, &
                       '  t,p,rho,dz=',t1d(k),p1d(k),rho(k),dz1d(k)
            call wrf_error_fatal('MICROWCP received NaN input (upstream breakage)')
          end if
        end do

        ! ---- advance the column ----
        call mp_microwcp_column(nz, t1d,p1d,rho,w1d,dz1d, &
             qv1d,qc1d,qr1d,qi1d,qs1d,qg1d, nc1d,nr1d,ni1d, nwfa1d,nifa1d, &
             has_aero, dt, pptrain,pptsnow,pptgraul,pptice)

        ! ---- guard: abort only if the column produced NaN/Inf (genuine breakage).
        !      Physical soft-bound violations are clamped on write-back below, so a
        !      stray roundoff/edge value never aborts a (long) production run. ----
        do k = kts, kte
          if ( t1d(k) /= t1d(k) .or. abs(t1d(k)) > 1.0e4_rk .or. &
               qv1d(k) /= qv1d(k) .or. qc1d(k) /= qc1d(k) .or. qr1d(k) /= qr1d(k) .or. &
               qi1d(k) /= qi1d(k) .or. qs1d(k) /= qs1d(k) .or. qg1d(k) /= qg1d(k) ) then
            write(0,*) 'MICROWCP NaN/Inf output  i,k,j=',i,k,j,' step=',itimestep
            call wrf_error_fatal('MICROWCP produced NaN/Inf')
          end if
        end do

        ! ---- write back (double -> single), T->th; clamp to physical bounds ----
        do k = kts, kte
          qv(i,k,j) = max(real(qv1d(k)), 0.0)
          qc(i,k,j) = min(max(real(qc1d(k)),0.0), 0.1)
          qr(i,k,j) = min(max(real(qr1d(k)),0.0), 0.1)
          qi(i,k,j) = min(max(real(qi1d(k)),0.0), 0.1)
          qs(i,k,j) = min(max(real(qs1d(k)),0.0), 0.1)
          qg(i,k,j) = min(max(real(qg1d(k)),0.0), 0.1)
          nr(i,k,j) = min(max(real(nr1d(k)),0.0), 1.0e13)
          ni(i,k,j) = min(max(real(ni1d(k)),0.0), 1.0e13)
          th(i,k,j) = real(t1d(k)/real(pii(i,k,j),rk))
          if (has_nc)   nc(i,k,j)   = min(max(real(nc1d(k)),0.0), 1.0e13)
          if (has_aero) then
            nwfa(i,k,j) = max(real(nwfa1d(k)),0.0); nifa(i,k,j) = max(real(nifa1d(k)),0.0)
          end if
        end do

        ! ---- optional diagnostics ----
        if (want_re .or. want_refl) then
          call column_diag(nz, t1d,rho, qc1d,qr1d,qi1d,qs1d,qg1d, &
                           nc1d,nr1d,ni1d, re_qc,re_qi,re_qs, dbz)
          if (want_re) then
            do k = kts, kte
              re_cloud(i,k,j) = real(re_qc(k))
              re_ice(i,k,j)   = real(re_qi(k))
              re_snow(i,k,j)  = real(re_qs(k))
            end do
          end if
          if (want_refl) then
            do k = kts, kte
              refl_10cm(i,k,j) = real(dbz(k))
            end do
          end if
        end if

        ! ---- precip: per-step buckets + accumulated totals [mm] ----
        RAINNCV(i,j) = real(pptrain+pptsnow+pptgraul+pptice)
        RAINNC(i,j)  = RAINNC(i,j) + RAINNCV(i,j)
        if (present(SNOWNCV))    SNOWNCV(i,j)    = real(pptsnow+pptice)
        if (present(SNOWNC) .and. present(SNOWNCV)) SNOWNC(i,j) = SNOWNC(i,j) + SNOWNCV(i,j)
        if (present(GRAUPELNCV)) GRAUPELNCV(i,j) = real(pptgraul)
        if (present(GRAUPELNC) .and. present(GRAUPELNCV)) &
                                 GRAUPELNC(i,j)  = GRAUPELNC(i,j) + GRAUPELNCV(i,j)
        SR(i,j) = real((pptsnow+pptgraul+pptice)/(pptrain+pptsnow+pptgraul+pptice+1.0e-12_rk))
      end do
    end do
  end subroutine mp_microwcp_driver

end module module_mp_microwcp
