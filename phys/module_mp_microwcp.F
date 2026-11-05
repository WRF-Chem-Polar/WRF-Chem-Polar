!=======================================================================
! module_mp_microwcp
! Author: Anderson Da Silva, EERL-EPFL, anderson.dasilva@epfl.ch
! Last modified: 5 Nov 2026
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
!          order of sub-processes does not matter (the "sum all tendencies,
!          then apply once" design; Thompson does the same with its
!          prr_*/prv_* accumulators).
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
! PARAMETER CATEGORIES  —  where to tune, and what not to touch
!-----------------------------------------------------------------------
! (A) PHYSICAL CONSTANTS — fundamental; NEVER change these to alter behaviour.
!     PI, CP, RD, RV, P00, L_VAP, L_SUB, L_FUS, E0, T0, T_FREEZE,
!     T_HOM (-38 C, homogeneous-freezing limit), RHO_W/I/G, D_V, K_A, NU_AIR,
!     MU_AIR, RHO_NOT (fall-speed reference air density), RHO_STD (standard air).
!
! (B) PSD / PARTICLE DEFINITION — set the assumed size distributions and fall
!     speeds. Semi-fundamental: change only with a physical basis, not to tune.
!       mu, bm, am  (mass = am*D^bm) ;  av, bv (V = av*D^bv*sqrt(RHO_NOT/rho)) ;
!       cloud-droplet shape nu_c(N_c) (cloud_nu) ;
!       N0_graupel ;  snow_N0 (Houze et al. 1979).
!     Gb1, G1, Gq, GN, G2B, VG2 and the BR_*/CL_* tables are DERIVED from these
!     in mp_microwcp_init (not knobs).
!
! (C) PHYSICS TUNING KNOBS — these SET CLOUD AMOUNTS; this is where to tune the
!     scheme against observations. Each lives in the named subroutine:
!       nucleation : CCN activation Twomey P=1.36, k=0.462, w floor 0.1 m/s
!                                                                  -> droplet number
!                    DeMott DM_A..DM_D and its humidity gate DEMOTT_*_GATE ;
!                    cirrus coeff 0.1, cap 1 cm-3, RHi_crit fit (Koop)
!                                                                  -> ice number
!                    homogeneous-freezing ramp: onset T_HOM+3, width 2 K
!       warm       : accretion 67 / 1.15 (KK2000) ; rain self-collection K_RR
!                    and breakup size D_BR (Thompson) ; rain-evaporation
!                    ventilation FV*_EVAP. Autoconversion (Berry-Reinhardt)
!                    has no free knob.                             -> rain, LWP
!       cold       : ice->snow autoconv. D_I2S=150 um, TAU_AGG=120 s,
!                      sticking E_agg=exp(0.09*Tc)                 -> cloud ice / IWP
!                    Bigg freezing B_BIGG, A_BIGG (rain and droplets)
!                                                                  -> supercooled water
!                    snow capacitance C_SNOW_* ; melting C_MELT_*, FV*_MELT
!       collection : riming efficiencies: cloud ice E_CI=0.7, snow and graupel
!                    from the Stokes number ; rain-ice E_RI, snow-rain E_SR
!                                                                  -> riming rates
!                    Hallett-Mossop HM_SPL=3.5e8 per kg of rime, -8..-3 C
!                                                                  -> secondary ice
!
! (D) SECURITY / NUMERICAL THRESHOLDS — guards, NOT physics. Changing them must
!     not alter a well-resolved solution; do NOT use them to tune cloud amounts.
!       Q_MIN=1e-12, N_MIN=1     : floors that prevent divide-by-zero
!       QC/QR/QI/QS/QG_MIN       : below this a category is treated as absent
!       M_ICE_MIN=1e-12          : seed mass of a freshly nucleated crystal
!       D_I2S, TAU_AGG live in (C) but also gate ice->snow — tune there
!       DT_MAX=2 K               : max latent-heat T change per process per step
!                                  (energy-conserving spike limiter; sits ~500x
!                                  above normal heating, so it only ever bites a
!                                  pathological cell that would crash the host)
!       limit_scale              : the same conserving scale-back also stops the
!                                  summed sinks of one process group from
!                                  overdrawing a mixing ratio
!       write-back clamps (driver): q in [0,0.03], N in [0,1e13]; aerosols >= 0
!       NaN/Inf guards (driver)  : fatal — genuine upstream/scheme breakage only
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
  real(rk), parameter :: RHO_G = 400.0_rk
  ! Fall-speed laws are fitted near the ground; in thinner air a particle of
  ! the same size falls faster, by sqrt(RHO_NOT/rho) (Foote & du Toit 1969).
  ! RHO_NOT is Thompson's reference density (rho_not, module_mp_thompson.F).
  real(rk), parameter :: RHO_NOT = 101325.0_rk/(287.05_rk*298.0_rk)
  ! Density of standard air (0 C, 1013.25 hPa): DeMott et al. (2010) count
  ! aerosol and INP per standard volume (Thompson's rho_not0).
  real(rk), parameter :: RHO_STD = 101325.0_rk/(287.05_rk*273.15_rk)
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
  ! ice -> snow autoconversion (aggregation): size threshold + timescale
  real(rk), parameter :: D_I2S   = 150.0e-6_rk   ! mean-diam. ice->snow [m]
  real(rk), parameter :: TAU_AGG = 120.0_rk      ! aggregation timescale [s]

  ! ---- DeMott (2010) humidity gate (see `nucleation`) --------------
  ! DeMott measures how many INP are ACTIVE at a temperature, not that they
  ! nucleate in any parcel. Firing it in every subsaturated cell pinned ni at
  ! N_INP(T) everywhere and each crystal then drained the supercooled liquid by
  ! deposition (Wegener-Bergeron-Findeisen). Thompson gates it to
  ! ssati >= 0.25 .or. (ssatw > 0 .and. T < 253.15) [module_mp_thompson.F:2745].
  ! Same structure here, but the water-saturated branch opens at -10 C rather
  ! than Thompson's -20 C: Thompson makes NO DeMott ice between 0 and -20 C in
  ! water-saturated cloud, a large part of why it under-produces ice, and
  ! MICROWCP deliberately wants more mixed-phase ice than Thompson.
  real(rk), parameter :: DEMOTT_SI_GATE = 1.25_rk    ! RHi threshold [-]
  real(rk), parameter :: DEMOTT_T_GATE  = 263.15_rk  ! -10 C          [K]
  ! DeMott et al. (2010, PNAS, eq. 1):
  !   N_INP [std L-1] = DM_A (T0-T)^DM_B n_a^(DM_C (T0-T) + DM_D),
  ! with n_a the aerosol > 0.5 um [std cm-3] and T0 = 273.16 K.
  real(rk), parameter :: DM_A = 5.94e-5_rk, DM_B = 3.33_rk
  real(rk), parameter :: DM_C = 0.0264_rk,  DM_D = 0.0033_rk

  ! ---- Bigg (1953) immersion freezing of CLOUD DROPLETS ------------
  ! Rain already had Bigg freezing; cloud droplets did not, so once DeMott is
  ! gated there would be no primary ice source left in warm mixed-phase cloud.
  ! Thompson carries this as the tpi_qcfz lookup [module_mp_thompson.F:2733].
  ! The D^6 weighting makes it negligible above ~-25 C and dominant by -35 C,
  ! which is what lets observed supercooled liquid survive to -25 C.
  real(rk), parameter :: B_BIGG = 100.0_rk     ! freezing rate per volume [m-3 s-1]
  real(rk), parameter :: A_BIGG = 0.66_rk      ! supercooling exponent    [K-1]
  ! (the cloud-PSD moments <D^3>, <D^6> use the N_c-dependent droplet shape,
  !  see cloud_nu and the CL_* tables)

  ! ---- snow capacitance (see dep_rate use in `cold`) ---------------
  ! Snow is not a sphere: its capacitance runs from ~0.15*D for the warm
  ! dendritic/aggregate habit to 0.5*D (the sphere value) below about -30 C
  ! [module_mp_thompson.F:2809]. Using the sphere value at ALL temperatures
  ! overestimated the vapour sink competing with the supercooled liquid.
  real(rk), parameter :: C_SNOW_PLATE = 0.15_rk
  real(rk), parameter :: C_SNOW_SPH   = 0.50_rk
  real(rk), parameter :: T_C_SNOW_LO  = -1.5_rk    ! [C] plate habit here
  real(rk), parameter :: T_C_SNOW_HI  = -30.0_rk   ! [C] sphere value by here

  ! ---- Stokes-number riming efficiencies (see `collection`) --------
  ! A droplet is only collected if its inertia carries it across the
  ! streamlines going round the collector, so the efficiency is a function of
  ! the Stokes number and collapses to ZERO for the small droplets of a polar
  ! cloud. Snow: Wang & Ji (2000) [module_mp_thompson.F:5218]; graupel:
  ! [module_mp_thompson.F:2548]. The previous fixed values (0.5 snow, 1.0
  ! graupel) rimed 9-12 um polar droplets several times too fast.
  real(rk), parameter :: MU_AIR   = 1.718e-5_rk  ! dynamic viscosity [kg m-1 s-1]
  real(rk), parameter :: D_S_RIME = 300.0e-6_rk  ! snow below this cannot rime [m]
  real(rk), parameter :: D_C_RIME = 6.0e-6_rk    ! droplets below this not collected [m]
  ! fixed collection efficiencies of the other paths
  real(rk), parameter :: E_CI = 0.7_rk           ! cloud droplets by cloud ice
  real(rk), parameter :: E_RI = 1.0_rk           ! cloud ice by rain
  real(rk), parameter :: E_SR = 1.0_rk           ! rain by snow
  ! Hallett & Mossop (1974) rime splintering: splinters per kg of rime
  ! (350 per mg), weighted 0 at -8 C, 1 at -5 C, 0 at -3 C (Thompson's tf)
  real(rk), parameter :: HM_SPL = 3.5e8_rk

  ! ---- rain: evaporation, self-collection, breakup (see `warm`) ----
  ! Evaporation is diffusional growth integrated over the rain spectrum with
  ! ventilation f_v = FV1 + FV2 Sc^(1/3) Re^(1/2) (Thompson t1_qr_ev, t2_qr_ev).
  real(rk), parameter :: FV1_EVAP = 0.78_rk, FV2_EVAP = 0.308_rk
  ! Self-collection dNr/dt = -K_RR E_rr Nr (rho qr), E_rr = 1 - exp(2300 (D_mvd -
  ! D_BR)): coalescence for small drops, breakup (E_rr < 0) above D_BR
  ! (Seifert 1994; Verlinde & Cotton 1993; Thompson pnr_rcr).
  real(rk), parameter :: K_RR = 2.0_rk           ! [m3 kg-1 s-1]
  real(rk), parameter :: D_BR = 1.95e-3_rk       ! [m]

  ! ---- melting of snow and graupel (see `cold`) --------------------
  ! Capacitance as in Thompson: wet snow 0.15*D (C_sqrd), graupel 0.5*D
  ! (sphere, C_cube); ventilation f_v = 0.86 + 0.28 Sc^(1/3) Re^(1/2).
  real(rk), parameter :: C_MELT_S = 0.15_rk, C_MELT_G = 0.50_rk
  real(rk), parameter :: FV1_MELT = 0.86_rk, FV2_MELT = 0.28_rk

  ! ---- Berry & Reinhardt (1974) autoconversion (see `warm`) --------
  ! Replaced KK2000, a stratocumulus-DRIZZLE fit (calibrated qc <~ 1 g/kg) that
  ! could not drain a thick frontal/orographic cloud: at qc = 1 g/kg it is 10x
  ! slower than B-R and at 2 g/kg 30-49x slower, so condensate piled up as
  ! cloud instead of raining out. B-R is SLOWER than KK below qc ~ 0.3 g/kg and
  ! much faster above it, i.e. regime-selective - it does not thin weak clouds.
  real(rk), parameter :: BR_D0C = 1.0e-6_rk    ! smallest cloud droplet  [m]
  real(rk), parameter :: BR_D0R = 50.0e-6_rk   ! smallest rain drop      [m]
  integer,  parameter :: BR_NU_MIN = 2, BR_NU_MAX = 15
  real(rk) :: BR_G2_G1(BR_NU_MIN:BR_NU_MAX)    ! gamma(nu+4)/gamma(nu+1)
  real(rk) :: BR_G3_G2(BR_NU_MIN:BR_NU_MAX)    ! gamma(nu+7)/gamma(nu+4)

  ! ---- cloud-droplet shape: N_c-dependent, as in Thompson ----------
  ! Cloud droplets do NOT use the fixed mu(IC): like Thompson (nu_c =
  ! MIN(15, NINT(1000.E6/Nc)+2), module_mp_thompson.F) the gamma shape is an
  ! integer 2..15 set by the droplet number (cloud_nu), narrow in clean air.
  ! EVERY droplet process uses it: slope, mean diameter, fall speeds, Bigg
  ! freezing, effective radius, reflectivity and the B-R autoconversion. The
  ! former fixed mu = 3 was broader than observed maritime cloud (Martin et al.
  ! 1994, k ~ 0.8 <-> nu ~ 12): r_e ~13% too large and Bigg freezing ~2.4x too
  ! fast at given qc, Nc, and it disagreed with the autoconversion's nu_c.
  ! Tables over nu (BR_G2_G1 = gamma(nu+4)/gamma(nu+1) above is also used):
  real(rk) :: CL_GQ(BR_NU_MIN:BR_NU_MAX)   ! gamma(nu+bm+bv+1)/gamma(nu+bm+1): mass-weighted V
  real(rk) :: CL_GN(BR_NU_MIN:BR_NU_MAX)   ! gamma(nu+bv+1)/gamma(nu+1):       number-weighted V
  real(rk) :: CL_G7_G1(BR_NU_MIN:BR_NU_MAX) ! gamma(nu+7)/gamma(nu+1):          <D^6> lam^6
  ! Robustness: max latent-heat T change a single process may apply per step [K].
  ! If a process would heat/cool a cell by more than this, its WHOLE tendency
  ! (mass + heat) is scaled back by the same factor -> conserving; the residual
  ! super/subsaturation relaxes on the next step. Stops a phase-change spike
  ! (e.g. depositing all available vapour in a single step) from handing the
  ! host dynamics a buoyancy shock it cannot integrate (CFL blow-up). Ordinary
  ! clouds heat far below this, so the limiter only bites pathological cells.
  real(rk), parameter :: DT_MAX = 2.0_rk

  ! ---- PSD mass-diameter params: m = am*D^bm,  shape mu ----------
  !                          cloud      rain     ice      snow     graupel
  ! cloud entry NOT used: the droplet shape is N_c-dependent (cloud_nu)
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
  ! (the cloud entries are not used: droplets use the nu_c tables above)
  real(rk) :: Gb1(NCAT)   ! gamma(mu+bm+1)                    mass-moment norm
  real(rk) :: G1(NCAT)    ! gamma(mu+1)                       number-moment norm
  real(rk) :: Gq(NCAT)    ! gamma(mu+bm+bv+1)/gamma(mu+bm+1)  mass-weighted V
  real(rk) :: GN(NCAT)    ! gamma(mu+bv+1)/gamma(mu+1)        number-weighted V
  real(rk) :: G2B(NCAT)   ! gamma(mu+2bm+1)/gamma(mu+1)       <m^2> (reflectivity)
  real(rk) :: VG2(NCAT)   ! gamma(mu+2+(bv+1)/2)/gamma(mu+1)  ventilated moment
  logical  :: initialized = .false.

  ! ---- diagnostics constants (reflectivity / effective radii) ----
  real(rk), parameter :: Nt_c   = 100.0e6_rk   ! fixed cloud number when no aerosol [m-3]
  real(rk), parameter :: K2_ICE = 0.176_rk     ! ice dielectric factor |K|^2
  real(rk), parameter :: K2_LIQ = 0.93_rk      ! liquid dielectric factor |K|^2
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
  !   Only the Fortran intrinsic gamma is used, once, at initialisation.
  !     Gb1 = gamma(mu+bm+1)                    (mass-moment normalisation)
  !     G1  = gamma(mu+1)                       (number-moment normalisation)
  !     Gq  = gamma(mu+bm+bv+1)/gamma(mu+bm+1)  (mass-weighted bulk fall speed)
  !     GN  = gamma(mu+bv+1)/gamma(mu+1)        (number-weighted bulk fall speed)
  !     G2B = gamma(mu+2bm+1)/gamma(mu+1)       (mean squared mass: reflectivity)
  !     VG2 = gamma(mu+2+(bv+1)/2)/gamma(mu+1)  (ventilated first moment)
  !=================================================================
  subroutine mp_microwcp_init()
    integer :: c, n
    ! Berry-Reinhardt cloud-droplet shape nu_c is an INTEGER 2..15, so the two
    ! gamma ratios it needs are a 14-entry table rather than a runtime gamma
    ! call on every grid point (this mirrors Thompson's ccg/ocg tables).
    do n = BR_NU_MIN, BR_NU_MAX
      BR_G2_G1(n) = gamma(real(n,rk)+4.0_rk) / gamma(real(n,rk)+1.0_rk)
      BR_G3_G2(n) = gamma(real(n,rk)+7.0_rk) / gamma(real(n,rk)+4.0_rk)
      CL_GQ(n)    = gamma(real(n,rk)+bm(IC)+bv(IC)+1.0_rk) / gamma(real(n,rk)+bm(IC)+1.0_rk)
      CL_GN(n)    = gamma(real(n,rk)+bv(IC)+1.0_rk) / gamma(real(n,rk)+1.0_rk)
      CL_G7_G1(n) = gamma(real(n,rk)+7.0_rk) / gamma(real(n,rk)+1.0_rk)
    end do
    do c = 1, NCAT
      Gb1(c) = gamma(mu(c) + bm(c) + 1.0_rk)
      G1(c)  = gamma(mu(c) + 1.0_rk)
      ! mass-weighted bulk fall speed: weight V(D)=av*D^bv by mass m(D)=am*D^bm,
      ! so the gamma ratio must include the mass-diameter exponent bm.
      ! (Previously used gamma(mu+bv+1)/gamma(mu+1), the NUMBER-weighted factor,
      !  which made q sediment ~3x too slow for rain -> surface accumulation.)
      Gq(c)  = gamma(mu(c) + bm(c) + bv(c) + 1.0_rk) / gamma(mu(c) + bm(c) + 1.0_rk)
      GN(c)  = gamma(mu(c) + bv(c) + 1.0_rk) / gamma(mu(c) + 1.0_rk)
      G2B(c) = gamma(mu(c) + 2.0_rk*bm(c) + 1.0_rk) / gamma(mu(c) + 1.0_rk)
      VG2(c) = gamma(mu(c) + 2.0_rk + 0.5_rk*(bv(c) + 1.0_rk)) / gamma(mu(c) + 1.0_rk)
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
  ! Thermodynamics.  All elemental: act on scalars or whole arrays.
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

  ! Mass-weighted bulk fall speed [m s-1] for category c in air of density rho:
  !   V_q = av*Gq/lam^bv * sqrt(RHO_NOT/rho)
  elemental real(rk) function fall_speed_q(lam, c, rho)
    real(rk), intent(in) :: lam, rho
    integer,  intent(in) :: c
    fall_speed_q = av(c)*Gq(c) / lam**bv(c) * sqrt(RHO_NOT/rho)
  end function

  ! Number-weighted bulk fall speed [m s-1]:  V_N = av*GN/lam^bv * sqrt(RHO_NOT/rho)
  elemental real(rk) function fall_speed_N(lam, c, rho)
    real(rk), intent(in) :: lam, rho
    integer,  intent(in) :: c
    fall_speed_N = av(c)*GN(c) / lam**bv(c) * sqrt(RHO_NOT/rho)
  end function

  ! Ventilated first moment of the size distribution of category c [m-2]:
  !   int D f_v(D) n(D) dD,   f_v = f1 + f2 Sc^(1/3) Re^(1/2),  Re = V(D) D / nu,
  ! with V(D) = av D^bv sqrt(RHO_NOT/rho). For n = N0 D^mu exp(-lam D),
  ! N0 = Nb lam^(mu+1)/gamma(mu+1), and h = (bv+1)/2:
  !   = Nb [ f1 (mu+1)/lam + f2 Sc^(1/3) (av sqrt(RHO_NOT/rho)/nu)^(1/2) VG2 / lam^(1+h) ]
  ! This is how Thompson integrates ventilation; used for rain evaporation and
  ! melting. (f1 alone, without ventilation, gives the capacitance integral
  ! 2*int (D/2) n dD.)
  elemental real(rk) function vent_moment(Nb, lam, c, rho, f1, f2)
    real(rk), intent(in) :: Nb, lam, rho, f1, f2
    integer,  intent(in) :: c
    real(rk) :: h
    h = 0.5_rk*(bv(c) + 1.0_rk)
    vent_moment = Nb*( f1*(mu(c) + 1.0_rk)/lam                              &
                     + f2*(NU_AIR/D_V)**(1.0_rk/3.0_rk)                       &
                       *sqrt(av(c)*sqrt(RHO_NOT/rho)/NU_AIR)*VG2(c)/lam**(1.0_rk + h) )
  end function

  ! Cloud-droplet gamma shape nu_c (integer 2..15) from the bulk droplet
  ! number [m-3], exactly as Thompson: NINT(1000e6/Nc)+2, capped at 15, with Nc
  ! floored at 1 m-3. The single definition every droplet process uses.
  elemental integer function cloud_nu(nc_bulk)
    real(rk), intent(in) :: nc_bulk
    cloud_nu = min(BR_NU_MAX, max(BR_NU_MIN, nint(1000.0e6_rk/max(nc_bulk,1.0_rk)) + 2))
  end function

  ! Cloud-droplet bulk fall speeds [m s-1] for shape nu (mass- and number-
  ! weighted), with the same air-density factor as the other categories.
  elemental real(rk) function cloud_vq(lam, nu, rho)
    real(rk), intent(in) :: lam, rho
    integer,  intent(in) :: nu
    cloud_vq = av(IC)*CL_GQ(nu) / lam**bv(IC) * sqrt(RHO_NOT/rho)
  end function

  elemental real(rk) function cloud_vn(lam, nu, rho)
    real(rk), intent(in) :: lam, rho
    integer,  intent(in) :: nu
    cloud_vn = av(IC)*CL_GN(nu) / lam**bv(IC) * sqrt(RHO_NOT/rho)
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
    integer  :: k, nu
    real(rk) :: ratio

    do k = 1, nz
      ! ---- 2-moment: cloud, rain, ice  (lam from N_bulk, q) ----
      ! cloud: shape nu_c from the droplet number (cloud_nu), so
      ! lam = [am*N*gamma(nu+4)/(rho*q*gamma(nu+1))]**(1/3)
      nu = cloud_nu(nc(k)*rho(k))
      Nbulk(k,IC) = nc(k)*rho(k)
      lam(k,IC) = (am(IC)*max(nc(k)*rho(k), N_MIN)*BR_G2_G1(nu) &
                   /(rho(k)*max(qc(k), Q_MIN)))**(1.0_rk/bm(IC))
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
      Dm(k,IC) = (real(nu,rk)+bm(IC)+1.0_rk)/lam(k,IC)
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
  !=================================================================
  subroutine upstream_flux(nz, q, V, rho, dz, dt, dq, flx)
    integer,  intent(in)  :: nz
    real(rk), intent(in)  :: q(nz), V(nz), rho(nz), dz(nz), dt
    real(rk), intent(out) :: dq(nz)
    ! optional: downward flux out of the bottom of each level, averaged over
    ! the sub-steps [kg m-2 s-1] (Thompson's precr/preci/precs/precg)
    real(rk), intent(out), optional :: flx(nz)
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
    if (present(flx)) flx = 0.0_rk
    do s = 1, nsub
      do k = 1, nz
        flux(k) = rho(k)*V(k)*qw(k)
      end do
      if (present(flx)) flx = flx + flux/real(nsub, rk)
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
  ! limit_scale
  !   Conserving limiter of one process group (nucleation, warm, cold,
  !   collection) in one cell: the largest s <= 1 such that, once all the
  !   group's tendencies are multiplied by s,
  !     (i)  its latent heating changes T by at most DT_MAX this step, and
  !     (ii) no mixing ratio x0 is driven below zero by the SUM of its
  !          sinks (each sub-process is capped at x0/dt on its own, but two
  !          sinks of the same category can together exceed x0).
  !   Scaling every tendency of the group by the same s keeps water and
  !   energy exact (Thompson rescales the sinks of an overdrawn category
  !   to the same end); what is left relaxes on the next step.
  !=================================================================
  pure real(rk) function limit_scale(aT, dt, x0, ax)
    real(rk), intent(in) :: aT, dt, x0(:), ax(:)
    integer :: i
    limit_scale = 1.0_rk
    if (abs(aT)*dt > DT_MAX) limit_scale = DT_MAX/(abs(aT)*dt)
    do i = 1, size(x0)
      if (x0(i) + ax(i)*dt < 0.0_rk) &
        limit_scale = min(limit_scale, max(x0(i),0.0_rk)/(-ax(i)*dt))
    end do
  end function limit_scale

  !=================================================================
  ! mp_microwcp_column   (<-> Thompson's mp_thompson)
  !   Advance ONE atmospheric column by one time step dt. This is the
  !   heart of the scheme: it runs the operator-split sequence
  !     1 nucleation -> 2 warm -> 3 cold -> 4 collection
  !     -> 5 sedimentation -> 6 melting of cloud ice above 0 C
  !     -> 7 saturation adjustment -> 8 number/mass consistency
  !   re-diagnosing the PSD (diagnose) before each group so particle sizes
  !   stay consistent with the updated state. Prognostic fields are 1-D
  !   (length nz, index 1 = lowest level) and updated in place. Surface
  !   precip is returned per category as depth [mm].
  !   Temperature is carried as FULL T here; the th<->T conversion lives
  !   in the driver (matching Thompson).
  !=================================================================
  subroutine mp_microwcp_column(nz, T, P, rho, w, dz, &
       qv,qc,qr,qi,qs,qg, nc,nr,ni, nwfa,nifa, has_aero, dt, &
       pptrain, pptsnow, pptgraul, pptice, &
       c2prec, rainprod, evapprod, precr, preci, precs, precg)
    integer,  intent(in)    :: nz
    real(rk), intent(inout) :: T(nz)
    real(rk), intent(in)    :: P(nz), rho(nz), w(nz), dz(nz)
    real(rk), intent(inout) :: qv(nz),qc(nz),qr(nz),qi(nz),qs(nz),qg(nz)
    real(rk), intent(inout) :: nc(nz),nr(nz),ni(nz),nwfa(nz),nifa(nz)
    logical,  intent(in)    :: has_aero
    real(rk), intent(in)    :: dt
    real(rk), intent(out)   :: pptrain, pptsnow, pptgraul, pptice
    ! optional WRF-Chem wet-scavenging diagnostics: the same processes as
    ! Thompson's (module_mp_thompson.F, WRF_CHEM blocks), in the units the
    ! WRF-Chem routines document (kg kg-1 s-1, as Morrison provides them).
    ! Thompson's own rainprod/c2prec are per unit volume (kg m-3 s-1), so its
    ! rainprod and qlsink are rho times these; and its evapprod also counts
    ! snow/graupel sublimation (used only for resuspension, off in GOCART).
    !   c2prec    cloud water -> precipitation (autoconversion, accretion,
    !             riming)                                       [kg kg-1 s-1]
    !   rainprod  c2prec + cloud ice -> precipitation            [kg kg-1 s-1]
    !   evapprod  rain evaporation                               [kg kg-1 s-1]
    !   precr/i/s/g  sedimentation flux out of each level         [kg m-2 s-1]
    real(rk), intent(out), optional :: c2prec(nz), rainprod(nz), evapprod(nz), &
                                       precr(nz), preci(nz), precs(nz), precg(nz)

    real(rk) :: lam(nz,NCAT), Dm(nz,NCAT), Nb(nz,NCAT), N0s(nz), N0g(nz)
    real(rk) :: dqr_sed(nz), dqi_sed(nz), dqs_sed(nz), dqg_sed(nz)
    real(rk) :: qc_b(nz), qr_b(nz), qi_b(nz), agg(nz)
    logical  :: want_scav

    if (.not. initialized) call mp_microwcp_init()

    ! ---- 1. nucleation ----
    if (has_aero) call nucleation(nz, T,P,rho,w, qv,qc,qi, nc,ni, nwfa,nifa, dt)

    want_scav = present(c2prec) .and. present(rainprod) .and. present(evapprod) &
                .and. present(precr) .and. present(preci) .and. present(precs) &
                .and. present(precg)

    ! ---- 2. warm phase ----
    !      warm changes qc only by autoconversion + accretion (both -> rain),
    !      and qr by those plus rain evaporation, so the scavenging rates
    !      follow from differencing the state around the call.
    call diagnose(nz,T,rho, qc,qr,qi,qs,qg, nc,nr,ni, lam,Dm,Nb,N0s,N0g)
    if (want_scav) then
      qc_b = qc; qr_b = qr
    end if
    call warm(nz, T,P,rho, qv,qc,qr, nc,nr, lam,Nb, dt)
    if (want_scav) then
      c2prec   = max((qc_b-qc)/dt, 0.0_rk)
      evapprod = max((qc_b-qc)/dt - (qr-qr_b)/dt, 0.0_rk)
    end if

    ! ---- 3. cold phase ----
    !      (cloud-droplet Bigg freezing makes cloud ICE, not precipitation, so
    !      only the ice -> snow aggregation is counted, in rainprod)
    call diagnose(nz,T,rho, qc,qr,qi,qs,qg, nc,nr,ni, lam,Dm,Nb,N0s,N0g)
    call cold(nz, T,P,rho, qv,qc,qr,qi,qs,qg, nc,nr,ni, lam,Dm,Nb, dt, agg)

    ! ---- 4. collection (+ Hallett-Mossop splintering) ----
    !      qc falls only by riming (-> snow/graupel), qi only by rain-ice
    !      collection (-> graupel); rime splinters add a little cloud ice back,
    !      which the max(.,0) below simply nets out
    call diagnose(nz,T,rho, qc,qr,qi,qs,qg, nc,nr,ni, lam,Dm,Nb,N0s,N0g)
    if (want_scav) then
      qc_b = qc; qi_b = qi
    end if
    call collection(nz, T,rho, qc,qr,qi,qs,qg, nc,nr,ni, lam,Dm,Nb, dt)
    if (want_scav) then
      c2prec   = c2prec + max((qc_b-qc)/dt, 0.0_rk)
      rainprod = c2prec + agg + max((qi_b-qi)/dt, 0.0_rk)
    end if

    ! ---- 5. sedimentation ----
    call diagnose(nz,T,rho, qc,qr,qi,qs,qg, nc,nr,ni, lam,Dm,Nb,N0s,N0g)
    call sedimentation(nz, rho,dz,dt, qc,qr,qi,qs,qg, nc,nr,ni, lam, &
                       dqr_sed,dqi_sed,dqs_sed,dqg_sed, precr,preci,precs,precg)

    ! ---- 6. cloud ice that sedimented or was advected above 0 C melts ----
    call melt_cloud_ice(nz, T, qc, qi, nc, ni)

    ! ---- 7. saturation adjustment ----
    call saturation_adjustment(nz, T,P, qv,qc,qi,qs,qg, nc)

    ! ---- 8. number<->mass consistency for frozen 2-moment categories ----
    !      zero the number where mass is negligible so advected number cannot
    !      strand in condensate-free cells (cloud droplets handled in satadj).
    where (qi < QI_MIN) ni = 0.0_rk
    where (qr < QR_MIN) nr = 0.0_rk

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
  !     - CCN activation (Twomey): updraft + water-friendly aerosol nwfa
  !       -> new cloud droplets (nc,qc), above -35 C.
  !     - Heterogeneous ice nucleation (DeMott 2010), humidity-gated:
  !       ice-friendly aerosol nifa, between -38 and 0 C -> new crystals
  !       (ni,qi); nifa is consumed.
  !     - Homogeneous freezing of droplets: a smooth ramp from -35 C,
  !       complete by about -40 C (qc,nc -> qi,ni).
  !     - Cirrus (homogeneous freezing of haze aerosol) below -35 C.
  !   Seed condensate is drawn from vapour (qv) so water is conserved.
  !   Only runs when aerosols are present (has_aero in the column routine;
  !   always the case in WRF, where nwfa/nifa are carried for MICROWCP).
  !-----------------------------------------------------------------
  subroutine nucleation(nz, T,P,rho,w, qv,qc,qi, nc,ni, nwfa,nifa, dt)
    integer,  intent(in)    :: nz
    real(rk), intent(in)    :: P(nz),rho(nz),w(nz),nwfa(nz),dt
    ! nifa is INOUT: an activated INP is locked inside an ice crystal and cannot
    ! nucleate again, so the reservoir must be consumed (Thompson :3100).
    ! Without depletion the deficit-based form below makes NIFA a floor that
    ! can never be exhausted.
    real(rk), intent(inout) :: T(nz),qv(nz),qc(nz),qi(nz),nc(nz),ni(nz),nifa(nz)
    ! All sub-process rates are evaluated from the FROZEN input state and
    ! accumulated, then applied once (matches stepper.step's tendency sum).
    integer  :: k
    real(rk) :: T0k,qv0,qc0,qi0,nc0,ni0
    real(rk) :: aqv,aqc,aqi,anc,ani,aT
    real(rk) :: qsat,Nact,dNc,m_drop,dqc,C_cm3,weff,f_hf
    real(rk) :: na_cm3,Tsc,Ninp_L,Ninp,dNi,dqi
    real(rk) :: Si,RHi_crit,N_hom,scal
    real(rk) :: nifa0,anifa,RHi,ssatw
    logical  :: fires
    m_drop = (4.0_rk/3.0_rk)*PI*(1.0e-6_rk)**3*1000.0_rk
    do k = 1, nz
      T0k=T(k); qv0=qv(k); qc0=qc(k); qi0=qi(k); nc0=nc(k); ni0=ni(k)
      nifa0=nifa(k)
      aqv=0.0_rk; aqc=0.0_rk; aqi=0.0_rk; anc=0.0_rk; ani=0.0_rk; aT=0.0_rk
      anifa=0.0_rk
      ! --- CCN activation (Twomey 1959, updraft-based) ---
      ! Activated number is set by the updraft peak supersaturation (scales with
      ! w), NOT the post-saturation-adjustment S (which is ~0, so Nc collapses ->
      ! giant droplets -> runaway autoconversion -> no liquid cloud). For a Twomey
      ! spectrum N=C*S^k:  N_act = P*C^(2/(k+2))*w^(3k/(2(k+2))), P=1.36, k=0.462.
      ! Do NOT activate in the homogeneous-freezing regime (T < T_HOM+3 ~ -35C):
      ! liquid cloud cannot persist there, so a droplet activated here is frozen
      ! by homogeneous freezing the same/next step, and because activation is
      ! deficit-based it immediately re-activates -> a vapour->tiny-droplet->
      ! tiny-ice-crystal NUMBER pump that drives ni to ~1e10/kg with crystals
      ! far below M_ICE_MIN (an unphysical ice-number runaway). Gating to T
      ! above the homogeneous onset breaks the conveyor; warm/mixed-phase
      ! activation (supercooled liquid) is unaffected.
      qsat = qsat_liquid(T0k,P(k))
      if ((qv0 >= 0.99_rk*qsat .or. qc0 > Q_MIN) .and. T0k > T_HOM+3.0_rk) then
        C_cm3 = max(nwfa(k)*rho(k)*1.0e-6_rk, 1.0e-3_rk)        ! CCN [cm-3]
        weff  = max(w(k), 0.1_rk)                              ! m/s; turbulent floor
        ! exponents written as the exact Twomey expressions (k = 0.462), not as
        ! rounded literals: 0.8127/0.2815 differed from 2/(k+2)=0.81235 and
        ! 3k/(2(k+2))=0.28148 enough to shift N_act by ~0.1-0.2% vs the prototype
        Nact  = 1.36_rk * C_cm3**(2.0_rk/(0.462_rk+2.0_rk)) &
                        * weff**(3.0_rk*0.462_rk/(2.0_rk*(0.462_rk+2.0_rk)))  ! [cm-3]
        Nact  = min(Nact*1.0e6_rk/max(rho(k),1.0e-3_rk), nwfa(k))  ! -> [kg-1]
        dNc   = max(Nact-nc0,0.0_rk)/dt
        dqc   = dNc*m_drop
        anc=anc+dNc; aqc=aqc+dqc; aqv=aqv-dqc; aT=aT+L_VAP/CP*dqc
      end if
      ! --- heterogeneous ice nucleation (DeMott 2010), humidity-gated ---
      !     Gate (see DEMOTT_*_GATE): nucleate only where the air is
      !     ice-supersaturated past RHi=1.25, or water-saturated and colder than
      !     -10 C. Ungated, this fired in every cell of the mixed-phase zone and
      !     pinned ni at N_INP(T) permanently, glaciating supercooled cloud.
      Tsc = max(T0+0.01_rk - T0k, 0.0_rk)            ! 273.16 - T
      if (T0k > T_HOM .and. T0k < T_FREEZE .and. Tsc > 0.0_rk) then
        RHi   = qv0/max(qsat_ice(T0k,P(k)),Q_MIN)
        ssatw = qv0/max(qsat_liquid(T0k,P(k)),Q_MIN) - 1.0_rk
        fires = (RHi >= DEMOTT_SI_GATE) .or. &
                (ssatw > 0.0_rk .and. T0k < DEMOTT_T_GATE)
        if (fires) then
          ! DeMott's aerosol and INP numbers are per STANDARD volume (0 C,
          ! 1013.25 hPa), as in Thompson: a kg of air fills 1/RHO_STD std m3
          na_cm3 = max(nifa0*RHO_STD*1.0e-6_rk, 1.0e-6_rk)       ! [std cm-3]
          Ninp_L = DM_A*(Tsc**DM_B)*(na_cm3**(DM_C*Tsc+DM_D))     ! [std L-1]
          Ninp   = Ninp_L*1.0e3_rk/RHO_STD                         ! [kg-1]
          dNi    = max(Ninp-ni0,0.0_rk)/dt
          dqi    = dNi*M_ICE_MIN
          ani=ani+dNi; aqi=aqi+dqi; aqv=aqv-dqi; aT=aT+L_SUB/CP*dqi
          anifa=anifa-dNi                         ! INP consumed by nucleation
        end if
      end if
      ! --- homogeneous freezing (qc->qi, nc->ni): smooth temperature ramp,
      !     not an instant cutoff. Onset ~-35 C, effectively complete ~-40 C. ---
      f_hf = 1.0_rk - exp(-max((T_HOM+3.0_rk-T0k)/2.0_rk, 0.0_rk))
      if (f_hf > 0.0_rk) then
        dqi = f_hf*qc0/dt; dNi = f_hf*nc0/dt
        aqi=aqi+dqi; ani=ani+dNi; aqc=aqc-dqi; anc=anc-dNi
        aT=aT+L_FUS/CP*dqi
      end if
      ! --- homogeneous aerosol (cirrus) freezing below T_HOM+3 (~-35 C) ---
      !     Solution/haze aerosol freezes once the ice supersaturation exceeds
      !     the Koop (2000) threshold (Ren & MacKenzie 2005 fit). Crystal number
      !     is updraft-scaled (Karcher & Lohmann 2002), capped, deficit-based and
      !     bounded by aerosol -> makes cirrus ice where droplet freezing cannot,
      !     WITHOUT a runaway. It starts where CCN activation stops (as Thompson's
      !     iceKoop, T < 238 K), so in the -35..-38 C band supersaturated
      !     droplet-free air freezes as haze. Aerosol not depleted (diagnostic
      !     field, as for CCN activation and DeMott).
      if (T0k < T_HOM+3.0_rk) then
        Si = qv0/max(qsat_ice(T0k,P(k)),Q_MIN)
        RHi_crit = min(max(2.418_rk - T0k/245.68_rk, 1.4_rk), 1.7_rk)
        if (Si > RHi_crit) then
          N_hom = min(0.1_rk*max(w(k),0.05_rk)**1.5_rk, 1.0_rk)       ! [cm-3]
          N_hom = min(N_hom*1.0e6_rk/max(rho(k),1.0e-3_rk), nwfa(k))  ! [kg-1]
          dNi = max(N_hom-ni0,0.0_rk)/dt
          dqi = dNi*M_ICE_MIN
          ani=ani+dNi; aqi=aqi+dqi; aqv=aqv-dqi; aT=aT+L_SUB/CP*dqi
        end if
      end if
      ! --- conserving limiter (see limit_scale): latent heating within
      !     DT_MAX, and no mixing ratio overdrawn by its summed sinks ---
      scal = limit_scale(aT, dt, (/qv0,qc0,qi0/), (/aqv,aqc,aqi/))
      if (scal < 1.0_rk) then
        aqv=aqv*scal; aqc=aqc*scal; aqi=aqi*scal
        anc=anc*scal; ani=ani*scal; aT=aT*scal
        anifa=anifa*scal
      end if
      ! --- apply once ---
      qc(k)=max(qc0+aqc*dt,0.0_rk); qv(k)=max(qv0+aqv*dt,0.0_rk)
      qi(k)=max(qi0+aqi*dt,0.0_rk)
      nc(k)=max(nc0+anc*dt,0.0_rk); ni(k)=max(ni0+ani*dt,0.0_rk)
      nifa(k)=max(nifa0+anifa*dt,0.0_rk)
      T(k) =T0k+aT*dt
    end do
  end subroutine nucleation

  !-----------------------------------------------------------------
  !-----------------------------------------------------------------
  ! warm
  !   Liquid-only processes that make and shape rain:
  !     - autoconversion (Berry & Reinhardt 1974, as Thompson): cloud
  !       droplets coalesce into rain.
  !     - accretion (Khairoutdinov & Kogan 2000): rain sweeps up droplets.
  !     - rain evaporation: diffusional, ventilated, in subsaturated air.
  !     - self-collection / breakup (Thompson): adjust the rain DROP NUMBER
  !       nr only (mass conserved).
  !-----------------------------------------------------------------
  subroutine warm(nz, T,P,rho, qv,qc,qr, nc,nr, lam,Nb, dt)
    integer,  intent(in)    :: nz
    real(rk), intent(in)    :: P(nz),rho(nz),lam(nz,NCAT),Nb(nz,NCAT),dt
    real(rk), intent(inout) :: T(nz),qv(nz),qc(nz),qr(nz),nc(nz),nr(nz)
    integer  :: k
    real(rk) :: qv0,qc0,qr0,nc0,nr0
    real(rk) :: aqv,aqc,aqr,anc,anr,aT
    real(rk) :: dqr,dNr,dNc
    real(rk) :: qsat,sw,Aw,Bw,dq_def,dqr_e,scal,mvd_r,E_rr
    real(rk) :: am_r,rc,rc_s,nc_bulk,lamc,xDc,Dc_g,Dc_b,zeta1,zeta,taud,tau
    real(rk) :: mvd_c,g2g1,g3g2
    integer  :: nu_c
    am_r = (PI/6.0_rk)*RHO_W
    do k = 1, nz
      qv0=qv(k); qc0=qc(k); qr0=qr(k); nc0=nc(k); nr0=nr(k)
      aqv=0.0_rk; aqc=0.0_rk; aqr=0.0_rk; anc=0.0_rk; anr=0.0_rk; aT=0.0_rk
      ! --- autoconversion: Berry & Reinhardt (1974) ---
      !     The rate scales with the CHARACTERISTIC droplet diameter, which
      !     grows with qc, so it switches on hard once the cloud thickens. The
      !     previous KK2000 form (1350 qc^2.47 Nc^-1.79) is a stratocumulus
      !     drizzle fit: 10x slower at qc = 1 g/kg and 30-49x slower at 2 g/kg,
      !     so thick maritime cloud could not drain and LWP piled up.
      rc = rho(k)*qc0                              ! cloud water per m3 [kg m-3]
      if (qc0 > QC_MIN .and. rc > 0.01e-3_rk) then
        rc_s    = max(rc, 1.0e-12_rk)
        nc_bulk = max(nc0*rho(k), 1.0_rk)
        ! droplet spectral shape: integer 2..15, gamma ratios from the table
        nu_c = cloud_nu(nc_bulk)
        g2g1 = BR_G2_G1(nu_c)
        g3g2 = BR_G3_G2(nu_c)
        lamc = (nc_bulk*am_r*g2g1/rc_s)**(1.0_rk/3.0_rk)
        ! characteristic diameters, in micrometres (as in Thompson)
        xDc  = max(BR_D0C*1.0e6_rk, (rc_s/(am_r*nc_bulk))**(1.0_rk/3.0_rk)*1.0e6_rk)
        Dc_g = g3g2**(1.0_rk/3.0_rk)/lamc*1.0e6_rk
        Dc_b = max(xDc**3*Dc_g**3 - xDc**6, 0.0_rk)**(1.0_rk/6.0_rk)
        zeta1 = max(6.25e-6_rk*xDc*Dc_b**3 - 0.4_rk, 0.0_rk)
        zeta  = 0.027_rk*rc_s*zeta1
        taud  = max(0.5_rk*Dc_b - 7.5_rk, 0.0_rk) + 1.0e-12_rk
        tau   = 3.72_rk/(rc_s*taud)
        dqr   = min(zeta/tau/rho(k), qc0/dt)       ! [kg m-3 s-1] -> [kg kg-1 s-1]
        if (dqr > 0.0_rk) then
          dNr   = dqr/(am_r*real(nu_c,rk)*10.0_rk*BR_D0R**3)
          mvd_c = min(max((3.0_rk+real(nu_c,rk)+0.672_rk)/lamc, BR_D0C), BR_D0R)
          dNc   = min(dqr/(am_r*mvd_c**3), nc0/dt)
          aqr=aqr+dqr; aqc=aqc-dqr; anr=anr+dNr; anc=anc-dNc
        end if
      end if
      ! --- accretion ---
      if (qc0 > QC_MIN .and. qr0 > QR_MIN) then
        dqr = min(67.0_rk*(qc0*qr0)**1.15_rk, qc0/dt)
        ! Droplet NUMBER is collected with the mass (Thompson pnc_rcw, :2325).
        ! Removing mass alone left the survivors ever smaller as the cloud
        ! drained, holding autoconversion down and keeping LWP high.
        dNc = min(nc0*dqr/max(qc0,Q_MIN), nc0/dt)
        aqr=aqr+dqr; aqc=aqc-dqr; anc=anc-dNc
      end if
      ! --- rain evaporation: diffusional growth over the rain spectrum ---
      !     dqr/dt = 2 pi (S_w - 1) int(D f_v n dD) / (rho (A_w + B_w)),
      !     A_w = L_v^2/(K_a R_v T^2), B_w = R_v T/(D_v e_sw) (Pruppacher &
      !     Klett 1997, eq. 13-28), the ventilation f_v integrated over the
      !     spectrum as in Thompson (Srivastava & Coen 1992, first order in
      !     S_w - 1). At most the rain present, and at most the vapour deficit
      !     that brings the air back to water saturation (latent cooling
      !     included). Drop number falls in proportion to mass.
      if (qr0 > QR_MIN) then
        qsat = qsat_liquid(T(k),P(k))
        sw   = qv0/max(qsat,Q_MIN)
        if (sw < 1.0_rk) then
          Aw     = L_VAP**2/(K_A*RV*T(k)**2)
          Bw     = RV*T(k)/(D_V*es_liquid(T(k)))
          dqr_e  = 2.0_rk*PI*(1.0_rk-sw)                                   &
                   *vent_moment(max(Nb(k,IR),N_MIN), lam(k,IR), IR, rho(k), &
                                FV1_EVAP, FV2_EVAP)/(rho(k)*(Aw+Bw))
          dq_def = (qsat-qv0)/(1.0_rk + L_VAP**2*qsat/(RV*CP*T(k)**2))
          dqr_e  = min(dqr_e, qr0/dt, dq_def/dt)
          aqr=aqr-dqr_e; aqv=aqv+dqr_e
          anr=anr - nr0*dqr_e/max(qr0,Q_MIN); aT=aT - L_VAP/CP*dqr_e
        end if
      end if
      ! --- rain self-collection and breakup (number only) ---
      !     dNr/dt = -K_RR E_rr Nr rho qr,  E_rr = 1 - exp(2300 (D_mvd - D_BR)),
      !     with the median-volume diameter D_mvd = (3.672 + mu_r)/lam_r
      !     (Thompson pnr_rcr; Seifert 1994, Verlinde & Cotton 1993). For
      !     drops above D_BR = 1.95 mm, E_rr < 0 turns it into breakup. As in
      !     Thompson, only for D_mvd > D0r = 50 um; capped at Nr/dt either way.
      if (qr0 > QR_MIN .and. nr0 > 0.0_rk) then
        mvd_r = (3.672_rk + mu(IR))/lam(k,IR)
        if (mvd_r > BR_D0R) then
          E_rr = 1.0_rk - exp(min(2300.0_rk*(mvd_r - D_BR), 50.0_rk))
          dNr  = K_RR*E_rr*nr0*rho(k)*qr0
          anr  = anr - min(max(dNr, -nr0/dt), nr0/dt)
        end if
      end if
      ! --- conserving limiter (see limit_scale) ---
      scal = limit_scale(aT, dt, (/qv0,qc0,qr0/), (/aqv,aqc,aqr/))
      if (scal < 1.0_rk) then
        aqv=aqv*scal; aqc=aqc*scal; aqr=aqr*scal
        anc=anc*scal; anr=anr*scal; aT=aT*scal
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
  !       diffusion-limited with a ventilation correction (dep_rate), the
  !       three together capped at the vapour excess over ice saturation.
  !     - Bigg immersion freezing: supercooled rain -> graupel, cloud
  !       droplets -> cloud ice.
  !     - ice -> snow autoconversion of crystals grown past D_I2S.
  !     - melting: snow/graupel -> rain above 0 C (cools; adds rain drops).
  !   (Hallett-Mossop splintering lives in `collection`, where the riming
  !   rate it scales with is computed.)
  !-----------------------------------------------------------------
  subroutine cold(nz, T,P,rho, qv,qc,qr,qi,qs,qg, nc,nr,ni, lam,Dm,Nb, dt, agg)
    integer,  intent(in)    :: nz
    real(rk), intent(in)    :: P(nz),rho(nz),lam(nz,NCAT),Dm(nz,NCAT),Nb(nz,NCAT),dt
    real(rk), intent(inout) :: T(nz),qv(nz),qc(nz),qr(nz),qi(nz),qs(nz),qg(nz)
    real(rk), intent(inout) :: nc(nz),nr(nz),ni(nz)
    ! ice -> snow aggregation rate actually applied [kg kg-1 s-1] (wet scavenging)
    real(rk), intent(out)   :: agg(nz)
    integer  :: k
    real(rk) :: T0k,qv0,qc0,qr0,qi0,qs0,qg0,nc0,nr0,ni0
    real(rk) :: aqv,aqc,aqr,aqi,aqs,aqg,anc,anr,ani,aT
    real(rk) :: dqi,dqs,dqg,dqr_f,dqs_a
    real(rk) :: heat,qsi,dq_ex,dsum,rate_s,rate_g
    real(rk) :: C_bigg, dTf, E_agg, scal
    real(rk) :: C_snow, fac_b, dqc_f, dNc_f, lamc
    integer  :: nu_b
    ! Bigg immersion freezing of rain -> graupel, bulk-integrated over the
    ! rain PSD: dqr/dt = C_bigg*(exp(A_BIGG*dT)-1)*nr*Dm^6. The Dm^6 (6th
    ! moment) makes small drizzle drops freeze slowly so supercooled rain
    ! persists. C_bigg = rho_w*(pi/6)^2 * B_BIGG * Gamma(7)/4^6 (720/4096).
    C_bigg=RHO_W*(PI/6.0_rk)**2*B_BIGG*(720.0_rk/4096.0_rk)
    do k = 1, nz
      T0k=T(k); qv0=qv(k); qc0=qc(k); qr0=qr(k)
      qi0=qi(k); qs0=qs(k); qg0=qg(k); nc0=nc(k); nr0=nr(k); ni0=ni(k)
      aqv=0.0_rk; aqc=0.0_rk; aqr=0.0_rk; aqi=0.0_rk; aqs=0.0_rk
      aqg=0.0_rk; anc=0.0_rk; anr=0.0_rk; ani=0.0_rk; aT=0.0_rk
      dqs_a=0.0_rk
      ! --- heat budget of a melting particle (T > 0 C): conduction from the
      !     air minus evaporative cooling of its wet surface when the air is
      !     below water saturation at 0 C. Positive = it melts (see below). ---
      heat = 0.0_rk
      if (T0k > T_FREEZE) heat = K_A*(T0k-T_FREEZE)                         &
             - L_VAP*D_V*rho(k)*max(qsat_liquid(T_FREEZE,P(k)) - qv0, 0.0_rk)
      ! --- deposition/sublimation on ice/snow/graupel ---
      !     Cloud ice only below 0 C (above, it melts at the end of the step,
      !     melt_cloud_ice). Above 0 C snow and graupel only sublimate, and
      !     only while they are not melting (a melting particle sits at 0 C),
      !     as in Thompson.
      dqi = 0.0_rk; dqs = 0.0_rk; dqg = 0.0_rk
      if (T0k < T_FREEZE) &
        dqi = dep_rate(qi0,Nb(k,II),Dm(k,II),lam(k,II),qv0,T0k,P(k),rho(k),QI_MIN,dt,II)
      if (T0k < T_FREEZE .or. heat <= 0.0_rk) then
        dqs = dep_rate(qs0,Nb(k,IS),Dm(k,IS),lam(k,IS),qv0,T0k,P(k),rho(k),QS_MIN,dt,IS)
        ! habit-dependent capacitance (see C_SNOW_*): dep_rate assumes the
        ! sphere value 0.5*D, so rescale to the temperature-appropriate habit.
        C_snow = C_SNOW_PLATE + ((T0k-T0)-T_C_SNOW_LO)*(C_SNOW_SPH-C_SNOW_PLATE) &
                                /(T_C_SNOW_HI-T_C_SNOW_LO)
        C_snow = min(max(C_snow, C_SNOW_PLATE), C_SNOW_SPH)
        dqs = dqs*(C_snow/C_SNOW_SPH)
        dqg = dep_rate(qg0,Nb(k,IG),Dm(k,IG),lam(k,IG),qv0,T0k,P(k),rho(k),QG_MIN,dt,IG)
      end if
      ! The three draw on one vapour reservoir: together they may at most
      ! bring the air to ice saturation within the step, latent heating
      ! included (Thompson caps the sum of its ice vapour terms likewise).
      qsi   = qsat_ice(T0k,P(k))
      dq_ex = (qv0 - qsi)/(1.0_rk + L_SUB**2*qsi/(RV*CP*T0k**2))
      dsum  = (dqi + dqs + dqg)*dt
      if ((dsum > 0.0_rk .and. dsum > dq_ex) .or. (dsum < 0.0_rk .and. dsum < dq_ex)) then
        scal = max(dq_ex/dsum, 0.0_rk)
        dqi = dqi*scal; dqs = dqs*scal; dqg = dqg*scal
      end if
      ! sublimating ice loses crystal NUMBER in proportion to mass
      if (dqi < 0.0_rk) ani = ani + ni0*dqi/max(qi0,Q_MIN)
      aqi=aqi+dqi; aqs=aqs+dqs; aqg=aqg+dqg
      aqv=aqv-(dqi+dqs+dqg); aT=aT+L_SUB/CP*(dqi+dqs+dqg)
      ! --- rain freezing (Bigg) -> graupel ---
      if (T0k < T_FREEZE .and. qr0 > QR_MIN) then
        dTf = max(T_FREEZE-T0k,0.0_rk)
        dqr_f = min(C_bigg*(exp(A_BIGG*dTf)-1.0_rk)*nr0*Dm(k,IR)**6, qr0/dt)
        aqr=aqr-dqr_f; aqg=aqg+dqr_f
        anr=anr - nr0*dqr_f/max(qr0,Q_MIN); aT=aT+L_FUS/CP*dqr_f
      end if
      ! --- Bigg (1953) immersion freezing of CLOUD DROPLETS -> ice ---
      !     Integrated over the cloud gamma PSD (shape nu_c = cloud_nu):
      !       dN/dt = B (exp(a dT)-1) (pi/6)   <D^3> N
      !       dq/dt = B (exp(a dT)-1) (pi/6)^2 (rho_w/rho) <D^6> N
      !     The D^6 weighting keeps this negligible above ~-25 C and makes it
      !     dominant by -35 C. It is the primary ice source that replaces the
      !     ungated DeMott in warm mixed-phase cloud.
      if (T0k < T_FREEZE .and. qc0 > QC_MIN) then
        lamc  = max(lam(k,IC), 1.0_rk)
        dTf   = max(T_FREEZE-T0k, 0.0_rk)
        fac_b = B_BIGG*(exp(A_BIGG*dTf)-1.0_rk)
        nu_b  = cloud_nu(nc0*rho(k))
        dqc_f = fac_b*(PI/6.0_rk)**2*(RHO_W/rho(k))*(nc0*rho(k))*CL_G7_G1(nu_b)/lamc**6
        dNc_f = fac_b*(PI/6.0_rk)*(nc0*rho(k))*BR_G2_G1(nu_b)/lamc**3/rho(k)
        dqc_f = min(dqc_f, qc0/dt)
        dNc_f = min(dNc_f, nc0/dt)
        aqc=aqc-dqc_f; aqi=aqi+dqc_f
        anc=anc-dNc_f; ani=ani+dNc_f
        aT=aT+L_FUS/CP*dqc_f
      end if
      ! --- ice -> snow autoconversion: only crystals grown past a snow-size
      !     threshold (Dm > D_I2S) convert; small pristine ice PERSISTS as
      !     suspended cloud ice (sustains IWP). Rate E_agg(T)*qi/TAU_AGG with a
      !     T-dependent sticking efficiency (weak in the cold). The old
      !     0.025*Nb(II)*qi form treated Nb (~1e4/m3) as a rate -> ~250/s, i.e.
      !     it flushed ALL cloud ice to snow every step whenever snow existed
      !     -> negligible cloud ice and near-zero ice water path. ---
      if (T0k < T_FREEZE .and. qi0 > QI_MIN .and. Dm(k,II) > D_I2S) then
        E_agg = min(max(exp(0.09_rk*(T0k-T_FREEZE)), 0.05_rk), 1.0_rk)
        dqs_a = min(E_agg*qi0/TAU_AGG, qi0/dt)
        aqs=aqs+dqs_a; aqi=aqi-dqs_a
        ani=ani - ni0*dqs_a/max(qi0,Q_MIN)
      end if
      ! --- melting: snow/graupel -> rain (T > 0 C) ---
      !     The heat budget above, with ventilation integrated over the
      !     spectrum (Rutledge & Hobbs 1983; Thompson prr_sml/prr_gml):
      !       dq/dt = -4 pi C' heat * int(D f_v n dD) / (rho L_f),
      !       heat  = K_a Tc - L_v D_v rho (q_sw(0 C) - q_v),
      !     C' = 0.15 for wet snow, 0.5 for graupel. Each melted particle
      !     becomes one rain drop: Nr gains N_x/q_x per unit mass melted.
      if (T0k > T_FREEZE .and. heat > 0.0_rk) then
        rate_s=0.0_rk; rate_g=0.0_rk
        if (qs0 > QS_MIN) rate_s = min(4.0_rk*PI*C_MELT_S*heat                 &
             *vent_moment(Nb(k,IS),lam(k,IS),IS,rho(k),FV1_MELT,FV2_MELT)        &
             /(rho(k)*L_FUS), qs0/dt)
        if (qg0 > QG_MIN) rate_g = min(4.0_rk*PI*C_MELT_G*heat                 &
             *vent_moment(Nb(k,IG),lam(k,IG),IG,rho(k),FV1_MELT,FV2_MELT)        &
             /(rho(k)*L_FUS), qg0/dt)
        aqs=aqs-rate_s; aqg=aqg-rate_g; aqr=aqr+(rate_s+rate_g)
        anr=anr + ( Nb(k,IS)*rate_s/max(qs0,Q_MIN)                              &
                  + Nb(k,IG)*rate_g/max(qg0,Q_MIN) )/rho(k)
        aT=aT - L_FUS/CP*(rate_s+rate_g)
      end if
      ! --- conserving limiter (see limit_scale) ---
      scal = limit_scale(aT, dt, (/qv0,qc0,qr0,qi0,qs0,qg0/), (/aqv,aqc,aqr,aqi,aqs,aqg/))
      if (scal < 1.0_rk) then
        aqv=aqv*scal; aqc=aqc*scal; aqr=aqr*scal; aqi=aqi*scal
        aqs=aqs*scal; aqg=aqg*scal; anc=anc*scal
        anr=anr*scal; ani=ani*scal; aT=aT*scal
        dqs_a=dqs_a*scal
      end if
      agg(k)=dqs_a
      ! --- apply once ---
      qv(k)=max(qv0+aqv*dt,0.0_rk); qc(k)=max(qc0+aqc*dt,0.0_rk)
      qr(k)=max(qr0+aqr*dt,0.0_rk); qi(k)=max(qi0+aqi*dt,0.0_rk)
      qs(k)=max(qs0+aqs*dt,0.0_rk); qg(k)=max(qg0+aqg*dt,0.0_rk)
      nc(k)=max(nc0+anc*dt,0.0_rk)
      nr(k)=max(nr0+anr*dt,0.0_rk); ni(k)=max(ni0+ani*dt,0.0_rk)
      T(k) =T0k+aT*dt
    end do
  contains
    ! dep_rate: diffusion-limited deposition(+)/sublimation(-) rate
    ! [kg kg-1 s-1] for one ice category c. Combines supersaturation over
    ! ice (Si-1), the particle population (Nbx,Dmx), and a ventilation
    ! enhancement fv from the fall-speed Reynolds number evaluated at the
    ! mean size. Sublimation is capped at the ice present; deposition is
    ! capped by the caller, jointly for the three categories.
    real(rk) function dep_rate(qx,Nbx,Dmx,lamx,qvx,Tx,Px,rhx,qmin,dtl,c)
      real(rk), intent(in) :: qx,Nbx,Dmx,lamx,qvx,Tx,Px,rhx,qmin,dtl
      integer,  intent(in) :: c
      real(rk) :: Vf,Re,xv,fv,Si,esi,Ad,Bd,den,r
      ! ventilation Reynolds number from the number-weighted fall speed and
      ! the mass-weighted mean diameter
      Vf = av(c)*GN(c)/max(lamx,1.0_rk)**bv(c)*sqrt(RHO_NOT/rhx)
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
      ! cap sublimation at the ice present
      if (r < 0.0_rk) r = max(r, -qx/dtl)
      dep_rate = r
    end function
  end subroutine cold

  !-----------------------------------------------------------------
  !-----------------------------------------------------------------
  ! collection  
  !   Mixed-phase collision-collection (riming and wet growth), active
  !   only below 0 C. A faster-falling "collector" sweeps up a slower
  !   "collected" species; the differential fall speed sets the rate:
  !     1 cloud droplets by cloud ice -> snow (the rime; crystals stay ice)
  !     2 cloud droplets by snow      -> snow
  !     3 cloud droplets by graupel   -> graupel
  !     4 rain + cloud ice            -> graupel
  !     5 rain by snow                -> graupel (snow is the collector)
  !   plus Hallett-Mossop splintering driven by the riming of 2 and 3.
  !   Each path freezes liquid onto ice, so it warms the air (L_fus).
  !-----------------------------------------------------------------
  subroutine collection(nz, T,rho, qc,qr,qi,qs,qg, nc,nr,ni, lam,Dm,Nb, dt)
    integer,  intent(in)    :: nz
    real(rk), intent(in)    :: rho(nz),lam(nz,NCAT),Dm(nz,NCAT),Nb(nz,NCAT),dt
    real(rk), intent(inout) :: T(nz),qc(nz),qr(nz),qi(nz),qs(nz),qg(nz),nc(nz),nr(nz),ni(nz)
    integer  :: k, nu_k
    real(rk) :: T0k,qc0,qr0,qi0,qs0,qg0,nc0,nr0,ni0
    real(rk) :: aqc,aqr,aqi,aqs,aqg,anc,anr,ani,aT
    real(rk) :: Vrel,rate,rate_qi,rate_qr,scal,sweep
    real(rk) :: E_cs,E_cg,vtc,vts,vtg
    real(rk) :: rate_s,rate_g,Tc,f_hm,dNi_h,dqi_h
    do k = 1, nz
      if (T(k) >= T_FREEZE) cycle
      T0k=T(k); qc0=qc(k); qr0=qr(k); qi0=qi(k); qs0=qs(k); qg0=qg(k)
      nc0=nc(k); nr0=nr(k); ni0=ni(k)
      nu_k = cloud_nu(nc0*rho(k))            ! droplet shape for the cloud fall speed
      vtc  = cloud_vq(lam(k,IC),nu_k,rho(k))
      aqc=0.0_rk; aqr=0.0_rk; aqi=0.0_rk; aqs=0.0_rk; aqg=0.0_rk
      anc=0.0_rk; anr=0.0_rk; ani=0.0_rk; aT=0.0_rk
      ! --- 1. riming of droplets by cloud ice: qc -> qs ---
      !     The rime is counted as (small rimed) snow; the collecting crystals
      !     keep their mass, so cloud-ice mass and number are unchanged.
      if (qc0 > QC_MIN .and. qi0 > QI_MIN) then
        Vrel = abs(fall_speed_q(lam(k,II),II,rho(k)) - vtc)
        rate = min(max(E_CI*(PI/4.0_rk)*Dm(k,II)**2*Nb(k,II)*Vrel*qc0,0.0_rk), qc0/dt)
        anc=anc + nc0*(-rate)/max(qc0,Q_MIN)
        aqc=aqc-rate; aqs=aqs+rate; aT=aT+L_FUS/CP*rate
      end if
      ! --- 2. riming of droplets by snow: qc -> qs (riming grows snow) ---
      !     Efficiency from the Stokes number (Wang & Ji 2000): exactly zero
      !     for droplets < 6 um or snow < 300 um.
      rate_s = 0.0_rk
      if (qc0 > QC_MIN .and. qs0 > QS_MIN) then
        vts  = fall_speed_q(lam(k,IS),IS,rho(k))
        Vrel = abs(vts - vtc)
        E_cs = Ef_snow_cloud(Dm(k,IC), Dm(k,IS), vts)
        rate_s = min(max(E_cs*(PI/4.0_rk)*Dm(k,IS)**2*Nb(k,IS)*Vrel*qc0,0.0_rk), qc0/dt)
        anc=anc + nc0*(-rate_s)/max(qc0,Q_MIN)
        aqc=aqc-rate_s; aqs=aqs+rate_s; aT=aT+L_FUS/CP*rate_s
      end if
      ! --- 3. riming of droplets by graupel: qc -> qg ---
      !     Stokes-number efficiency (zero below St = 0.4).
      rate_g = 0.0_rk
      if (qc0 > QC_MIN .and. qg0 > QG_MIN) then
        vtg  = fall_speed_q(lam(k,IG),IG,rho(k))
        Vrel = abs(vtg - vtc)
        E_cg = Ef_graupel_cloud(Dm(k,IC), Dm(k,IG), vtg)
        rate_g = min(max(E_cg*(PI/4.0_rk)*Dm(k,IG)**2*Nb(k,IG)*Vrel*qc0,0.0_rk), qc0/dt)
        anc=anc + nc0*(-rate_g)/max(qc0,Q_MIN)
        aqc=aqc-rate_g; aqg=aqg+rate_g; aT=aT+L_FUS/CP*rate_g
      end if
      ! --- Hallett-Mossop rime splintering (secondary ice) ---
      !     HM_SPL splinters per kg of rime collected by snow and graupel
      !     (Hallett & Mossop 1974), weighted 0 at -8 C, 1 at -5 C, 0 at -3 C
      !     (Thompson's tf). Each splinter is a seed crystal of M_ICE_MIN taken
      !     from the rimers, so water is conserved and no latent heat is due.
      Tc = T0k - T_FREEZE
      if (Tc > -8.0_rk .and. Tc < -3.0_rk .and. rate_s + rate_g > 0.0_rk) then
        if (Tc >= -5.0_rk) then
          f_hm = 0.5_rk*(-3.0_rk - Tc)
        else
          f_hm = (8.0_rk + Tc)/3.0_rk
        end if
        dNi_h = HM_SPL*f_hm*(rate_s + rate_g)
        dqi_h = M_ICE_MIN*dNi_h
        ani = ani + dNi_h; aqi = aqi + dqi_h
        aqs = aqs - dqi_h*rate_s/(rate_s + rate_g)
        aqg = aqg - dqi_h*rate_g/(rate_s + rate_g)
      end if
      ! --- 4. rain-ice collection: qr + qi -> qg ---
      !     A drop and a crystal that collide freeze together. The collision
      !     rate per volume, E pi/4 D_r^2 |dV| N_r N_i, removes each crystal at
      !     the rate E pi/4 D_r^2 |dV| N_r and each drop at E pi/4 D_r^2 |dV| N_i.
      if (qi0 > QI_MIN .and. qr0 > QR_MIN) then
        Vrel  = abs(fall_speed_q(lam(k,IR),IR,rho(k)) - fall_speed_q(lam(k,II),II,rho(k)))
        sweep = E_RI*(PI/4.0_rk)*Dm(k,IR)**2*Vrel
        rate_qi = min(sweep*Nb(k,IR)*qi0, qi0/dt)
        rate_qr = min(sweep*Nb(k,II)*qr0, qr0/dt)
        ani=ani - ni0*rate_qi/max(qi0,Q_MIN)
        anr=anr - nr0*rate_qr/max(qr0,Q_MIN)
        aqi=aqi-rate_qi; aqr=aqr-rate_qr; aqg=aqg+(rate_qi+rate_qr)
        aT=aT+L_FUS/CP*rate_qr
      end if
      ! --- 5. snow collects rain: qr -> qg (rain freezes; snow is collector) ---
      if (qr0 > QR_MIN .and. qs0 > QS_MIN) then
        Vrel = abs(fall_speed_q(lam(k,IR),IR,rho(k)) - fall_speed_q(lam(k,IS),IS,rho(k)))
        rate_qr = min(max(E_SR*(PI/4.0_rk)*Dm(k,IS)**2*Nb(k,IS)*Vrel*qr0,0.0_rk), qr0/dt)
        anr=anr + nr0*(-rate_qr)/max(qr0,Q_MIN)
        aqr=aqr-rate_qr; aqg=aqg+rate_qr
        aT=aT+L_FUS/CP*rate_qr
      end if
      ! --- conserving limiter (see limit_scale) ---
      scal = limit_scale(aT, dt, (/qc0,qr0,qi0,qs0,qg0/), (/aqc,aqr,aqi,aqs,aqg/))
      if (scal < 1.0_rk) then
        aqc=aqc*scal; aqr=aqr*scal; aqi=aqi*scal; aqs=aqs*scal; aqg=aqg*scal
        anc=anc*scal; anr=anr*scal; ani=ani*scal; aT=aT*scal
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
  ! Ef_snow_cloud
  !   Snow-collects-droplet efficiency, Wang & Ji (2000), transliterated from
  !   Thompson's table_Efsw (module_mp_thompson.F:5218) but evaluated at the
  !   mean diameters instead of tabulated over the spectrum. Returns exactly
  !   zero when the droplet is smaller than D_C_RIME, the snow smaller than
  !   D_S_RIME, the size ratio p > 0.25, or the snow is not falling.
  !     Dc  mean cloud-droplet diameter [m]
  !     Ds  mean snow diameter          [m]
  !     vts snow fall speed             [m s-1]
  !-----------------------------------------------------------------
  real(rk) function Ef_snow_cloud(Dc, Ds, vts)
    real(rk), intent(in) :: Dc, Ds, vts
    real(rk) :: Ds_m, p, stokes, reyn, F, G, K0, z, H, yc0, Ef
    Ef_snow_cloud = 0.0_rk
    if (Ds < D_S_RIME .or. Dc < D_C_RIME .or. vts < 1.0e-3_rk) return
    ! melted-equivalent snow diameter (Wang & Ji's effective cross-section)
    Ds_m = (am(IS)*Ds**bm(IS)/am(IR))**(1.0_rk/3.0_rk)
    p    = Dc/max(Ds_m,1.0e-12_rk)
    if (p > 0.25_rk) return
    p    = max(p, 1.0e-6_rk)
    stokes = max(Dc*Dc*vts*RHO_W/(9.0_rk*MU_AIR*max(Ds_m,1.0e-12_rk)), 1.0e-12_rk)
    reyn   = max(9.0_rk*stokes/(p*p*RHO_W), 1.0e-12_rk)
    F   = log(reyn)
    G   = -0.1007_rk - 0.358_rk*F + 0.0261_rk*F*F
    K0  = exp(G)
    z   = log(stokes/(K0+1.0e-15_rk))
    H   = 0.1465_rk + 1.302_rk*z - 0.607_rk*z*z + 0.293_rk*z*z*z
    yc0 = 2.0_rk/PI*atan(H)
    Ef  = (yc0+p)**2/(1.0_rk+p)**2
    Ef_snow_cloud = min(max(Ef, 0.0_rk), 0.95_rk)
  end function Ef_snow_cloud

  !-----------------------------------------------------------------
  ! Ef_graupel_cloud
  !   Graupel-collects-droplet efficiency from the Stokes number
  !   (module_mp_thompson.F:2548): zero below St = 0.4, 0.55*log10(2.51 St)
  !   between 0.4 and 10, and 0.77 above.
  !-----------------------------------------------------------------
  real(rk) function Ef_graupel_cloud(Dc, Dg, vtg)
    real(rk), intent(in) :: Dc, Dg, vtg
    real(rk) :: St
    St = Dc*Dc*vtg*RHO_W/(9.0_rk*MU_AIR*max(Dg,1.0e-12_rk))
    if (St < 0.4_rk) then
      Ef_graupel_cloud = 0.0_rk
    else if (St > 10.0_rk) then
      Ef_graupel_cloud = 0.77_rk
    else
      Ef_graupel_cloud = 0.55_rk*log10(max(2.51_rk*St, 1.0e-12_rk))
    end if
    Ef_graupel_cloud = min(max(Ef_graupel_cloud, 0.0_rk), 0.95_rk)
  end function Ef_graupel_cloud

  !-----------------------------------------------------------------
  !-----------------------------------------------------------------
  ! sedimentation  
  !   Gravitational settling of every category through the column via the
  !   upwind scheme above. 2-moment species (cloud,rain,ice) settle mass
  !   and number separately (mass- vs number-weighted fall speed); snow
  !   and graupel settle mass only. Fall speeds carry the air-density
  !   factor sqrt(RHO_NOT/rho). Number is passed to upstream_flux as a
  !   specific quantity [kg-1], like the mixing ratios, so the scheme
  !   conserves it exactly. Returns the per-step settling tendencies for
  !   rain/ice/snow/graupel (dq*_sed), which the column routine integrates
  !   over the column to get surface precip.
  !-----------------------------------------------------------------
  subroutine sedimentation(nz, rho,dz,dt, qc,qr,qi,qs,qg, nc,nr,ni, lam, &
                           dqr_sed,dqi_sed,dqs_sed,dqg_sed, flr,fli,fls,flg)
    integer,  intent(in)    :: nz
    real(rk), intent(in)    :: rho(nz),dz(nz),dt,lam(nz,NCAT)
    real(rk), intent(inout) :: qc(nz),qr(nz),qi(nz),qs(nz),qg(nz),nc(nz),nr(nz),ni(nz)
    real(rk), intent(out)   :: dqr_sed(nz),dqi_sed(nz),dqs_sed(nz),dqg_sed(nz)
    ! optional precipitation fluxes out of each level [kg m-2 s-1]
    real(rk), intent(out), optional :: flr(nz),fli(nz),fls(nz),flg(nz)
    real(rk) :: Vq(nz),Vn(nz),dq(nz),dN(nz)
    integer  :: k, nu
    ! cloud (2-moment)
    do k=1,nz
      nu = cloud_nu(nc(k)*rho(k))
      Vq(k)=cloud_vq(lam(k,IC),nu,rho(k)); Vn(k)=cloud_vn(lam(k,IC),nu,rho(k))
    end do
    call upstream_flux(nz,qc,Vq,rho,dz,dt,dq); qc=max(qc+dq*dt,0.0_rk)
    call upstream_flux(nz,nc,Vn,rho,dz,dt,dN); nc=max(nc+dN*dt,0.0_rk)
    ! rain
    do k=1,nz
      Vq(k)=fall_speed_q(lam(k,IR),IR,rho(k)); Vn(k)=fall_speed_N(lam(k,IR),IR,rho(k))
    end do
    call upstream_flux(nz,qr,Vq,rho,dz,dt,dqr_sed,flr); qr=max(qr+dqr_sed*dt,0.0_rk)
    call upstream_flux(nz,nr,Vn,rho,dz,dt,dN); nr=max(nr+dN*dt,0.0_rk)
    ! ice
    do k=1,nz
      Vq(k)=fall_speed_q(lam(k,II),II,rho(k)); Vn(k)=fall_speed_N(lam(k,II),II,rho(k))
    end do
    call upstream_flux(nz,qi,Vq,rho,dz,dt,dqi_sed,fli); qi=max(qi+dqi_sed*dt,0.0_rk)
    call upstream_flux(nz,ni,Vn,rho,dz,dt,dN); ni=max(ni+dN*dt,0.0_rk)
    ! snow (1-moment)
    do k=1,nz; Vq(k)=fall_speed_q(lam(k,IS),IS,rho(k)); end do
    call upstream_flux(nz,qs,Vq,rho,dz,dt,dqs_sed,fls); qs=max(qs+dqs_sed*dt,0.0_rk)
    ! graupel (1-moment)
    do k=1,nz; Vq(k)=fall_speed_q(lam(k,IG),IG,rho(k)); end do
    call upstream_flux(nz,qg,Vq,rho,dz,dt,dqg_sed,flg); qg=max(qg+dqg_sed*dt,0.0_rk)
  end subroutine sedimentation

  !-----------------------------------------------------------------
  !-----------------------------------------------------------------
  ! melt_cloud_ice
  !   Cloud ice cannot survive above 0 C: whatever sedimented or was
  !   advected into such a cell melts into cloud droplets (one droplet per
  !   crystal), as in Thompson's end-of-step clean-up ("Instantly melt any
  !   cloud ice into cloud water if above 0C"). The melt is limited by the
  !   heat available (the cell is not cooled below 0 C) and by DT_MAX like
  !   every other process; any rest melts on the next step.
  !-----------------------------------------------------------------
  subroutine melt_cloud_ice(nz, T, qc, qi, nc, ni)
    integer,  intent(in)    :: nz
    real(rk), intent(inout) :: T(nz), qc(nz), qi(nz), nc(nz), ni(nz)
    integer  :: k
    real(rk) :: f
    do k = 1, nz
      if (T(k) > T_FREEZE .and. qi(k) > 0.0_rk) then
        f = min(1.0_rk, min(DT_MAX, T(k)-T_FREEZE)*CP/(L_FUS*qi(k)))
        qc(k) = qc(k) + f*qi(k);  nc(k) = nc(k) + f*ni(k)
        T(k)  = T(k) - L_FUS/CP*f*qi(k)
        qi(k) = (1.0_rk-f)*qi(k); ni(k) = (1.0_rk-f)*ni(k)
      end if
    end do
  end subroutine melt_cloud_ice

  !-----------------------------------------------------------------
  !-----------------------------------------------------------------
  ! saturation_adjustment
  !   Final liquid-phase cleanup: where T>T_HOM (i.e. wherever liquid water
  !   can exist, SUPERCOOLED LIQUID INCLUDED, not just T>=0 C), condense excess
  !   vapour to cloud (or evaporate cloud) so qv relaxes to qsat_liquid,
  !   updating T for the latent heat. Three Newton iterations handle qsat's own
  !   T-dependence. Running this below 0 C lets supercooled liquid cloud form
  !   (pins RHw<=1 and feeds the Bergeron process via cold's deposition); the
  !   old T>=0 C guard suppressed supercooled liquid and let RHw/RHi run away.
  !   Below T_HOM there is no liquid: homogeneous
  !   freezing has already glaciated any qc. Between T_HOM and T_HOM+3 K no
  !   new cloud water is made in a droplet-free cell (activation is off there).
  !   The ice phase is deliberately NOT adjusted here -- it is left to the
  !   diffusion-limited rates in `cold`.
  !-----------------------------------------------------------------
  subroutine saturation_adjustment(nz, T,P, qv,qc,qi,qs,qg, nc)
    integer,  intent(in)    :: nz
    real(rk), intent(in)    :: P(nz)
    real(rk), intent(inout) :: T(nz),qv(nz),qc(nz),qi(nz),qs(nz),qg(nz),nc(nz)
    integer  :: k, it
    real(rk) :: ql,excess,denom,dq,qc_old,dq_max,dq_cum
    dq_max = DT_MAX*CP/L_VAP        ! max |cumulative| condensation for DT_MAX
    do k = 1, nz
      dq_cum = 0.0_rk
      do it = 1, 3
        if (T(k) > T_HOM) then
          ql = qsat_liquid(T(k),P(k))
          excess = qv(k) - ql
          denom  = 1.0_rk + L_VAP**2*ql/(RV*CP*T(k)**2)
          dq = excess/denom
          dq = max(dq, -qc(k))
          ! no NEW cloud water without droplets in the homogeneous-freezing
          ! band (T <= T_HOM+3, -35..-38 C), where CCN activation is off: it
          ! would be numberless liquid. The vapour is left to the ice
          ! processes; existing droplets still grow or evaporate.
          if (T(k) <= T_HOM+3.0_rk .and. nc(k) <= 0.0_rk) dq = min(dq, 0.0_rk)
          ! latent-heat limiter: cap cumulative condensation so the whole
          ! adjustment changes T by at most DT_MAX (residual relaxes next step)
          dq = min(max(dq, -dq_max-dq_cum), dq_max-dq_cum)
          dq_cum = dq_cum + dq
          qc_old = qc(k)
          qc(k) = qc(k) + dq
          qv(k) = qv(k) - dq
          T(k)  = T(k) + L_VAP/CP*dq
          ! shed droplet number in proportion to cloud water evaporated (dq<0);
          ! condensation only grows existing droplets, so nc is unchanged.
          ! Stops activation-created nc from being stranded (as spurious cloud
          ! droplet number in cloud-free air) where the cloud has evaporated.
          if (dq < 0.0_rk .and. qc_old > Q_MIN) &
               nc(k) = nc(k)*max(qc(k),0.0_rk)/qc_old
        end if
        qc(k)=max(qc(k),0.0_rk); qi(k)=max(qi(k),0.0_rk)
        qs(k)=max(qs(k),0.0_rk); qg(k)=max(qg(k),0.0_rk); qv(k)=max(qv(k),0.0_rk)
      end do
      ! number<->mass consistency: no droplets where there is no cloud water
      if (qc(k) < QC_MIN) nc(k) = 0.0_rk
    end do
  end subroutine saturation_adjustment

  !-----------------------------------------------------------------
  ! column_diag  
  !   Optional radiation/radar diagnostics for one column:
  !     re_qc,re_qi,re_qs  effective radii [m] (clamped to WRF bounds)
  !                        Reff = (mu+3)/(2*lam); re_qs uses the same form as
  !                        cloud & ice so RRTMG also gets a snow radius.
  !     dbz                equivalent 10-cm reflectivity [dBZ], Rayleigh: the
  !                        linear sum of the liquid 6th moments and, for
  !                        the frozen categories, of the reflectivity of
  !                        solid-ice spheres of the same mass, scaled by
  !                        |K_ice|^2/|K_water|^2 (Smith 1984).
  !-----------------------------------------------------------------
  subroutine column_diag(nz, T,rho, qc,qr,qi,qs,qg, nc,nr,ni, &
                         re_qc,re_qi,re_qs, dbz)
    integer,  intent(in)  :: nz
    real(rk), intent(in)  :: T(nz),rho(nz)
    real(rk), intent(in)  :: qc(nz),qr(nz),qi(nz),qs(nz),qg(nz),nc(nz),nr(nz),ni(nz)
    real(rk), intent(out) :: re_qc(nz),re_qi(nz),re_qs(nz),dbz(nz)
    real(rk) :: lam(nz,NCAT),Dm(nz,NCAT),Nb(nz,NCAT),N0s(nz),N0g(nz)
    real(rk) :: ze
    integer  :: k, nu
    call diagnose(nz,T,rho, qc,qr,qi,qs,qg, nc,nr,ni, lam,Dm,Nb,N0s,N0g)
    do k = 1, nz
      nu = cloud_nu(nc(k)*rho(k))
      re_qc(k) = min(max((real(nu,rk)+3.0_rk)/(2.0_rk*lam(k,IC)), REQC_MIN), REQC_MAX)
      re_qi(k) = min(max((mu(II)+3.0_rk)/(2.0_rk*lam(k,II)), REQI_MIN), REQI_MAX)
      re_qs(k) = min(max((mu(IS)+3.0_rk)/(2.0_rk*lam(k,IS)), REQS_MIN), REQS_MAX)
      ze = 1.0e18_rk*Nb(k,IC)*CL_G7_G1(nu)/lam(k,IC)**6  &
         + 1.0e18_rk*Nb(k,IR)*G2B(IR)/lam(k,IR)**6       &
         + ze_ice(Nb(k,II),lam(k,II),II)                  &
         + ze_ice(Nb(k,IS),lam(k,IS),IS)                  &
         + ze_ice(Nb(k,IG),lam(k,IG),IG)
      dbz(k) = 10.0_rk*log10(max(ze, 1.0e-4_rk))
    end do
  contains
    ! linear equivalent reflectivity [mm6 m-3] of a frozen category: a
    ! particle of mass m = am*D^bm reflects like the solid-ice sphere of
    ! the same mass, D_i^6 = (6 m/(pi rho_i))^2, so
    !   Ze = |K_i|^2/|K_w|^2 (6 am/(pi rho_i))^2 int D^(2 bm) n dD
    ! (for snow this uses its own mass law, not a fixed density).
    real(rk) function ze_ice(Nbx,lamx,c)
      real(rk), intent(in) :: Nbx,lamx
      integer,  intent(in) :: c
      ze_ice = 1.0e18_rk*(K2_ICE/K2_LIQ)*(6.0_rk*am(c)/(PI*RHO_I))**2 &
               *Nbx*G2B(c)/lamx**(2.0_rk*bm(c))
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
       wetscav_on, rainprod, evapprod, qlsink,                           &
       precr, preci, precs, precg,                                       &
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
    ! ---- optional WRF-Chem wet-scavenging outputs (as Thompson's) ----
    logical, optional, intent(in) :: wetscav_on
    real, dimension(ims:ime,kms:kme,jms:jme), optional, intent(inout) :: &
                          rainprod, evapprod, qlsink, precr, preci, precs, precg

    ! ---- locals : DOUBLE precision column work arrays ----
    integer  :: i,j,k, nz, i_end, j_end
    real(rk) :: dt
    real(rk) :: t1d(kts:kte), p1d(kts:kte), w1d(kts:kte), dz1d(kts:kte), rho(kts:kte)
    real(rk) :: qv1d(kts:kte),qc1d(kts:kte),qr1d(kts:kte),qi1d(kts:kte),qs1d(kts:kte),qg1d(kts:kte)
    real(rk) :: nc1d(kts:kte),nr1d(kts:kte),ni1d(kts:kte),nwfa1d(kts:kte),nifa1d(kts:kte)
    real(rk) :: re_qc(kts:kte), re_qi(kts:kte), re_qs(kts:kte), dbz(kts:kte)
    real(rk) :: c2p1d(kts:kte), rp1d(kts:kte), ep1d(kts:kte)
    real(rk) :: pr1d(kts:kte), pi1d(kts:kte), ps1d(kts:kte), pg1d(kts:kte)
    logical  :: want_scav
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
    want_scav = present(rainprod) .and. present(evapprod) .and. present(qlsink) .and. &
                present(precr) .and. present(preci) .and. present(precs) .and. present(precg)
    if (present(wetscav_on)) want_scav = want_scav .and. wetscav_on

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
        if (want_scav) then
          call mp_microwcp_column(nz, t1d,p1d,rho,w1d,dz1d, &
               qv1d,qc1d,qr1d,qi1d,qs1d,qg1d, nc1d,nr1d,ni1d, nwfa1d,nifa1d, &
               has_aero, dt, pptrain,pptsnow,pptgraul,pptice, &
               c2p1d, rp1d, ep1d, pr1d, pi1d, ps1d, pg1d)
        else
          call mp_microwcp_column(nz, t1d,p1d,rho,w1d,dz1d, &
               qv1d,qc1d,qr1d,qi1d,qs1d,qg1d, nc1d,nr1d,ni1d, nwfa1d,nifa1d, &
               has_aero, dt, pptrain,pptsnow,pptgraul,pptice)
        end if

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
          qc(i,k,j) = min(max(real(qc1d(k)),0.0), 0.03)
          qr(i,k,j) = min(max(real(qr1d(k)),0.0), 0.03)
          qi(i,k,j) = min(max(real(qi1d(k)),0.0), 0.03)
          qs(i,k,j) = min(max(real(qs1d(k)),0.0), 0.03)
          qg(i,k,j) = min(max(real(qg1d(k)),0.0), 0.03)
          nr(i,k,j) = min(max(real(nr1d(k)),0.0), 1.0e13)
          ni(i,k,j) = min(max(real(ni1d(k)),0.0), 1.0e13)
          th(i,k,j) = real(t1d(k)/real(pii(i,k,j),rk))
          if (has_nc)   nc(i,k,j)   = min(max(real(nc1d(k)),0.0), 1.0e13)
          if (has_aero) then
            nwfa(i,k,j) = max(real(nwfa1d(k)),0.0); nifa(i,k,j) = max(real(nifa1d(k)),0.0)
          end if
        end do

        ! ---- WRF-Chem wet scavenging: the fields Thompson provides, units as
        !      documented (qlsink = fractional cloud-water loss to precipitation,
        !      s-1; see mp_microwcp_column for the difference from Thompson) ----
        if (want_scav) then
          do k = kts, kte
            rainprod(i,k,j) = real(rp1d(k))
            evapprod(i,k,j) = real(ep1d(k))
            precr(i,k,j)    = real(pr1d(k))
            preci(i,k,j)    = real(pi1d(k))
            precs(i,k,j)    = real(ps1d(k))
            precg(i,k,j)    = real(pg1d(k))
            if (qc1d(k) > 1.0e-10_rk) then
              qlsink(i,k,j) = real(min(1.0_rk, max(0.0_rk, c2p1d(k)/qc1d(k))))
            else
              qlsink(i,k,j) = 0.0
            end if
          end do
        end if

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
