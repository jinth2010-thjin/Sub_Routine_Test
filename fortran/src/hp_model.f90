!===============================================================================
! hp_model : 축방향 그루브(axial groove) 암모니아 히트파이프 물리 모델
!            CCHP (정전도 히트파이프)와 VCHP (가변전도 히트파이프, NCG flat-front)
!
! 열해석 솔버(노드-전도체 네트워크)와의 연결 개념
!
!   [증발부 벽 노드 i] --G_evap(i)--+                  +--G_cond(j)*f(j)-- [응축부 벽 노드 j]
!                                   [증기노드 Tve]--G_vap--[증기노드 Tvc]
!                                                              (NCG 저장조 노드 Tres)
!
!   - 솔버가 에너지 보존(노드 온도)을 풀고, 이 모듈은 매 반복마다 전도체 값만 갱신한다.
!   - G_evap, G_cond : 벽 반경방향 전도 + 그루브 액체층 전도 (증발/응축 막 저항)
!   - G_vap          : 증기 유동 압력강하 -> Clausius-Clapeyron 온도강하
!   - f(j)           : VCHP 응축부 구간 활성 비율 (NCG 가 막은 구간은 0)
!   - 수송 한계(모세관/음속/비말동반/비등/점성)를 넘으면 증발부 전도체를 축소해
!     수송 열량을 Q_max 로 제한 (증발부 윅 dry-out 근사). 응축부는 축소하지 않으므로
!     증기 온도는 응축부 쪽에 머물고, 남는 열은 증발부(장비) 온도 상승으로 나타난다.
!
! 단위: SI (m, K, W, Pa, kg, mol)
!===============================================================================
module hp_model
  use hp_props
  implicit none
  private

  real(dp), parameter :: PI = 3.14159265358979323846_dp
  real(dp), parameter :: RE_TRANS = 2300.0_dp

  integer, parameter, public :: LIM_NONE = 0, LIM_CAPILLARY = 1, LIM_SONIC = 2, &
                                LIM_ENTRAIN = 3, LIM_BOILING = 4, LIM_VISCOUS = 5
  character(len=10), parameter, public :: LIMIT_NAME(0:5) = [character(len=10) :: &
       'none', 'capillary', 'sonic', 'entrain', 'boiling', 'viscous']

  !-----------------------------------------------------------------------------
  ! 입력(설정) : 형상 / 운용 / VCHP / 수치 파라미터
  !-----------------------------------------------------------------------------
  type, public :: hp_config
     ! --- 길이 [m]
     real(dp) :: L_evap = 0.30_dp
     real(dp) :: L_adia = 0.40_dp
     real(dp) :: L_cond = 0.60_dp
     ! --- 단면 (원형 튜브 + 사각 축방향 그루브로 단순화) [m]
     real(dp) :: r_vap    = 4.0e-3_dp    ! 증기 코어 반경
     real(dp) :: groove_w = 0.6e-3_dp    ! 그루브 폭
     real(dp) :: groove_d = 1.2e-3_dp    ! 그루브 깊이
     integer  :: n_groove = 24           ! 그루브 개수
     real(dp) :: r_out    = 6.35e-3_dp   ! 튜브 외반경
     real(dp) :: k_wall   = 200.0_dp     ! 벽 열전도도 [W/m.K] (Al 6063)
     ! --- 운용 조건
     real(dp) :: contact_angle = 0.0_dp  ! 접촉각 [rad]
     real(dp) :: elevation = 0.0_dp      ! 증발부가 응축부보다 높은 높이 [m] (지상시험 역경사, 궤도 0)
     real(dp) :: g_accel = 9.80665_dp    ! 중력가속도 [m/s2]
     real(dp) :: r_nuc = 2.54e-7_dp      ! 비등 한계용 기포핵 반경 [m]
     ! 비등 한계: 고전식(Chi)은 축방향 그루브에서 매우 보수적이라 기본은 끔.
     ! 제작사가 증발부 최대 열유속을 주면 q_flux_max [W/m2, 외표면 기준] 로 대신 적용
     logical  :: use_boiling_formula = .false.
     real(dp) :: q_flux_max = 0.0_dp
     ! --- 시험 보정 계수 (열진공시험 상관관계용)
     real(dp) :: mult_evap = 1.0_dp
     real(dp) :: mult_cond = 1.0_dp
     ! --- VCHP (n_gas = 0 이면 CCHP)
     real(dp) :: n_gas = 0.0_dp          ! NCG 충전량 [mol]
     real(dp) :: v_res = 0.0_dp          ! 저장조 체적 [m3]
     real(dp) :: blocked_frac = 0.0_dp   ! NCG 차단 구간의 잔류 전도 비율 (확산 등, 0~1)
     ! --- 수치
     real(dp) :: relax = 0.5_dp          ! 전도체 갱신 완화계수 (0<relax<=1)
     real(dp) :: g_vap_max = 1.0e5_dp    ! 증기 전도체 상한 [W/K] (행렬 조건수 보호)
     logical  :: limit_transport = .true.! 수송 한계 적용 여부
  end type hp_config

  !-----------------------------------------------------------------------------
  ! 상태(결과) : 전도체 값과 진단 정보. 히트파이프 1개당 1개 유지
  !-----------------------------------------------------------------------------
  type, public :: hp_state
     integer :: ne = 0, nc = 0
     real(dp), allocatable :: dx_e(:), dx_c(:)     ! 구간 길이
     real(dp), allocatable :: g_evap(:), g_cond(:) ! 전도체 [W/K]
     real(dp), allocatable :: f_act(:)             ! 응축부 활성 비율 (완화 적용 후)
     real(dp), allocatable :: r_prev(:)            ! 직전 front 잔차 (Aitken 완화용)
     real(dp) :: omega = 0.5_dp                    ! front 완화계수 (Aitken 적응)
     real(dp) :: g_vap = 0.0_dp
     real(dp) :: t_vap = 0.0_dp          ! 물성 평가 온도 (증발부 증기)
     real(dp) :: q_transport = 0.0_dp    ! 현재 수송 열량 [W]
     real(dp) :: q_max = 0.0_dp          ! 수송 한계 (최소값) [W]
     real(dp) :: q_lim(5) = 0.0_dp       ! 개별 한계 [W]
     integer  :: limit_id = LIM_NONE     ! 지배 한계
     real(dp) :: L_active = 0.0_dp       ! 활성 응축 길이 [m]
     real(dp) :: L_eff = 0.0_dp          ! 유효 길이 [m]
     real(dp) :: dp_cap_max = 0.0_dp, dp_liq = 0.0_dp, dp_vap = 0.0_dp, dp_grav = 0.0_dp
     real(dp) :: re_vap = 0.0_dp, mach_vap = 0.0_dp
     real(dp) :: scale = 1.0_dp          ! 수송 한계 축소 계수 (증발부 전도체에 적용)
     logical  :: limited = .false.       ! 수송 한계 도달 (dry-out)
     logical  :: shut_off = .false.      ! VCHP 응축부 전체 차단
     logical  :: out_of_range = .false.  ! 물성 테이블 범위 이탈
     logical  :: initialized = .false.
  end type hp_state

  public :: hp_init, hp_update, hp_check_config
  public :: groove_keff, conductance_per_length, liquid_dp, vapor_dp, &
            capillary_limit, other_limits, ncg_front, vchp_charge, shah_fre, vapor_flow

contains

  !=============================================================================
  ! 설정 점검. 문제가 있으면 ierr /= 0, msg 에 사유
  !=============================================================================
  subroutine hp_check_config(cfg, ierr, msg)
    type(hp_config), intent(in) :: cfg
    integer, intent(out) :: ierr
    character(len=*), intent(out) :: msg
    real(dp) :: w_fin
    ierr = 0; msg = 'ok'
    w_fin = 2.0_dp * PI * cfg%r_vap / cfg%n_groove - cfg%groove_w
    if (cfg%L_evap <= 0 .or. cfg%L_cond <= 0 .or. cfg%L_adia < 0) then
       ierr = 1; msg = 'lengths must be positive'
    else if (w_fin <= 0) then
       ierr = 2; msg = 'groove too wide: fin width <= 0 (check n_groove, groove_w, r_vap)'
    else if (cfg%r_out <= cfg%r_vap + cfg%groove_d) then
       ierr = 3; msg = 'r_out must exceed r_vap + groove_d'
    else if (cfg%n_gas > 0 .and. cfg%v_res <= 0) then
       ierr = 4; msg = 'VCHP needs v_res > 0'
    else if (cfg%relax <= 0 .or. cfg%relax > 1) then
       ierr = 5; msg = 'relax must be in (0,1]'
    end if
  end subroutine hp_check_config

  !=============================================================================
  ! 초기화: 구간 길이(솔버 노드 분할)를 등록
  !=============================================================================
  subroutine hp_init(cfg, dx_e, dx_c, st, ierr, msg)
    type(hp_config), intent(in) :: cfg
    real(dp), intent(in) :: dx_e(:), dx_c(:)
    type(hp_state), intent(inout) :: st
    integer, intent(out) :: ierr
    character(len=*), intent(out) :: msg
    call hp_check_config(cfg, ierr, msg)
    if (ierr /= 0) return
    if (abs(sum(dx_e) - cfg%L_evap) > 1.0e-6_dp * cfg%L_evap .or. &
        abs(sum(dx_c) - cfg%L_cond) > 1.0e-6_dp * cfg%L_cond) then
       ierr = 10; msg = 'sum of segment lengths must equal L_evap / L_cond'
       return
    end if
    st%ne = size(dx_e); st%nc = size(dx_c)
    if (allocated(st%dx_e)) deallocate(st%dx_e, st%dx_c, st%g_evap, st%g_cond, st%f_act, st%r_prev)
    allocate(st%dx_e(st%ne), st%dx_c(st%nc), st%g_evap(st%ne), st%g_cond(st%nc), &
             st%f_act(st%nc), st%r_prev(st%nc))
    st%dx_e = dx_e; st%dx_c = dx_c
    st%g_evap = 0; st%g_cond = 0; st%f_act = 1; st%g_vap = 0
    st%r_prev = 0; st%omega = cfg%relax
    st%scale = 1; st%initialized = .false.
  end subroutine hp_init

  !=============================================================================
  ! 매 반복(또는 시간 스텝) 호출: 현재 노드 온도로 전도체 값 갱신
  !
  !   T_e(ne)  : 증발부 벽 노드 온도      T_ve : 증발부 증기 노드 온도
  !   T_c(nc)  : 응축부 벽 노드 온도      T_vc : 응축부 증기 노드 온도
  !              (j=1 이 단열부 쪽, j=nc 가 저장조 쪽)
  !   T_res    : NCG 저장조 온도 (CCHP 이면 아무 값)
  !=============================================================================
  subroutine hp_update(cfg, st, T_e, T_ve, T_vc, T_c, T_res)
    type(hp_config), intent(in) :: cfg
    type(hp_state), intent(inout) :: st
    real(dp), intent(in) :: T_e(:), T_ve, T_vc, T_c(:), T_res

    type(fluid_props) :: p
    real(dp) :: gpl, q_cur, target_scale, dTv, g_vap_new, denom
    real(dp) :: ge_new(st%ne), f_new(st%nc), r(st%nc), dr(st%nc)
    logical :: shut

    p = get_props(T_ve)
    st%t_vap = T_ve
    st%out_of_range = p%out_of_range

    ! 1) 현재 수송 열량 (직전 전도체 기준). 첫 호출은 비축소 전도체로 추정
    gpl = conductance_per_length(cfg, p)
    if (st%initialized) then
       q_cur = sum(st%g_evap * (T_e - T_ve))
    else
       q_cur = sum(cfg%mult_evap * gpl * st%dx_e * (T_e - T_ve))
    end if
    q_cur = abs(q_cur)
    st%q_transport = q_cur

    ! 2) VCHP : NCG flat-front -> 응축부 구간별 활성 비율
    if (cfg%n_gas > 0.0_dp) then
       call ncg_front(cfg%n_gas, cfg%v_res, PI * cfg%r_vap**2, T_vc, T_res, T_c, st%dx_c, &
                      f_new, shut)
    else
       f_new = 1.0_dp; shut = .false.
    end if
    st%shut_off = shut

    ! front 위치는 노드 온도에 매우 민감(강성)하므로 Aitken(Irons-Tuck) 적응 완화로 갱신.
    ! 고정 완화로는 front 가 구간 사이를 왕복하며 수렴하지 않는 경우가 많다.
    if (st%initialized) then
       r = f_new - st%f_act
       dr = r - st%r_prev
       denom = dot_product(dr, dr)
       if (denom > 1.0e-24_dp) then
          ! 증가는 직전 값의 2배까지만 허용 (구간 경계의 꺾임에서 생기는 limit cycle 방지)
          st%omega = min(-st%omega * dot_product(st%r_prev, dr) / denom, 2.0_dp * st%omega)
          st%omega = min(max(st%omega, 0.02_dp), cfg%relax)
       end if
       st%f_act = min(max(st%f_act + st%omega * r, 0.0_dp), 1.0_dp)
       st%r_prev = r
    else
       st%f_act = f_new
    end if
    st%L_active = sum(st%f_act * st%dx_c)
    st%L_eff = cfg%L_adia + 0.5_dp * (cfg%L_evap + st%L_active)

    ! 3) 수송 한계
    st%q_lim(LIM_CAPILLARY) = capillary_limit(cfg, p, st%L_eff)
    call other_limits(cfg, p, st%L_eff, st%q_lim(LIM_SONIC), st%q_lim(LIM_ENTRAIN), &
                      st%q_lim(LIM_BOILING), st%q_lim(LIM_VISCOUS))
    if (.not. cfg%use_boiling_formula) st%q_lim(LIM_BOILING) = huge(1.0_dp)
    if (cfg%q_flux_max > 0.0_dp) st%q_lim(LIM_BOILING) = min(st%q_lim(LIM_BOILING), &
         cfg%q_flux_max * 2.0_dp * PI * cfg%r_out * cfg%L_evap)
    st%limit_id = minloc(st%q_lim, dim=1)
    st%q_max = st%q_lim(st%limit_id)

    ! 4) 진단용 압력 수지 (현재 Q 기준)
    st%dp_cap_max = 2.0_dp * p%sigma * cos(cfg%contact_angle) / cfg%groove_w
    st%dp_liq  = liquid_dp(cfg, p, q_cur, st%L_eff)
    st%dp_vap  = vapor_dp(cfg, p, q_cur, st%L_eff)
    call vapor_flow(cfg, p, q_cur, dTv, st%re_vap, st%mach_vap)
    st%dp_grav = p%rho_l * cfg%g_accel * cfg%elevation

    ! 5) 수송 한계 초과 시 전도체 축소 계수 (고정점: Q = Q_max)
    target_scale = 1.0_dp
    if (cfg%limit_transport) then
       if (st%q_max <= 0.0_dp) then
          target_scale = 1.0e-6_dp
       else if (q_cur > 0.0_dp) then
          target_scale = min(1.0_dp, max(1.0e-6_dp, st%scale * st%q_max / q_cur))
       end if
    end if
    st%scale = target_scale
    st%limited = (st%scale < 0.999_dp)

    ! 6) 전도체
    ge_new = cfg%mult_evap * gpl * st%dx_e * st%scale

    ! 증기 전도체: G_v = Q / dT_v,  dT_v = T dP_v / (rho_v h_fg)
    if (q_cur > 1.0e-9_dp) then
       dTv = T_ve * st%dp_vap / (p%rho_v * p%h_fg)
       g_vap_new = q_cur / max(dTv, 1.0e-30_dp)
    else   ! Q -> 0 : 층류 기울기
       dTv = T_ve * vapor_dp(cfg, p, 1.0_dp, st%L_eff) / (p%rho_v * p%h_fg)
       g_vap_new = 1.0_dp / max(dTv, 1.0e-30_dp)
    end if
    g_vap_new = min(g_vap_new, cfg%g_vap_max)

    ! 7) 완화 갱신 (응축부는 이미 완화된 f_act 를 사용)
    if (st%initialized) then
       st%g_evap = st%g_evap + cfg%relax * (ge_new - st%g_evap)
       st%g_vap  = st%g_vap  + cfg%relax * (g_vap_new - st%g_vap)
    else
       st%g_evap = ge_new; st%g_vap = g_vap_new
       st%initialized = .true.
    end if
    st%g_cond = cfg%mult_cond * gpl * st%dx_c &
                * (st%f_act + (1.0_dp - st%f_act) * cfg%blocked_frac)
  end subroutine hp_update

  !=============================================================================
  ! 물리 모델 세부 함수 (단위 테스트/참조 비교용으로 공개)
  !=============================================================================

  !> 사각 축방향 그루브의 유효 반경방향 열전도도 (Chi, 1976)
  pure real(dp) function groove_keff(cfg, k_l) result(k_eff)
    type(hp_config), intent(in) :: cfg
    real(dp), intent(in) :: k_l
    real(dp) :: w, wf, d, kw
    w = cfg%groove_w; d = cfg%groove_d; kw = cfg%k_wall
    wf = 2.0_dp * PI * cfg%r_vap / cfg%n_groove - w
    k_eff = (wf * k_l * kw * d + w * k_l * (0.185_dp * wf * kw + d * k_l)) &
          / ((w + wf) * (0.185_dp * wf * kw + d * k_l))
  end function groove_keff

  !> 단위 길이당 반경방향 전도도 [W/m.K] : 벽 + 그루브 액체층 (직렬)
  pure real(dp) function conductance_per_length(cfg, p) result(gpl)
    type(hp_config), intent(in) :: cfg
    type(fluid_props), intent(in) :: p
    real(dp) :: r_g
    r_g = cfg%r_vap + cfg%groove_d
    gpl = 2.0_dp * PI / ( log(cfg%r_out / r_g) / cfg%k_wall &
                        + log(r_g / cfg%r_vap) / groove_keff(cfg, p%k_l) )
  end function conductance_per_length

  !> 사각 덕트 층류 Fanning fRe (Shah & London). alpha = 짧은변/긴변
  pure real(dp) function shah_fre(alpha)
    real(dp), intent(in) :: alpha
    shah_fre = 24.0_dp * (1.0_dp - 1.3553_dp*alpha + 1.9467_dp*alpha**2 &
               - 1.7012_dp*alpha**3 + 0.9564_dp*alpha**4 - 0.2537_dp*alpha**5)
  end function shah_fre

  !> 그루브 액체 귀환 압력강하 [Pa]
  !  자유표면(전단 0)을 대칭면으로 보고 폭 w, 높이 2d 의 사각 덕트로 등가화
  pure real(dp) function liquid_dp(cfg, p, Q, L_eff) result(dpl)
    type(hp_config), intent(in) :: cfg
    type(fluid_props), intent(in) :: p
    real(dp), intent(in) :: Q, L_eff
    real(dp) :: w, d, A_l, D_h, alpha
    w = cfg%groove_w; d = cfg%groove_d
    A_l = cfg%n_groove * w * d
    D_h = 4.0_dp * w * d / (w + 2.0_dp * d)
    alpha = min(w, 2.0_dp*d) / max(w, 2.0_dp*d)
    dpl = 2.0_dp * shah_fre(alpha) * p%mu_l * L_eff * (Q / p%h_fg) / (p%rho_l * A_l * D_h**2)
  end function liquid_dp

  !> 증기 유동 압력강하 [Pa] (층류 Hagen-Poiseuille / 난류 Blasius, 비압축)
  pure real(dp) function vapor_dp(cfg, p, Q, L_eff) result(dpv)
    type(hp_config), intent(in) :: cfg
    type(fluid_props), intent(in) :: p
    real(dp), intent(in) :: Q, L_eff
    real(dp) :: A_v, D, u, re, f
    A_v = PI * cfg%r_vap**2; D = 2.0_dp * cfg%r_vap
    call vapor_flow(cfg, p, Q, u, re)
    if (re < RE_TRANS) then
       dpv = 32.0_dp * p%mu_v * u * L_eff / D**2
    else
       f = 0.079_dp * re**(-0.25_dp)
       dpv = 2.0_dp * f * p%rho_v * u**2 * L_eff / D
    end if
  end function vapor_dp

  !> 증기 평균속도 [m/s], 레이놀즈수, 마하수
  pure subroutine vapor_flow(cfg, p, Q, u, re, mach)
    type(hp_config), intent(in) :: cfg
    type(fluid_props), intent(in) :: p
    real(dp), intent(in) :: Q
    real(dp), intent(out) :: u, re
    real(dp), intent(out), optional :: mach
    real(dp) :: A_v, mdot
    A_v = PI * cfg%r_vap**2
    mdot = Q / p%h_fg
    u = mdot / (p%rho_v * A_v)
    re = mdot * 2.0_dp * cfg%r_vap / (A_v * p%mu_v)
    if (present(mach)) mach = u / sqrt(p%gamma_v * R_UNIV / M_NH3 * p%T)
  end subroutine vapor_flow

  !> 모세관 한계 [W] : 2 sigma cos(theta)/w = dP_l(Q) + dP_v(Q) + dP_g  를 이분법으로 풂
  pure real(dp) function capillary_limit(cfg, p, L_eff) result(qcap)
    type(hp_config), intent(in) :: cfg
    type(fluid_props), intent(in) :: p
    real(dp), intent(in) :: L_eff
    real(dp) :: dpc, dpg, qlo, qhi, qm
    integer :: it
    dpc = 2.0_dp * p%sigma * cos(cfg%contact_angle) / cfg%groove_w
    dpg = p%rho_l * cfg%g_accel * cfg%elevation
    qcap = 0.0_dp
    if (dpg >= dpc) return
    qlo = 0.0_dp; qhi = 1.0_dp
    do it = 1, 60
       if (balance(qhi) > 0.0_dp) exit
       qlo = qhi; qhi = 2.0_dp * qhi
    end do
    do it = 1, 80
       qm = 0.5_dp * (qlo + qhi)
       if (balance(qm) > 0.0_dp) then
          qhi = qm
       else
          qlo = qm
       end if
       if (qhi - qlo < 1.0e-9_dp * qhi) exit
    end do
    qcap = 0.5_dp * (qlo + qhi)
  contains
    pure real(dp) function balance(q)
      real(dp), intent(in) :: q
      balance = liquid_dp(cfg, p, q, L_eff) + vapor_dp(cfg, p, q, L_eff) + dpg - dpc
    end function balance
  end function capillary_limit

  !> 음속 / 비말동반 / 비등 / 점성 한계 [W]
  pure subroutine other_limits(cfg, p, L_eff, q_sonic, q_entr, q_boil, q_visc)
    type(hp_config), intent(in) :: cfg
    type(fluid_props), intent(in) :: p
    real(dp), intent(in) :: L_eff
    real(dp), intent(out) :: q_sonic, q_entr, q_boil, q_visc
    real(dp) :: A_v, Rg, r_g, dpc
    A_v = PI * cfg%r_vap**2
    Rg = R_UNIV / M_NH3
    r_g = cfg%r_vap + cfg%groove_d
    dpc = 2.0_dp * p%sigma * cos(cfg%contact_angle) / cfg%groove_w
    ! 음속 (Levy/Busse)
    q_sonic = A_v * p%rho_v * p%h_fg * sqrt(p%gamma_v * Rg * p%T / (2.0_dp * (p%gamma_v + 1.0_dp)))
    ! 비말동반 : 액체 표면 수력반경 = w/2
    q_entr = A_v * p%h_fg * sqrt(p%sigma * p%rho_v / (2.0_dp * 0.5_dp * cfg%groove_w))
    ! 비등 (증발부 그루브 액체층)
    q_boil = 2.0_dp * PI * cfg%L_evap * groove_keff(cfg, p%k_l) * p%T &
             / (p%h_fg * p%rho_v * log(r_g / cfg%r_vap)) &
             * max(2.0_dp * p%sigma / cfg%r_nuc - dpc, 0.0_dp)
    ! 점성 (Busse) : 저온 시동 시 지배
    q_visc = A_v * cfg%r_vap**2 * p%h_fg * p%rho_v * p%psat / (16.0_dp * p%mu_v * L_eff)
  end subroutine other_limits

  !=============================================================================
  ! VCHP NCG flat-front (구간별 온도와 증기 분압 고려)
  !
  !  전체 압력 = P_sat(T_vc). NCG 구간의 증기 분압 = 벽 온도 포화압 P_sat(T_wall).
  !  NCG 몰수: 저장조 + 저장조 쪽 끝에서부터 차단된 응축부 구간의 합이 n_gas 가 되도록
  !  front 위치를 찾는다.  n_k = (P_sat(T_vc) - P_sat(T_k)) * V_k / (R T_k)
  !=============================================================================
  pure subroutine ncg_front(n_gas, v_res, a_vap, T_vc, T_res, T_c, dx_c, f, shut)
    real(dp), intent(in) :: n_gas, v_res, a_vap, T_vc, T_res, T_c(:), dx_c(:)
    real(dp), intent(out) :: f(size(T_c))
    logical, intent(out) :: shut
    real(dp) :: pv, remain, nk
    integer :: j
    pv = psat_of_T(T_vc)
    f = 1.0_dp; shut = .false.
    remain = n_gas - gas_moles(pv, T_res, v_res)
    if (remain <= 0.0_dp) return
    do j = size(T_c), 1, -1
       nk = gas_moles(pv, T_c(j), a_vap * dx_c(j))
       if (remain >= nk) then
          f(j) = 0.0_dp
          remain = remain - nk
       else
          f(j) = 1.0_dp - remain / nk
          remain = 0.0_dp
          exit
       end if
    end do
    shut = (remain > 0.0_dp)
  contains
    pure real(dp) function gas_moles(pv_, T_, V_)
      real(dp), intent(in) :: pv_, T_, V_
      real(dp) :: pg
      pg = max(pv_ - psat_of_T(min(T_, T_vc)), 1.0e-6_dp * pv_)
      gas_moles = pg * V_ / (R_UNIV * T_)
    end function gas_moles
  end subroutine ncg_front

  !> VCHP NCG 충전량 설계 [mol]
  !  증기온도가 T_full 이고 저장조가 T_res 일 때 front 가 응축부 끝(=완전 개방)에
  !  오도록 하는 충전량.  (P_sat(T_full) - P_sat(T_res)) V_res / (R T_res)
  pure real(dp) function vchp_charge(T_full, T_res, v_res)
    real(dp), intent(in) :: T_full, T_res, v_res
    vchp_charge = (psat_of_T(T_full) - psat_of_T(min(T_res, T_full))) * v_res / (R_UNIV * T_res)
  end function vchp_charge

end module hp_model
