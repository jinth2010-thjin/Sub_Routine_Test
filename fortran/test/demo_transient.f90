!===============================================================================
! demo_transient : 위성 장비 - 히트파이프 - 라디에이터 과도 열해석 데모
!
!  Space Systems Thermal 이 하는 일(노드 네트워크 + 암시적 시간적분)을 아주 작게
!  흉내 낸 독립 솔버입니다. 실제 솔버에서는 hp_update 를 USER1 서브루틴에서
!  호출하게 됩니다. 여기서는 hp_model 이 노드 네트워크 안에서 수렴하는지,
!  CCHP/VCHP 거동이 물리적으로 타당한지 확인하는 용도입니다.
!
!  노드 구성
!    1        장비 박스 (발열 Q_box)
!    2..3     증발부 벽 (2구간)
!    4, 5     증기 노드 (증발부 / 응축부)
!    6..15    응축부 벽 (10구간, 6 이 단열부 쪽)
!    16..25   라디에이터 (10구간, 우주로 복사 + 궤도 흡수열)
!    26       VCHP NCG 저장조 (라디에이터 끝단에 부착, 냉저장조)
!
!  케이스 (명령행 인자):
!    cchp_hot | vchp_hot   : 고온 케이스 (Q_box 60 W, 흡수 플럭스 250 W/m2)
!    cchp_cold | vchp_cold : 저온 케이스 (Q_box 30 W, 흡수 플럭스 150 W/m2)
!    cchp_pulse            : 고온 + 매 궤도 15분간 400 W 펄스 -> 모세관 한계 초과(dry-out) 후 회복
!    cchp_tilt             : 고온 + 지상 열진공시험 역경사 10 mm (중력 수두가 모세관 압력을 잠식)
!  출력: build/demo_<case>.csv
!===============================================================================
program demo_transient
  use hp_props
  use hp_model
  implicit none

  integer, parameter :: NE = 2, NC = 10, NN = 26
  integer, parameter :: IBOX = 1, IE0 = 2, IVE = 4, IVC = 5, IC0 = 6, IR0 = 16, IRES = 26
  real(dp), parameter :: SB = 5.670374419e-8_dp, PI = 3.14159265358979_dp
  real(dp), parameter :: ORBIT = 5700.0_dp, SUNLIT = 3600.0_dp   ! 95분 궤도, 60분 일조
  real(dp), parameter :: DT = 30.0_dp
  integer, parameter :: N_ORBIT = 12

  type(hp_config) :: cfg
  type(hp_state) :: st
  character(len=32) :: casename
  real(dp) :: C(NN), Gfix(NN,NN), area_rad(NN), eps(NN), T(NN), Told(NN), Tit(NN)
  real(dp) :: A(NN,NN), b(NN), q_box, q_env_peak, tnow, q_abs, a_wall, g_ax, q_now
  real(dp) :: dx_e(NE), dx_c(NC)
  integer :: i, j, n, it, u, ierr, nsteps, max_it_used, n_limited, n_nonconv
  character(len=200) :: msg
  logical :: is_vchp

  casename = 'vchp_cold'
  if (command_argument_count() >= 1) call get_command_argument(1, casename)

  ! ---------------- 케이스 정의
  is_vchp = (casename(1:4) == 'vchp')
  select case (trim(casename))
  case ('cchp_hot', 'vchp_hot');   q_box = 60.0_dp;  q_env_peak = 250.0_dp
  case ('cchp_cold', 'vchp_cold'); q_box = 30.0_dp;  q_env_peak = 150.0_dp
  case ('cchp_pulse');             q_box = 60.0_dp;  q_env_peak = 250.0_dp
  case ('cchp_tilt');              q_box = 60.0_dp;  q_env_peak = 250.0_dp
     cfg%elevation = 0.010_dp
  case default
     print *, 'unknown case: ', trim(casename); stop 1
  end select

  ! ---------------- 히트파이프 설정 (기본값 = 12.7 mm 급 Al/암모니아 그루브 HP)
  dx_e = cfg%L_evap / NE
  dx_c = cfg%L_cond / NC
  if (is_vchp) then
     cfg%v_res = 30.0e-6_dp                                   ! 30 cm3 저장조
     cfg%n_gas = vchp_charge(293.15_dp, 263.15_dp, cfg%v_res) ! T_v=20C, T_res=-10C 에서 완전 개방
  end if
  call hp_init(cfg, dx_e, dx_c, st, ierr, msg)
  if (ierr /= 0) then
     print *, 'hp_init error: ', trim(msg); stop 1
  end if

  ! ---------------- 노드 네트워크
  C = 0; Gfix = 0; area_rad = 0; eps = 0
  C(IBOX) = 6000.0_dp
  C(IE0:IE0+NE-1) = 25.0_dp
  C(IVE) = 1.0_dp; C(IVC) = 1.0_dp
  C(IC0:IC0+NC-1) = 15.0_dp
  C(IR0:IR0+NC-1) = 400.0_dp
  C(IRES) = 40.0_dp
  ! 박스 - 증발부 (접촉 + 새들)
  do i = 0, NE-1
     call link(IBOX, IE0+i, 10.0_dp)
  end do
  ! 박스 - 라디에이터 패널 직접 경로 (장착 구조물, 히트파이프와 병렬)
  call link(IBOX, IR0, 0.5_dp)
  ! 응축부 벽 - 라디에이터 (새들 + 패널)
  do i = 0, NC-1
     call link(IC0+i, IR0+i, 3.0_dp)
  end do
  ! 라디에이터 패널 측면 전도
  do i = 0, NC-2
     call link(IR0+i, IR0+i+1, 0.5_dp)
  end do
  ! 히트파이프 벽 축방향 전도 (작동유체와 무관한 순수 금속 전도 경로)
  a_wall = PI * (cfg%r_out**2 - (cfg%r_vap + cfg%groove_d)**2)
  call link(IE0, IE0+1, cfg%k_wall * a_wall / dx_e(1))
  do i = 0, NC-2
     call link(IC0+i, IC0+i+1, cfg%k_wall * a_wall / dx_c(1))
  end do
  g_ax = cfg%k_wall * a_wall / (cfg%L_adia + 0.5_dp*(dx_e(1) + dx_c(1)))
  call link(IE0+NE-1, IC0, g_ax)
  ! 저장조 - 라디에이터 끝단
  call link(IRES, IR0+NC-1, 0.3_dp)
  ! 복사면
  area_rad(IR0:IR0+NC-1) = 0.06_dp; eps(IR0:IR0+NC-1) = 0.85_dp   ! 총 0.6 m2
  area_rad(IRES) = 0.004_dp;        eps(IRES) = 0.85_dp

  T = 273.15_dp
  nsteps = nint(N_ORBIT * ORBIT / DT)
  max_it_used = 0; n_limited = 0; n_nonconv = 0

  open(newunit=u, file='build/demo_'//trim(casename)//'.csv', status='replace')
  write(u,'(a)') 'time_min,T_box,T_ve,T_vc,T_res,T_rad_mean,T_rad_end,q_box,Q_hp,Q_max,' // &
                 'L_active,scale,limit,iters'

  ! ---------------- 시간 적분 (후진 오일러 + Picard 반복)
  do n = 1, nsteps
     tnow = n * DT
     q_abs = merge(q_env_peak, 0.0_dp, mod(tnow, ORBIT) < SUNLIT)
     q_now = q_box
     if (trim(casename) == 'cchp_pulse' .and. mod(tnow, ORBIT) >= 1200.0_dp &
         .and. mod(tnow, ORBIT) < 2100.0_dp) q_now = 400.0_dp
     Told = T
     do it = 1, 200
        Tit = T
        ! >>> 실제 Space Systems Thermal 에서는 이 부분이 USER1 서브루틴 호출 <<<
        call hp_update(cfg, st, T(IE0:IE0+NE-1), T(IVE), T(IVC), T(IC0:IC0+NC-1), T(IRES))
        call assemble()
        call solve(A, b, T)
        if (maxval(abs(T - Tit)) < 1.0e-4_dp .and. it > 2) exit
     end do
     if (it > 200) n_nonconv = n_nonconv + 1
     max_it_used = max(max_it_used, it)
     if (st%limited) n_limited = n_limited + 1
     if (mod(n, 4) == 0) then
        write(u,'(f10.2,9(",",f12.5),",",f10.5,",",f8.5,",",a,",",i4)') tnow/60.0_dp, &
             T(IBOX)-273.15_dp, T(IVE)-273.15_dp, T(IVC)-273.15_dp, T(IRES)-273.15_dp, &
             sum(T(IR0:IR0+NC-1))/NC-273.15_dp, T(IR0+NC-1)-273.15_dp, q_now, &
             st%q_transport, st%q_max, st%L_active, st%scale, trim(LIMIT_NAME(st%limit_id)), it
     end if
  end do
  close(u)

  print '(a,a12,a,f7.2,a,f7.2,a,i4,a,i5,a,i5)', 'case=', trim(casename), &
       '  last-orbit T_box min/max [C] =', minval_last(), ' /', maxval_last(), &
       '  max iters=', max_it_used, '  limited steps=', n_limited, '  non-converged steps=', n_nonconv

contains

  subroutine link(i1, i2, g)
    integer, intent(in) :: i1, i2
    real(dp), intent(in) :: g
    Gfix(i1,i2) = Gfix(i1,i2) + g
    Gfix(i2,i1) = Gfix(i2,i1) + g
  end subroutine link

  subroutine assemble()
    real(dp) :: G(NN,NN), hr
    integer :: k
    G = Gfix
    ! 히트파이프 전도체 (hp_model 이 매 반복 갱신)
    do k = 1, NE
       G(IE0+k-1, IVE) = G(IE0+k-1, IVE) + st%g_evap(k)
       G(IVE, IE0+k-1) = G(IVE, IE0+k-1) + st%g_evap(k)
    end do
    do k = 1, NC
       G(IC0+k-1, IVC) = G(IC0+k-1, IVC) + st%g_cond(k)
       G(IVC, IC0+k-1) = G(IVC, IC0+k-1) + st%g_cond(k)
    end do
    G(IVE, IVC) = G(IVE, IVC) + st%g_vap
    G(IVC, IVE) = G(IVC, IVE) + st%g_vap

    A = 0; b = 0
    do k = 1, NN
       A(k,k) = C(k) / DT
       b(k) = C(k) / DT * Told(k)
       do j = 1, NN
          if (j /= k .and. G(k,j) > 0) then
             A(k,k) = A(k,k) + G(k,j)
             A(k,j) = A(k,j) - G(k,j)
          end if
       end do
       if (area_rad(k) > 0) then   ! 복사 선형화: T^4 ~ 4 T*^3 T - 3 T*^4
          hr = eps(k) * SB * area_rad(k)
          A(k,k) = A(k,k) + 4.0_dp * hr * T(k)**3
          b(k) = b(k) + 3.0_dp * hr * T(k)**4
          if (k >= IR0 .and. k < IR0+NC) b(k) = b(k) + q_abs * area_rad(k)
       end if
    end do
    b(IBOX) = b(IBOX) + q_now
  end subroutine assemble

  subroutine solve(M, r, x)   ! 부분 피벗 가우스 소거
    real(dp), intent(in) :: M(NN,NN), r(NN)
    real(dp), intent(out) :: x(NN)
    real(dp) :: W(NN,NN+1), tmp(NN+1), fct
    integer :: k, p, ii
    W(:,1:NN) = M; W(:,NN+1) = r
    do k = 1, NN-1
       p = maxloc(abs(W(k:NN,k)), dim=1) + k - 1
       if (p /= k) then
          tmp = W(k,:); W(k,:) = W(p,:); W(p,:) = tmp
       end if
       do ii = k+1, NN
          fct = W(ii,k) / W(k,k)
          W(ii,k:NN+1) = W(ii,k:NN+1) - fct * W(k,k:NN+1)
       end do
    end do
    do k = NN, 1, -1
       x(k) = (W(k,NN+1) - dot_product(W(k,k+1:NN), x(k+1:NN))) / W(k,k)
    end do
  end subroutine solve

  ! 마지막 궤도의 박스 온도 최소/최대 (CSV 를 다시 읽지 않도록 간단히 재계산 대신 파일에서)
  real(dp) function minval_last()
    minval_last = last_orbit_stat(.true.)
  end function minval_last
  real(dp) function maxval_last()
    maxval_last = last_orbit_stat(.false.)
  end function maxval_last
  real(dp) function last_orbit_stat(want_min)
    logical, intent(in) :: want_min
    integer :: uu, ios
    character(len=400) :: line
    real(dp) :: tm, tb
    last_orbit_stat = merge(1.0e9_dp, -1.0e9_dp, want_min)
    open(newunit=uu, file='build/demo_'//trim(casename)//'.csv', status='old')
    read(uu,'(a)') line
    do
       read(uu,'(a)', iostat=ios) line
       if (ios /= 0) exit
       read(line,*) tm, tb
       if (tm >= (N_ORBIT-1) * ORBIT / 60.0_dp) then
          if (want_min) then
             last_orbit_stat = min(last_orbit_stat, tb)
          else
             last_orbit_stat = max(last_orbit_stat, tb)
          end if
       end if
    end do
    close(uu)
  end function last_orbit_stat

end program demo_transient
