!===============================================================================
! hp_driver : 여러 개의 히트파이프를 입력 파일(namelist)로 정의하고, 솔버의
!             노드 온도를 읽어 전도체 값을 써 넣는 "솔버 독립" 드라이버
!
!  솔버와의 연결은 두 개의 프로시저로만 이루어집니다.
!     get_temp(node_id)          -> 노드 온도 [K]
!     set_cond(cond_id, value)   -> 전도체 값 [W/K] 설정
!  Space Systems Thermal(TMG) 에서는 USER1 서브루틴 안에서 TMG API 로 이 두
!  동작을 구현해 넘겨주면 됩니다 (fortran/src/user1_template.f90 참고).
!
!  입력 파일 예: examples/hp_input.nml
!===============================================================================
module hp_driver
  use hp_props
  use hp_model
  implicit none
  private

  integer, parameter, public :: MAX_HP = 50, MAX_SEG = 100

  abstract interface
     function get_temp_if(node_id) result(T)
       import :: dp
       integer, intent(in) :: node_id
       real(dp) :: T
     end function get_temp_if
     subroutine set_cond_if(cond_id, value)
       import :: dp
       integer, intent(in) :: cond_id
       real(dp), intent(in) :: value
     end subroutine set_cond_if
  end interface
  public :: get_temp_if, set_cond_if

  type :: hp_item
     character(len=32) :: name = ''
     type(hp_config) :: cfg
     type(hp_state)  :: st
     integer :: ne = 0, nc = 0
     integer :: evap_nodes(MAX_SEG) = 0, cond_nodes(MAX_SEG) = 0
     integer :: evap_conds(MAX_SEG) = 0, cond_conds(MAX_SEG) = 0
     integer :: vap_node_e = 0, vap_node_c = 0, res_node = 0, vap_cond = 0
  end type hp_item

  type(hp_item), allocatable, save :: hps(:)
  integer, save :: n_hp = 0
  logical, save :: ready = .false.

  public :: hp_driver_init, hp_driver_update, hp_driver_report, hp_driver_ready, &
            hp_driver_count, hp_driver_state

contains

  logical function hp_driver_ready()
    hp_driver_ready = ready
  end function hp_driver_ready

  integer function hp_driver_count()
    hp_driver_count = n_hp
  end function hp_driver_count

  !> 진단/테스트용: k 번째 히트파이프 상태 복사본
  function hp_driver_state(k) result(st)
    integer, intent(in) :: k
    type(hp_state) :: st
    st = hps(k)%st
  end function hp_driver_state

  !=============================================================================
  ! 입력 파일 읽기
  !=============================================================================
  subroutine hp_driver_init(filename, ierr, msg)
    character(len=*), intent(in) :: filename
    integer, intent(out) :: ierr
    character(len=*), intent(out) :: msg

    ! --- namelist 변수 (hp_config 기본값을 그대로 사용)
    type(hp_config) :: d
    character(len=32) :: name
    integer :: ne, nc, evap_nodes(MAX_SEG), cond_nodes(MAX_SEG), evap_conds(MAX_SEG), &
               cond_conds(MAX_SEG), vap_node_e, vap_node_c, res_node, vap_cond
    real(dp) :: dx_e(MAX_SEG), dx_c(MAX_SEG)
    real(dp) :: L_evap, L_adia, L_cond, r_vap, groove_w, groove_d, r_out, k_wall, &
                contact_angle, elevation, g_accel, r_nuc, mult_evap, mult_cond, &
                n_gas, v_res, blocked_frac, relax, g_vap_max, q_flux_max, &
                charge_T_full, charge_T_res
    integer :: n_groove
    logical :: limit_transport, use_boiling_formula
    integer :: nhp, u, k, ios
    character(len=256) :: iomsg_

    namelist /hp_global/ nhp
    namelist /heat_pipe/ name, ne, nc, dx_e, dx_c, evap_nodes, cond_nodes, evap_conds, &
         cond_conds, vap_node_e, vap_node_c, res_node, vap_cond, &
         L_evap, L_adia, L_cond, r_vap, groove_w, groove_d, n_groove, r_out, k_wall, &
         contact_angle, elevation, g_accel, r_nuc, mult_evap, mult_cond, &
         n_gas, v_res, blocked_frac, relax, g_vap_max, limit_transport, &
         use_boiling_formula, q_flux_max, charge_T_full, charge_T_res

    ierr = 0; msg = 'ok'
    open(newunit=u, file=filename, status='old', action='read', iostat=ios, iomsg=iomsg_)
    if (ios /= 0) then
       ierr = 1; msg = 'cannot open '//trim(filename)//': '//trim(iomsg_); return
    end if
    nhp = 0
    read(u, nml=hp_global, iostat=ios, iomsg=iomsg_)
    if (ios /= 0 .or. nhp < 1 .or. nhp > MAX_HP) then
       ierr = 2; msg = 'bad &hp_global: '//trim(iomsg_); close(u); return
    end if
    if (allocated(hps)) deallocate(hps)
    allocate(hps(nhp))
    n_hp = nhp

    do k = 1, nhp
       ! 기본값으로 리셋
       name = ''; ne = 0; nc = 0; dx_e = 0; dx_c = 0
       evap_nodes = 0; cond_nodes = 0; evap_conds = 0; cond_conds = 0
       vap_node_e = 0; vap_node_c = 0; res_node = 0; vap_cond = 0
       L_evap = d%L_evap; L_adia = d%L_adia; L_cond = d%L_cond
       r_vap = d%r_vap; groove_w = d%groove_w; groove_d = d%groove_d
       n_groove = d%n_groove; r_out = d%r_out; k_wall = d%k_wall
       contact_angle = d%contact_angle; elevation = d%elevation; g_accel = d%g_accel
       r_nuc = d%r_nuc; mult_evap = d%mult_evap; mult_cond = d%mult_cond
       n_gas = d%n_gas; v_res = d%v_res; blocked_frac = d%blocked_frac
       relax = d%relax; g_vap_max = d%g_vap_max; limit_transport = d%limit_transport
       use_boiling_formula = d%use_boiling_formula; q_flux_max = d%q_flux_max
       charge_T_full = 0; charge_T_res = 0

       read(u, nml=heat_pipe, iostat=ios, iomsg=iomsg_)
       if (ios /= 0) then
          ierr = 3; write(msg,'(a,i0,a,a)') 'bad &heat_pipe #', k, ': ', trim(iomsg_)
          close(u); return
       end if
       if (ne < 1 .or. nc < 1 .or. ne > MAX_SEG .or. nc > MAX_SEG) then
          ierr = 4; write(msg,'(a,a)') 'ne/nc out of range in ', trim(name); close(u); return
       end if
       if (vap_node_e == 0) then
          ierr = 5; write(msg,'(a,a)') 'vap_node_e required in ', trim(name); close(u); return
       end if
       if (vap_node_c == 0) vap_node_c = vap_node_e    ! 증기 노드 1개 모델

       associate (h => hps(k), c => hps(k)%cfg)
         h%name = name; h%ne = ne; h%nc = nc
         h%evap_nodes = evap_nodes; h%cond_nodes = cond_nodes
         h%evap_conds = evap_conds; h%cond_conds = cond_conds
         h%vap_node_e = vap_node_e; h%vap_node_c = vap_node_c
         h%res_node = res_node; h%vap_cond = vap_cond
         c%L_evap = L_evap; c%L_adia = L_adia; c%L_cond = L_cond
         c%r_vap = r_vap; c%groove_w = groove_w; c%groove_d = groove_d
         c%n_groove = n_groove; c%r_out = r_out; c%k_wall = k_wall
         c%contact_angle = contact_angle; c%elevation = elevation; c%g_accel = g_accel
         c%r_nuc = r_nuc; c%mult_evap = mult_evap; c%mult_cond = mult_cond
         c%v_res = v_res; c%blocked_frac = blocked_frac
         c%relax = relax; c%g_vap_max = g_vap_max; c%limit_transport = limit_transport
         c%use_boiling_formula = use_boiling_formula; c%q_flux_max = q_flux_max
         ! NCG 충전량: 직접 입력(n_gas) 또는 설계점(charge_T_full, charge_T_res)으로 계산
         if (charge_T_full > 0 .and. charge_T_res > 0) then
            c%n_gas = vchp_charge(charge_T_full, charge_T_res, v_res)
         else
            c%n_gas = n_gas
         end if
         if (c%n_gas > 0 .and. res_node == 0) then
            ierr = 6; write(msg,'(a,a)') 'VCHP needs res_node in ', trim(name); close(u); return
         end if
         call hp_init(c, dx_e(1:ne), dx_c(1:nc), h%st, ierr, msg)
         if (ierr /= 0) then
            msg = trim(name)//': '//trim(msg); close(u); return
         end if
       end associate
    end do
    close(u)
    ready = .true.
  end subroutine hp_driver_init

  !=============================================================================
  ! 매 반복 호출: 온도 읽기 -> 모델 갱신 -> 전도체 쓰기
  !=============================================================================
  subroutine hp_driver_update(get_temp, set_cond)
    procedure(get_temp_if) :: get_temp
    procedure(set_cond_if) :: set_cond
    integer :: k, i
    real(dp) :: Te(MAX_SEG), Tc(MAX_SEG), Tve, Tvc, Tres
    if (.not. ready) return
    do k = 1, n_hp
       associate (h => hps(k))
         do i = 1, h%ne
            Te(i) = get_temp(h%evap_nodes(i))
         end do
         do i = 1, h%nc
            Tc(i) = get_temp(h%cond_nodes(i))
         end do
         Tve = get_temp(h%vap_node_e)
         Tvc = get_temp(h%vap_node_c)
         Tres = Tvc
         if (h%res_node /= 0) Tres = get_temp(h%res_node)

         call hp_update(h%cfg, h%st, Te(1:h%ne), Tve, Tvc, Tc(1:h%nc), Tres)

         do i = 1, h%ne
            call set_cond(h%evap_conds(i), h%st%g_evap(i))
         end do
         do i = 1, h%nc
            call set_cond(h%cond_conds(i), h%st%g_cond(i))
         end do
         if (h%vap_cond /= 0) call set_cond(h%vap_cond, h%st%g_vap)
       end associate
    end do
  end subroutine hp_driver_update

  !=============================================================================
  ! 진단 출력 (시간 스텝 끝 등에서 호출)
  !=============================================================================
  subroutine hp_driver_report(unit, time)
    integer, intent(in) :: unit
    real(dp), intent(in) :: time
    integer :: k
    do k = 1, n_hp
       associate (s => hps(k)%st)
         write(unit,'(a,es12.4,1x,a12,a,f8.2,a,f8.2,a,f8.2,a,a,a,f7.3,a,f6.3,3(a,l1))') &
              't=', time, trim(hps(k)%name), ' Tv[C]=', s%t_vap - 273.15_dp, &
              ' Q[W]=', s%q_transport, ' Qmax[W]=', s%q_max, &
              ' (', trim(LIMIT_NAME(s%limit_id)), ') Lact[m]=', s%L_active, &
              ' scale=', s%scale, ' limited=', s%limited, ' shut=', s%shut_off, &
              ' oor=', s%out_of_range
       end associate
    end do
  end subroutine hp_driver_report

end module hp_driver
