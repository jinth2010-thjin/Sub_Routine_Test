!===============================================================================
! test_driver : hp_driver (namelist 입력 + 노드/전도체 ID 매핑) 점검
!  가짜 솔버(노드 ID -> 온도 배열)로 hp_driver_update 를 호출하고, 같은 입력으로
!  hp_update 를 직접 호출한 결과와 전도체 값이 같은지 비교한다.
!===============================================================================
module fake_solver
  use hp_props, only: dp
  implicit none
  real(dp) :: temp(9999) = 280.0_dp
  real(dp) :: cond(9999) = -1.0_dp
contains
  function get_temp(id) result(T)
    integer, intent(in) :: id
    real(dp) :: T
    T = temp(id)
  end function get_temp
  subroutine set_cond(id, v)
    integer, intent(in) :: id
    real(dp), intent(in) :: v
    cond(id) = v
  end subroutine set_cond
end module fake_solver

program test_driver
  use hp_props
  use hp_model
  use hp_driver
  use fake_solver
  implicit none
  integer :: ierr, i, nfail
  character(len=256) :: msg
  type(hp_config) :: cfg
  type(hp_state) :: st
  real(dp) :: dx_c(10), Tc(10)

  call hp_driver_init('examples/hp_input.nml', ierr, msg)
  if (ierr /= 0) then
     print *, 'init failed: ', trim(msg); stop 1
  end if
  print '(a,i0)', 'heat pipes loaded: ', hp_driver_count()

  ! 온도 설정 (VCHP: 증발부 295 K, 응축부 271~262 K, 저장조 262 K -> front 가 응축부 중간)
  temp(1001:1002) = 300.0_dp; temp(1100) = 297.0_dp; temp(1011:1014) = 290.0_dp
  temp(2001:2002) = 295.0_dp; temp(2100) = 292.0_dp; temp(2101) = 291.9_dp
  Tc = [(272.0_dp - i, i = 1, 10)]
  temp(2011:2020) = Tc
  temp(2200) = 262.0_dp

  call hp_driver_update(get_temp, set_cond)
  call hp_driver_report(6, 0.0_dp)

  ! 같은 조건으로 hp_model 을 직접 호출
  dx_c = 0.06_dp
  cfg%v_res = 30.0e-6_dp
  cfg%n_gas = vchp_charge(293.15_dp, 263.15_dp, cfg%v_res)
  call hp_init(cfg, [0.15_dp, 0.15_dp], dx_c, st, ierr, msg)
  call hp_update(cfg, st, [295.0_dp, 295.0_dp], 292.0_dp, 291.9_dp, Tc, 262.0_dp)

  nfail = 0
  do i = 1, 10
     if (abs(cond(6010+i) - st%g_cond(i)) > 1.0e-12_dp) nfail = nfail + 1
  end do
  if (abs(cond(6001) - st%g_evap(1)) > 1.0e-12_dp) nfail = nfail + 1
  if (abs(cond(6100) - st%g_vap) > 1.0e-9_dp) nfail = nfail + 1
  if (any(cond(5001:5002) <= 0) .or. any(cond(5011:5014) <= 0)) nfail = nfail + 1
  if (cond(0+1) /= -1.0_dp) nfail = nfail + 1   ! 관련 없는 ID 는 건드리지 않아야 함

  print '(a,10f7.3)', 'VCHP f_act    : ', st%f_act
  print '(a,10f7.3)', 'VCHP G_cond   : ', cond(6011:6020)
  if (nfail == 0) then
     print '(a)', 'RESULT: PASS'
  else
     print '(a,i0,a)', 'RESULT: FAIL (', nfail, ')'; stop 1
  end if
end program test_driver
