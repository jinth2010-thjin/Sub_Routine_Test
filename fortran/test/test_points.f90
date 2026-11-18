!===============================================================================
! test_points : 물리 함수 점검용 출력 (python/check_fortran.py 가 참조 모델과 비교)
!  출력: build/points_hp.csv, build/points_front.csv
!===============================================================================
program test_points
  use hp_props
  use hp_model
  implicit none
  type(hp_config) :: cfg
  type(fluid_props) :: p
  real(dp), parameter :: Ts(6) = [220.0_dp, 250.0_dp, 273.15_dp, 300.0_dp, 330.0_dp, 360.0_dp]
  real(dp), parameter :: elev(2) = [0.0_dp, 0.01_dp]
  real(dp) :: L_eff, qs, qe, qb, qv, dx(10), Tc(10), f(10), ngas
  integer :: i, k, u, c
  logical :: shut

  open(newunit=u, file='build/points_hp.csv', status='replace')
  write(u,'(a)') 'T,elev,psat,rho_l,rho_v,h_fg,mu_l,mu_v,k_l,sigma,gamma_v,keff,gpl,' // &
                 'dp_l100,dp_v100,dp_v5000,qcap,qsonic,qentr,qboil,qvisc'
  L_eff = cfg%L_adia + 0.5_dp * (cfg%L_evap + cfg%L_cond)
  do k = 1, 2
     cfg%elevation = elev(k)
     do i = 1, size(Ts)
        p = get_props(Ts(i))
        call other_limits(cfg, p, L_eff, qs, qe, qb, qv)
        write(u,'(f8.3,",",f6.3,19(",",es20.12))') Ts(i), elev(k), p%psat, p%rho_l, p%rho_v, &
             p%h_fg, p%mu_l, p%mu_v, p%k_l, p%sigma, p%gamma_v, groove_keff(cfg, p%k_l), &
             conductance_per_length(cfg, p), liquid_dp(cfg, p, 100.0_dp, L_eff), &
             vapor_dp(cfg, p, 100.0_dp, L_eff), vapor_dp(cfg, p, 5000.0_dp, L_eff), &
             capillary_limit(cfg, p, L_eff), qs, qe, qb, qv
     end do
  end do
  close(u)

  ! VCHP front: 저장조 30 cm3, 설계 충전 (T_full=303.15 K, T_res=253.15 K)
  dx = 0.06_dp
  ngas = vchp_charge(303.15_dp, 253.15_dp, 30.0e-6_dp)
  open(newunit=u, file='build/points_front.csv', status='replace')
  write(u,'(a)') 'case,n_gas,T_vc,T_res,Tc1,Tc2,Tc3,Tc4,Tc5,Tc6,Tc7,Tc8,Tc9,Tc10,' // &
                 'f1,f2,f3,f4,f5,f6,f7,f8,f9,f10,shut'
  do c = 1, 6
     select case (c)
     case (1); Tc = 253.15_dp;  call front(293.15_dp, 253.15_dp)
     case (2); Tc = 253.15_dp;  call front(298.15_dp, 253.15_dp)
     case (3); Tc = 253.15_dp;  call front(305.15_dp, 253.15_dp)
     case (4); Tc = [(283.15_dp - 4.0_dp*i, i=1,10)]; call front(296.15_dp, 240.0_dp)
     case (5); Tc = 233.15_dp;  call front(288.15_dp, 233.15_dp)
     case (6); Tc = 270.0_dp;   call front(272.0_dp, 270.0_dp)
     end select
  end do
  close(u)
  print '(a)', 'wrote build/points_hp.csv, build/points_front.csv'
contains
  subroutine front(Tvc, Tres)
    real(dp), intent(in) :: Tvc, Tres
    call ncg_front(ngas, 30.0e-6_dp, 3.14159265358979_dp * cfg%r_vap**2, Tvc, Tres, Tc, dx, f, shut)
    write(u,'(i2,",",es20.12,22(",",es20.12),",",l1)') c, ngas, Tvc, Tres, Tc, f, shut
  end subroutine front
end program test_points
