!===============================================================================
! hp_props : 암모니아 포화 물성 (테이블 보간)
!
!  - 테이블: hp_props_data.inc (tools/gen_ammonia_table.py 로 CoolProp 에서 생성)
!  - 보간: 포화압력은 ln(P) vs 1/T 선형 (Clausius-Clapeyron 형태라 정확도 높음),
!          나머지 물성은 T 에 대한 선형 보간
!  - 범위 밖 온도는 테이블 끝값으로 고정(clip)하고 out_of_range = .true. 반환
!===============================================================================
module hp_props
  implicit none
  private
  integer, parameter, public :: dp = kind(1.0d0)

  include 'hp_props_data.inc'

  real(dp), parameter, public :: R_UNIV = 8.314462618_dp        ! [J/mol.K]
  real(dp), parameter, public :: M_NH3  = 0.01703052_dp         ! [kg/mol]

  type, public :: fluid_props
     real(dp) :: T      = 0.0_dp   ! [K]
     real(dp) :: psat   = 0.0_dp   ! [Pa]
     real(dp) :: rho_l  = 0.0_dp   ! [kg/m3]
     real(dp) :: rho_v  = 0.0_dp   ! [kg/m3]
     real(dp) :: h_fg   = 0.0_dp   ! [J/kg]
     real(dp) :: mu_l   = 0.0_dp   ! [Pa.s]
     real(dp) :: mu_v   = 0.0_dp   ! [Pa.s]
     real(dp) :: k_l    = 0.0_dp   ! [W/m.K]
     real(dp) :: sigma  = 0.0_dp   ! [N/m]
     real(dp) :: gamma_v = 0.0_dp  ! [-]  증기 비열비
     logical  :: out_of_range = .false.
  end type fluid_props

  public :: get_props, psat_of_T, props_tmin, props_tmax

contains

  pure real(dp) function props_tmin()
    props_tmin = TAB_T(1)
  end function props_tmin

  pure real(dp) function props_tmax()
    props_tmax = TAB_T(NTAB)
  end function props_tmax

  ! 보간 구간 인덱스와 가중치
  pure subroutine locate(T, i, w, Tc, oor)
    real(dp), intent(in)  :: T
    integer,  intent(out) :: i
    real(dp), intent(out) :: w, Tc
    logical,  intent(out) :: oor
    oor = (T < TAB_T(1)) .or. (T > TAB_T(NTAB))
    Tc = min(max(T, TAB_T(1)), TAB_T(NTAB))
    i = int((Tc - TAB_TMIN) / TAB_DT) + 1
    i = min(max(i, 1), NTAB - 1)
    w = (Tc - TAB_T(i)) / (TAB_T(i+1) - TAB_T(i))
  end subroutine locate

  pure real(dp) function lin(tab, i, w)
    real(dp), intent(in) :: tab(:), w
    integer,  intent(in) :: i
    lin = tab(i) + w * (tab(i+1) - tab(i))
  end function lin

  pure real(dp) function psat_interp(Tc, i)
    real(dp), intent(in) :: Tc
    integer,  intent(in) :: i
    real(dp) :: x, x1, x2
    x  = 1.0_dp / Tc
    x1 = 1.0_dp / TAB_T(i)
    x2 = 1.0_dp / TAB_T(i+1)
    psat_interp = exp(log(TAB_PSAT(i)) + (x - x1) / (x2 - x1) &
                      * (log(TAB_PSAT(i+1)) - log(TAB_PSAT(i))))
  end function psat_interp

  !> 포화압력 [Pa]
  pure real(dp) function psat_of_T(T)
    real(dp), intent(in) :: T
    integer :: i
    real(dp) :: w, Tc
    logical :: oor
    call locate(T, i, w, Tc, oor)
    psat_of_T = psat_interp(Tc, i)
  end function psat_of_T

  !> 전체 포화 물성
  pure function get_props(T) result(p)
    real(dp), intent(in) :: T
    type(fluid_props) :: p
    integer :: i
    real(dp) :: w, Tc
    call locate(T, i, w, Tc, p%out_of_range)
    p%T       = Tc
    p%psat    = psat_interp(Tc, i)
    p%rho_l   = lin(TAB_RHO_L, i, w)
    p%rho_v   = lin(TAB_RHO_V, i, w)
    p%h_fg    = lin(TAB_H_FG, i, w)
    p%mu_l    = lin(TAB_MU_L, i, w)
    p%mu_v    = lin(TAB_MU_V, i, w)
    p%k_l     = lin(TAB_K_L, i, w)
    p%sigma   = lin(TAB_SIGMA, i, w)
    p%gamma_v = lin(TAB_GAMMA_V, i, w)
  end function get_props

end module hp_props
