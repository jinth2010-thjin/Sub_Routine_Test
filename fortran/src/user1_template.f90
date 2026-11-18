!===============================================================================
! user1_template : Space Systems Thermal (TMG 솔버) USER1 서브루틴 연결 템플릿
!
! !!! 반드시 읽어 주세요 !!!
!  이 파일은 "연결 지점"만 정의한 템플릿입니다. TMG USER1 의 정확한 인자 목록,
!  호출 시점 플래그, 노드 온도를 읽고 전도체 값을 바꾸는 API 이름은 이 저장소를
!  만들 때 문서를 확인하지 못했습니다(help.mayahtt.com 접근 불가).
!  NX2027 설치본의 도움말(Thermal API / User-written subroutine example)에서
!  아래 [TODO-1] ~ [TODO-3] 을 실제 내용으로 채워야 합니다.
!
!  확인된 사항 (Siemens/Maya HTT 공개 문서 요약 기준)
!   - USER1 : 해석 도중 열 모델을 수정하는 사용자 Fortran 서브루틴
!   - 고급 파라미터 ACTIVATE USER DLL SCHEME 로 동적 링크(DLL) 방식 활성화
!   - Solution 대화상자 > Thermal 탭 > Include File 그룹에서 파일 지정
!
! 사용 방법
!   1) 실제 USER1 서브루틴(문서 예제의 시그니처 그대로) 안에서
!         call hp_user1_entry(phase, time)
!      를 호출합니다. phase 는 아래 PHASE_* 중 하나로 변환해 넘깁니다.
!   2) tmg_get_temp / tmg_set_cond 의 본문을 TMG API 호출로 구현합니다.
!   3) hp_props, hp_model, hp_driver, 이 파일을 함께 컴파일해 DLL 로 링크합니다.
!   4) 해석 작업 폴더에 hp_input.nml 을 둡니다 (또는 HP_INPUT_FILE 수정).
!===============================================================================
module hp_user1_bridge
  use hp_props, only: dp
  use hp_driver
  implicit none
  private

  integer, parameter, public :: PHASE_INIT = 1       ! 해석 시작 (모델 로드 후 1회)
  integer, parameter, public :: PHASE_ITERATION = 2  ! 매 반복 (전도체 갱신)
  integer, parameter, public :: PHASE_STEP_END = 3   ! 시간 스텝 / 정상상태 수렴 후
  integer, parameter, public :: PHASE_FINAL = 4      ! 해석 종료

  character(len=*), parameter :: HP_INPUT_FILE = 'hp_input.nml'
  character(len=*), parameter :: HP_LOG_FILE   = 'hp_subroutine.log'
  integer, save :: log_unit = -1

  public :: hp_user1_entry

contains

  recursive subroutine hp_user1_entry(phase, time)
    integer, intent(in) :: phase
    real(dp), intent(in) :: time
    integer :: ierr
    character(len=256) :: msg

    select case (phase)
    case (PHASE_INIT)
       open(newunit=log_unit, file=HP_LOG_FILE, status='replace', action='write')
       call hp_driver_init(HP_INPUT_FILE, ierr, msg)
       if (ierr /= 0) then
          write(log_unit,'(a)') 'HP SUBROUTINE INIT ERROR: '//trim(msg)
          ! [TODO] 솔버에 오류를 알리고 중단하는 TMG 방식이 있으면 사용
          error stop 'hp_input.nml error - see hp_subroutine.log'
       end if
       write(log_unit,'(a,i0)') 'heat pipes loaded: ', hp_driver_count()

    case (PHASE_ITERATION)
       if (.not. hp_driver_ready()) call hp_user1_entry(PHASE_INIT, time)
       call hp_driver_update(tmg_get_temp, tmg_set_cond)

    case (PHASE_STEP_END)
       if (log_unit > 0) call hp_driver_report(log_unit, time)

    case (PHASE_FINAL)
       if (log_unit > 0) close(log_unit)
    end select
  end subroutine hp_user1_entry

  !-----------------------------------------------------------------------------
  ! [TODO-2] 노드 ID -> 온도 [K]
  !  TMG API 로 노드(요소) 번호의 현재 반복 온도를 읽어 반환하도록 구현.
  !  온도 단위가 C 이면 +273.15 하세요.
  !-----------------------------------------------------------------------------
  function tmg_get_temp(node_id) result(T)
    integer, intent(in) :: node_id
    real(dp) :: T
    T = 0.0_dp
    error stop 'tmg_get_temp not implemented: fill in TMG API call (user1_template.f90)'
    if (node_id < 0) T = -1.0_dp   ! (미사용 경고 방지)
  end function tmg_get_temp

  !-----------------------------------------------------------------------------
  ! [TODO-3] 전도체 ID 에 값 [W/K] 설정
  !  TMG API 로 해당 전도체(conductance) 값을 덮어쓰도록 구현.
  !  솔버가 전도체 값을 반복마다 다시 읽는지(행렬 재조립) 문서에서 꼭 확인하세요.
  !-----------------------------------------------------------------------------
  subroutine tmg_set_cond(cond_id, value)
    integer, intent(in) :: cond_id
    real(dp), intent(in) :: value
    error stop 'tmg_set_cond not implemented: fill in TMG API call (user1_template.f90)'
    if (cond_id < 0 .and. value < 0) return
  end subroutine tmg_set_cond

end module hp_user1_bridge

!-------------------------------------------------------------------------------
! [TODO-1] 실제 USER1 서브루틴
!  NX2027 도움말의 USER1 예제에서 시그니처(인자 목록)를 그대로 복사하고,
!  호출 위치 플래그를 PHASE_* 로 변환해 hp_user1_entry 를 호출하세요.
!  예 (가상의 형태 - 실제 인자 이름/의미는 문서 기준으로 교체):
!
!     SUBROUTINE USER1(...)
!       USE hp_user1_bridge
!       ...
!       IF (<해석 시작 플래그>)  CALL hp_user1_entry(PHASE_INIT, TIME)
!       IF (<반복 플래그>)       CALL hp_user1_entry(PHASE_ITERATION, TIME)
!       IF (<스텝 종료 플래그>)  CALL hp_user1_entry(PHASE_STEP_END, TIME)
!       IF (<해석 종료 플래그>)  CALL hp_user1_entry(PHASE_FINAL, TIME)
!     END SUBROUTINE USER1
!-------------------------------------------------------------------------------
