# Sub_Routine_Test
NX의 Subroutine 기능을 알아봅니다.

## CCHP / VCHP 히트파이프 사용자 서브루틴 (Space Systems Thermal)

작동유체(암모니아)의 유동·상변화 현상을 반영한 축방향 그루브 히트파이프 모델을 Fortran 사용자 서브루틴으로 구현했습니다.

- 그루브 액체층 전도(증발/응축 막 저항), 증기 압력강하, 모세관·음속·비말동반·점성 한계, dry-out 처리
- VCHP 비응축가스(NCG) flat-front 모델 (구간별 온도, 증기 분압 반영)
- 입력 파일(namelist)로 여러 개의 히트파이프와 노드/전도체 ID를 정의

자세한 이론, 검증, 연결 절차는 **[docs/heat_pipe_model.md](docs/heat_pipe_model.md)**에 있습니다.

```
fortran/src/hp_props.f90        암모니아 포화 물성 (테이블 보간)
fortran/src/hp_model.f90        CCHP/VCHP 물리 모델 (솔버 독립)
fortran/src/hp_driver.f90       namelist 입력 + 노드/전도체 매핑
fortran/src/user1_template.f90  TMG USER1 연결 템플릿 (TODO 3곳 구현 필요)
fortran/test/                   점검 프로그램, 독립 과도해석 데모
python/                         참조 모델(교차 검증), 그래프
examples/hp_input.nml           입력 파일 예제
```

### 빌드와 검증 (gfortran, Python 3 + matplotlib)

```bash
make test
```

```bash
make demo
```
