# 히트파이프 사용자 서브루틴 - 로컬 빌드/검증 (gfortran 기준)
#   make            : 모듈 + 테스트 + 데모 빌드
#   make test       : Fortran vs Python 교차 검증 + 드라이버 테스트
#   make demo       : 과도해석 데모 6케이스 실행 + 그래프 (docs/img)
#   make table      : 암모니아 물성 테이블 재생성 (CoolProp 필요)
FC      ?= gfortran
FFLAGS  ?= -O2 -Wall -Wextra -std=f2008 -fcheck=all
B       := fortran/build
SRC     := fortran/src
MODS    := $(B)/hp_props.o $(B)/hp_model.o $(B)/hp_driver.o
PY      ?= python3
CASES   := cchp_hot vchp_hot cchp_cold vchp_cold cchp_pulse cchp_tilt

all: $(B)/test_points $(B)/test_driver $(B)/demo_transient $(B)/user1_template.o

$(B):
	mkdir -p $(B) build

$(B)/hp_props.o: $(SRC)/hp_props.f90 $(SRC)/hp_props_data.inc | $(B)
	$(FC) $(FFLAGS) -I$(SRC) -J$(B) -c $< -o $@
$(B)/hp_model.o: $(SRC)/hp_model.f90 $(B)/hp_props.o
	$(FC) $(FFLAGS) -I$(B) -J$(B) -c $< -o $@
$(B)/hp_driver.o: $(SRC)/hp_driver.f90 $(B)/hp_model.o
	$(FC) $(FFLAGS) -I$(B) -J$(B) -c $< -o $@
$(B)/user1_template.o: $(SRC)/user1_template.f90 $(B)/hp_driver.o
	$(FC) $(FFLAGS) -Wno-unused -I$(B) -J$(B) -c $< -o $@

$(B)/%: fortran/test/%.f90 $(MODS)
	$(FC) $(FFLAGS) -I$(B) -J$(B) $(MODS) $< -o $@

test: $(B)/test_points $(B)/test_driver
	./$(B)/test_points
	$(PY) python/check_fortran.py
	./$(B)/test_driver

demo: $(B)/demo_transient
	@for c in $(CASES); do ./$(B)/demo_transient $$c; done
	$(PY) python/plot_demo.py

table:
	$(PY) tools/gen_ammonia_table.py

clean:
	rm -rf $(B) build

.PHONY: all test demo table clean
