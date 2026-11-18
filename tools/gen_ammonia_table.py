"""암모니아 포화 물성 테이블 생성기 (CoolProp 기반).

생성 파일:
  - fortran/src/hp_props_data.inc : Fortran DATA 배열 (hp_props 모듈이 include)
  - python/ammonia_table.csv      : Python 참조 모델이 읽는 동일한 테이블

실행:  python3 tools/gen_ammonia_table.py
(CoolProp 이 필요합니다: pip install CoolProp)
"""
from pathlib import Path

import CoolProp
import CoolProp.CoolProp as CP

FLUID = "Ammonia"
T_MIN, T_MAX, DT = 200.0, 380.0, 5.0   # [K]  (삼중점 195.5 K, 임계점 405.6 K)

ROOT = Path(__file__).resolve().parents[1]
COLS = ["T", "Psat", "rho_l", "rho_v", "h_fg", "mu_l", "mu_v", "k_l", "sigma", "gamma_v"]
UNITS = ["K", "Pa", "kg/m3", "kg/m3", "J/kg", "Pa.s", "Pa.s", "W/m.K", "N/m", "-"]


def row(T):
    p = lambda key, q: CP.PropsSI(key, "T", T, "Q", q, FLUID)
    return [
        T,
        p("P", 0),
        p("D", 0),
        p("D", 1),
        p("H", 1) - p("H", 0),
        p("V", 0),
        p("V", 1),
        p("L", 0),
        p("I", 0),
        p("Cpmass", 1) / p("Cvmass", 1),
    ]


def main():
    n = int(round((T_MAX - T_MIN) / DT)) + 1
    rows = [row(T_MIN + i * DT) for i in range(n)]

    csv = ROOT / "python" / "ammonia_table.csv"
    with csv.open("w") as f:
        f.write(f"# Saturated ammonia, CoolProp {CoolProp.__version__}\n")
        f.write("# units: " + ",".join(UNITS) + "\n")
        f.write(",".join(COLS) + "\n")
        for r in rows:
            f.write(",".join(f"{v:.10e}" for v in r) + "\n")

    inc = ROOT / "fortran" / "src" / "hp_props_data.inc"
    with inc.open("w") as f:
        f.write(f"! 자동 생성 파일 - 직접 수정하지 마세요 (tools/gen_ammonia_table.py)\n")
        f.write(f"! Saturated ammonia, CoolProp {CoolProp.__version__}, "
                f"T = {T_MIN:.0f} .. {T_MAX:.0f} K, dT = {DT:.0f} K\n")
        f.write(f"integer, parameter :: NTAB = {n}\n")
        f.write(f"real(dp), parameter :: TAB_TMIN = {T_MIN:.1f}_dp, TAB_DT = {DT:.1f}_dp\n")
        for j, (name, unit) in enumerate(zip(COLS, UNITS)):
            f.write(f"! {name} [{unit}]\n")
            f.write(f"real(dp), parameter :: TAB_{name.upper()}(NTAB) = [ &\n")
            vals = [f"{r[j]:.10e}_dp" for r in rows]
            for k in range(0, n, 3):
                chunk = ", ".join(vals[k:k + 3])
                tail = ", &" if k + 3 < n else " ]"
                f.write(f"    {chunk}{tail}\n")
    print(f"wrote {csv.relative_to(ROOT)} and {inc.relative_to(ROOT)} ({n} rows)")


if __name__ == "__main__":
    main()
