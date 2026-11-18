"""Fortran 결과(build/points_*.csv)를 Python 참조 모델과 비교.

실행: python3 python/check_fortran.py   (먼저 make test_points 실행)
종료코드 0 = 통과
"""
import csv
import math
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
import hp_reference as ref  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
RTOL = 1e-6
fails = 0


def check(name, got, exp, rtol=RTOL, atol=1e-12):
    global fails
    ok = abs(got - exp) <= atol + rtol * abs(exp)
    if not ok:
        fails += 1
        print(f"  FAIL {name}: fortran={got:.10g} python={exp:.10g}")
    return ok


def main():
    c = ref.Config()
    n = 0
    with open(ROOT / "build" / "points_hp.csv") as f:
        for r in csv.DictReader(f, skipinitialspace=True):
            r = {k: float(v) for k, v in r.items()}
            c.elevation = r["elev"]
            p = ref.props(r["T"])
            L = c.L_eff_full
            tag = f"T={r['T']:.2f},elev={r['elev']}"
            for k in ["psat", "rho_l", "rho_v", "h_fg", "mu_l", "mu_v", "k_l", "sigma", "gamma_v"]:
                check(f"{tag} {k}", r[k], getattr(p, k))
            check(f"{tag} keff", r["keff"], ref.keff(c, p.k_l))
            check(f"{tag} gpl", r["gpl"], ref.gpl(c, p))
            check(f"{tag} dp_l100", r["dp_l100"], 100 * ref.liquid_dp_per_W(c, p, L))
            check(f"{tag} dp_v100", r["dp_v100"], ref.vapor_dp(c, p, 100, L))
            check(f"{tag} dp_v5000", r["dp_v5000"], ref.vapor_dp(c, p, 5000, L))
            check(f"{tag} qcap", r["qcap"], ref.capillary_limit(c, p, L), rtol=1e-7)
            for k, v in zip(["qsonic", "qentr", "qboil", "qvisc"], ref.other_limits(c, p, L)):
                check(f"{tag} {k}", r[k], v)
            n += 1
    print(f"points_hp: {n} rows checked")

    a_vap = math.pi * c.r_vap**2
    v_res = 30e-6
    n = 0
    with open(ROOT / "build" / "points_front.csv") as f:
        for r in csv.DictReader(f, skipinitialspace=True):
            shut = r.pop("shut").strip() == "T"
            r = {k: float(v) for k, v in r.items()}
            Tc = [r[f"Tc{i}"] for i in range(1, 11)]
            dx = [0.06] * 10
            fr, sh = ref.ncg_front(r["n_gas"], v_res, a_vap, r["T_vc"], r["T_res"], Tc, dx)
            for i in range(10):
                check(f"case{int(r['case'])} f{i+1}", r[f"f{i+1}"], fr[i])
            if sh != shut:
                print(f"  FAIL case{int(r['case'])} shut flag")
                globals()["fails"] += 1
            # 균일 온도 케이스는 해석해와 비교
            if all(abs(t - r["T_res"]) < 1e-9 for t in Tc) and not shut:
                lb = ref.blocked_length_uniform(r["n_gas"], v_res, a_vap, r["T_vc"], r["T_res"])
                lb_f = sum((1 - r[f"f{i}"]) * 0.06 for i in range(1, 11))
                check(f"case{int(r['case'])} blocked length vs analytic", lb_f, max(lb, 0.0), atol=1e-9)
                print(f"  case{int(r['case'])}: blocked length {lb_f*100:.3f} cm (analytic {max(lb,0)*100:.3f} cm)")
            n += 1
    print(f"points_front: {n} cases checked")

    print("RESULT:", "PASS" if fails == 0 else f"FAIL ({fails})")
    return 0 if fails == 0 else 1


if __name__ == "__main__":
    sys.exit(main())
