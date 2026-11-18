"""히트파이프 물리 모델 Python 참조 구현 (Fortran hp_model.f90 검증용).

Fortran 과 같은 식을 독립적으로 다시 구현한 것입니다. 모세관 한계는 층류
구간에서 닫힌 해(선형식)로, VCHP front 는 균일 온도일 때의 해석해로도 계산해
Fortran 의 수치해(이분법, 구간 행진)를 별도로 확인합니다.
"""
from __future__ import annotations

import math
from dataclasses import dataclass
from pathlib import Path

R_UNIV = 8.314462618
M_NH3 = 0.01703052
PI = math.pi

_TABLE = Path(__file__).with_name("ammonia_table.csv")


def _load_table():
    rows = [l for l in _TABLE.read_text().splitlines() if l and not l.startswith("#")]
    cols = rows[0].split(",")
    data = {c: [] for c in cols}
    for r in rows[1:]:
        for c, v in zip(cols, r.split(",")):
            data[c].append(float(v))
    return data


TAB = _load_table()


@dataclass
class Props:
    T: float
    psat: float
    rho_l: float
    rho_v: float
    h_fg: float
    mu_l: float
    mu_v: float
    k_l: float
    sigma: float
    gamma_v: float


def _locate(T):
    t = TAB["T"]
    Tc = min(max(T, t[0]), t[-1])
    i = min(max(int((Tc - t[0]) / (t[1] - t[0])), 0), len(t) - 2)
    return Tc, i, (Tc - t[i]) / (t[i + 1] - t[i])


def psat(T):
    Tc, i, _ = _locate(T)
    t, p = TAB["T"], TAB["Psat"]
    x, x1, x2 = 1 / Tc, 1 / t[i], 1 / t[i + 1]
    return math.exp(math.log(p[i]) + (x - x1) / (x2 - x1) * (math.log(p[i + 1]) - math.log(p[i])))


def props(T) -> Props:
    Tc, i, w = _locate(T)
    lin = lambda k: TAB[k][i] + w * (TAB[k][i + 1] - TAB[k][i])
    return Props(Tc, psat(T), lin("rho_l"), lin("rho_v"), lin("h_fg"), lin("mu_l"),
                 lin("mu_v"), lin("k_l"), lin("sigma"), lin("gamma_v"))


@dataclass
class Config:
    L_evap: float = 0.30
    L_adia: float = 0.40
    L_cond: float = 0.60
    r_vap: float = 4.0e-3
    groove_w: float = 0.6e-3
    groove_d: float = 1.2e-3
    n_groove: int = 24
    r_out: float = 6.35e-3
    k_wall: float = 200.0
    contact_angle: float = 0.0
    elevation: float = 0.0
    g_accel: float = 9.80665
    r_nuc: float = 2.54e-7

    @property
    def L_eff_full(self):
        return self.L_adia + 0.5 * (self.L_evap + self.L_cond)


def keff(c: Config, k_l):
    w, d, kw = c.groove_w, c.groove_d, c.k_wall
    wf = 2 * PI * c.r_vap / c.n_groove - w
    return (wf * k_l * kw * d + w * k_l * (0.185 * wf * kw + d * k_l)) / (
        (w + wf) * (0.185 * wf * kw + d * k_l))


def gpl(c: Config, p: Props):
    rg = c.r_vap + c.groove_d
    return 2 * PI / (math.log(c.r_out / rg) / c.k_wall + math.log(rg / c.r_vap) / keff(c, p.k_l))


def fre_rect(a):
    return 24 * (1 - 1.3553 * a + 1.9467 * a**2 - 1.7012 * a**3 + 0.9564 * a**4 - 0.2537 * a**5)


def liquid_dp_per_W(c: Config, p: Props, L_eff):
    w, d = c.groove_w, c.groove_d
    A = c.n_groove * w * d
    Dh = 4 * w * d / (w + 2 * d)
    a = min(w, 2 * d) / max(w, 2 * d)
    return 2 * fre_rect(a) * p.mu_l * L_eff / (p.h_fg * p.rho_l * A * Dh**2)


def vapor_dp(c: Config, p: Props, Q, L_eff):
    A, D = PI * c.r_vap**2, 2 * c.r_vap
    mdot = Q / p.h_fg
    u = mdot / (p.rho_v * A)
    re = mdot * D / (A * p.mu_v)
    if re < 2300:
        return 32 * p.mu_v * u * L_eff / D**2
    return 2 * 0.079 * re**-0.25 * p.rho_v * u**2 * L_eff / D


def capillary_limit(c: Config, p: Props, L_eff):
    """층류면 닫힌 해, 난류면 이분법."""
    dpc = 2 * p.sigma * math.cos(c.contact_angle) / c.groove_w
    dpg = p.rho_l * c.g_accel * c.elevation
    if dpg >= dpc:
        return 0.0
    kv = vapor_dp(c, p, 1.0, L_eff)            # 층류 기울기 [Pa/W]
    q = (dpc - dpg) / (liquid_dp_per_W(c, p, L_eff) + kv)
    re = (q / p.h_fg) * 2 * c.r_vap / (PI * c.r_vap**2 * p.mu_v)
    if re < 2300:
        return q
    lo, hi = 0.0, q
    f = lambda Q: liquid_dp_per_W(c, p, L_eff) * Q + vapor_dp(c, p, Q, L_eff) + dpg - dpc
    while f(hi) < 0:
        hi *= 2
    for _ in range(200):
        m = 0.5 * (lo + hi)
        lo, hi = (lo, m) if f(m) > 0 else (m, hi)
    return 0.5 * (lo + hi)


def other_limits(c: Config, p: Props, L_eff):
    A = PI * c.r_vap**2
    Rg = R_UNIV / M_NH3
    rg = c.r_vap + c.groove_d
    dpc = 2 * p.sigma * math.cos(c.contact_angle) / c.groove_w
    sonic = A * p.rho_v * p.h_fg * math.sqrt(p.gamma_v * Rg * p.T / (2 * (p.gamma_v + 1)))
    entr = A * p.h_fg * math.sqrt(p.sigma * p.rho_v / c.groove_w)
    boil = (2 * PI * c.L_evap * keff(c, p.k_l) * p.T / (p.h_fg * p.rho_v * math.log(rg / c.r_vap))
            * max(2 * p.sigma / c.r_nuc - dpc, 0.0))
    visc = A * c.r_vap**2 * p.h_fg * p.rho_v * p.psat / (16 * p.mu_v * L_eff)
    return sonic, entr, boil, visc


def vchp_charge(T_full, T_res, v_res):
    return (psat(T_full) - psat(min(T_res, T_full))) * v_res / (R_UNIV * T_res)


def ncg_front(n_gas, v_res, a_vap, T_vc, T_res, Tc, dx):
    pv = psat(T_vc)

    def moles(T, V):
        pg = max(pv - psat(min(T, T_vc)), 1e-6 * pv)
        return pg * V / (R_UNIV * T)

    f = [1.0] * len(Tc)
    remain = n_gas - moles(T_res, v_res)
    if remain <= 0:
        return f, False
    for j in range(len(Tc) - 1, -1, -1):
        nk = moles(Tc[j], a_vap * dx[j])
        if remain >= nk:
            f[j] = 0.0
            remain -= nk
        else:
            f[j] = 1 - remain / nk
            remain = 0.0
            break
    return f, remain > 0


def blocked_length_uniform(n_gas, v_res, a_vap, T_vc, T_s):
    """응축부와 저장조가 모두 T_s 일 때의 해석해 (이론식)."""
    pg = psat(T_vc) - psat(T_s)
    return (n_gas * R_UNIV * T_s / pg - v_res) / a_vap
