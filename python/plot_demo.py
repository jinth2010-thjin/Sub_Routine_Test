"""데모 과도해석 결과(build/demo_*.csv) 그래프 -> docs/img/*.png

실행: python3 python/plot_demo.py   (먼저 make demo)
"""
import csv
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "docs" / "img"
BLUE, ORANGE, AQUA = "#2a78d6", "#eb6834", "#1baf7a"
INK, MUTED, GRID = "#1f1f1e", "#6b6a63", "#e4e3dd"
LAST_ORBITS = 3 * 95.0  # 마지막 3궤도만 표시 [min]

plt.rcParams.update({
    "font.size": 10, "axes.edgecolor": MUTED, "axes.labelcolor": INK,
    "xtick.color": MUTED, "ytick.color": MUTED, "axes.grid": True,
    "grid.color": GRID, "grid.linewidth": 0.8, "axes.spines.top": False,
    "axes.spines.right": False, "lines.linewidth": 2, "legend.frameon": False,
})


def load(case):
    with open(ROOT / "build" / f"demo_{case}.csv") as f:
        rows = list(csv.DictReader(f, skipinitialspace=True))
    t_end = float(rows[-1]["time_min"])
    rows = [r for r in rows if float(r["time_min"]) >= t_end - LAST_ORBITS]
    t0 = float(rows[0]["time_min"])
    col = lambda k: [float(r[k]) for r in rows]
    d = {k: col(k) for k in rows[0] if k != "limit"}
    d["t"] = [x - t0 for x in d["time_min"]]
    return d


def fig_box_temperature():
    fig, axes = plt.subplots(1, 2, figsize=(10, 3.8), sharey=True)
    for ax, env, title in [(axes[0], "hot", "Hot case"), (axes[1], "cold", "Cold case")]:
        for kind, color, label in [("cchp", BLUE, "CCHP"), ("vchp", ORANGE, "VCHP")]:
            d = load(f"{kind}_{env}")
            ax.plot(d["t"], d["T_box"], color=color, label=label)
            ax.annotate(label, (d["t"][-1], d["T_box"][-1]), xytext=(4, 0),
                        textcoords="offset points", va="center", color=INK, fontsize=9)
        ax.set_title(title, loc="left", color=INK)
        ax.set_xlabel("Time [min] (last 3 orbits)")
    axes[0].set_ylabel("Equipment box temperature [°C]")
    axes[0].legend(loc="lower left")
    fig.suptitle("Equipment temperature: CCHP vs VCHP", x=0.01, ha="left", color=INK)
    fig.tight_layout()
    fig.savefig(OUT / "box_temperature.png", dpi=150)


def fig_vchp_front():
    fig, axes = plt.subplots(2, 1, figsize=(10, 5.2), sharex=True)
    for env, color, label in [("hot", BLUE, "VCHP hot"), ("cold", ORANGE, "VCHP cold")]:
        d = load(f"vchp_{env}")
        axes[0].plot(d["t"], [100 * x for x in d["L_active"]], color=color, label=label)
        axes[1].plot(d["t"], d["T_vc"], color=color, label=label)
    axes[0].axhline(60, color=MUTED, lw=1, ls="--")
    axes[0].text(2, 61.5, "full condenser (60 cm)", color=MUTED, fontsize=8)
    axes[0].set_ylabel("Active condenser\nlength [cm]")
    axes[1].set_ylabel("Vapor temperature [°C]")
    axes[1].set_xlabel("Time [min] (last 3 orbits)")
    axes[0].legend(loc="center right")
    axes[0].set_title("VCHP: NCG front moves to hold vapor temperature", loc="left", color=INK)
    fig.tight_layout()
    fig.savefig(OUT / "vchp_front.png", dpi=150)


def fig_transport_limit():
    fig, axes = plt.subplots(2, 2, figsize=(10, 5.6), sharex=True)
    for col, case, title in [(0, "cchp_pulse", "400 W pulse (orbit)"),
                             (1, "cchp_tilt", "Ground test, 10 mm adverse tilt")]:
        d = load(case)
        ax = axes[0][col]
        ax.plot(d["t"], d["Q_hp"], color=BLUE, label="Transported Q")
        ax.plot(d["t"], d["Q_max"], color=ORANGE, label="Transport limit Q_max")
        ax.plot(d["t"], d["q_box"], color=MUTED, lw=1, ls="--", label="Box dissipation")
        ax.set_title(title, loc="left", color=INK)
        ax2 = axes[1][col]
        ax2.plot(d["t"], d["T_box"], color=BLUE, label="Box")
        ax2.plot(d["t"], d["T_ve"], color=AQUA, label="Vapor")
        ax2.set_xlabel("Time [min] (last 3 orbits)")
    axes[0][0].set_ylabel("Heat [W]")
    axes[1][0].set_ylabel("Temperature [°C]")
    h1, l1 = axes[0][0].get_legend_handles_labels()
    h2, l2 = axes[1][0].get_legend_handles_labels()
    fig.legend(h1 + h2, l1 + [f"T {x}" for x in l2], loc="upper right", ncol=5, fontsize=8)
    fig.suptitle("Transport limit (dry-out) handling", x=0.01, ha="left", color=INK)
    fig.tight_layout(rect=(0, 0, 1, 0.97))
    fig.savefig(OUT / "transport_limit.png", dpi=150)


if __name__ == "__main__":
    OUT.mkdir(parents=True, exist_ok=True)
    fig_box_temperature()
    fig_vchp_front()
    fig_transport_limit()
    print("wrote", ", ".join(p.name for p in sorted(OUT.glob("*.png"))))
