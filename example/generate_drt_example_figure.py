#!/usr/bin/env python3
"""
Generator obrazku DRT spektra pro doc/DRT_INTUITION.md
======================================================

Vytvori cisty, reprodukovatelny obrazek realneho DRT spektra ze syntetickych dat
(R_inf + dva RC procesy, tau ~ 1 ms a ~ 100 ms) - presne dvoupikovy priklad
pouzity v dokumentu. Vystup se uklada do doc/images/drt_example.png.

Pouziti:
    python3 example/generate_drt_example_figure.py
"""

from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from eis_analysis.drt import calculate_drt

# Vystupni cesta odvozena od umisteni skriptu (nezavisle na cwd)
REPO_ROOT = Path(__file__).resolve().parent.parent
OUT_PATH = REPO_ROOT / "doc" / "images" / "drt_example.png"

# --- Synteticka ground truth: R_inf + (R1||C1) + (R2||C2) ---
R_inf = 10.0
R1, tau1 = 50.0, 1e-3      # rychly proces, tau = 1 ms
R2, tau2 = 80.0, 1e-1      # pomaly proces, tau = 100 ms


def main() -> None:
    f = np.logspace(5, -2, 70)         # 100 kHz ... 10 mHz
    w = 2 * np.pi * f
    Z = R_inf + R1 / (1 + 1j * w * tau1) + R2 / (1 + 1j * w * tau2)

    # maly reprodukovatelny merici sum
    rng = np.random.default_rng(0)
    Z = Z * (1 + 0.003 * (rng.standard_normal(Z.shape)
                          + 1j * rng.standard_normal(Z.shape)))

    res = calculate_drt(f, Z, n_tau=200, auto_lambda=True, peak_method="gmm")
    if not res.success:
        raise SystemExit("DRT vypocet selhal")

    tau, gamma = res.tau, res.gamma

    fig, ax = plt.subplots(figsize=(8, 4.5))
    ax.semilogx(tau, gamma, color="#1f4e8c", lw=2.2, label="DRT gamma(tau)")
    ax.fill_between(tau, 0, gamma, color="#1f4e8c", alpha=0.18)

    # vyznac skutecne casove konstanty
    for tg, lbl in [(tau1, "tau = 1 ms"), (tau2, "tau = 100 ms")]:
        ax.axvline(tg, color="#c0392b", ls="--", lw=1.0, alpha=0.7)
        ax.text(tg, ax.get_ylim()[1] * 0.92, "  " + lbl, color="#c0392b",
                fontsize=9, ha="left", va="top")

    ax.set_xlabel("tau [s]  (relaxacni casova konstanta)")
    ax.set_ylabel("gamma(tau) [Ohm]")
    ax.set_title("DRT spektrum: dva RC procesy  (plocha piku = odpor)")
    ax.grid(True, which="both", alpha=0.25)
    ax.legend(loc="upper left")

    txt = (f"R_inf = {res.R_inf:.1f} Ohm\n"
           f"R_pol = {res.R_pol:.1f} Ohm\n"
           f"lambda = {res.lambda_used:.2e}\n"
           f"piku nalezeno: {len(res.peaks) if res.peaks else 0}")
    ax.text(0.985, 0.95, txt, transform=ax.transAxes, fontsize=8.5,
            ha="right", va="top", family="monospace",
            bbox=dict(boxstyle="round", fc="white", ec="#888", alpha=0.85))

    fig.tight_layout()
    OUT_PATH.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT_PATH, dpi=130)
    print(f"saved {OUT_PATH}")
    print("R_inf=%.2f R_pol=%.2f lambda=%.3e peaks=%s"
          % (res.R_inf, res.R_pol, res.lambda_used,
             len(res.peaks) if res.peaks else 0))


if __name__ == "__main__":
    main()
