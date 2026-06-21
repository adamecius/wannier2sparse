#!/usr/bin/env python3
"""Plot a Hall-type conductivity sigma(E) produced by tools/hr_exactdiag.py and
make the quantization (or its absence) visible.

Two json producers feed this one plotter:
  * `hr_exactdiag.py ahc <seed>`  -> key "sigma_xy_e2h"  (charge AHC, e^2/h;
        a flat integer plateau in a gap IS the Chern number);
  * `hr_exactdiag.py shc <seed>`  -> key "shc"           (intrinsic spin Hall,
        natural / calibrated units; a trivial insulator shows NO plateau).

The single pedagogical point of the figure is the contrast between a quantized
plateau (Haldane: an exact integer step across the gap) and a non-quantized
response (PdSe2: a large peak but no flat integer in the gap). The shaded band
marks the energy window claimed to be "the gap"; read whether sigma is flat and
integer there (quantized) or not (trivial).

Style: PRL/APS single-column per agent_template/plotting_style.md. LaTeX is used
when available, else Computer-Modern mathtext (the path actually used is printed).
No title; the caption travels with the figure (printed to stdout and a sidecar
.caption.txt). Units and every visual encoding are defined in that caption.

Usage:
    python3 plot_hall.py JSON --out fig.png [--ef E_F] [--ylabel LBL]
        [--gap LO HI] [--plateau] [--ymax Y]
"""
import sys, os, json, argparse, shutil
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

OKABE_ITO = ["#0072B2", "#D55E00", "#009E73", "#CC79A7", "#E69F00", "#56B4E9", "#000000"]


def style(usetex):
    """Reproduce the agent_template/plotting_style.md single-column preset exactly
    (no prl_style.py is shipped). Width is locked at export; only aspect is free."""
    matplotlib.rcParams.update({
        "text.usetex": usetex,
        "font.family": "serif",
        "font.serif": ["cmr10", "Computer Modern Roman", "DejaVu Serif"],
        "mathtext.fontset": "cm",
        "axes.unicode_minus": False,
        "font.size": 9,
        "axes.labelsize": 10, "axes.titlesize": 10,
        "xtick.labelsize": 8, "ytick.labelsize": 8,
        "legend.fontsize": 8, "legend.frameon": False,
        "lines.linewidth": 1.4, "lines.markersize": 5,
        "axes.linewidth": 0.6,
        "xtick.direction": "in", "ytick.direction": "in",
        "xtick.top": True, "ytick.right": True,
        "xtick.minor.visible": True, "ytick.minor.visible": True,
        "figure.dpi": 600, "savefig.dpi": 600,
    })


def main():
    p = argparse.ArgumentParser()
    p.add_argument("json")
    p.add_argument("--out", default="hall.png")
    p.add_argument("--ef", type=float, default=0.0, help="shift x by -E_F (eV)")
    p.add_argument("--ylabel", default=r"$\sigma_{xy}\;[e^2/h]$")
    p.add_argument("--gap", nargs=2, type=float, default=None, metavar=("LO", "HI"),
                   help="shade [LO,HI] (post-shift eV); the window read for a plateau")
    p.add_argument("--gap2", nargs=2, type=float, default=None, metavar=("LO", "HI"),
                   help="shade a second window (e.g. a topological gap), distinct fill")
    p.add_argument("--plateau", action="store_true",
                   help="annotate the mean value in the shaded window (Chern read-off)")
    p.add_argument("--ymax", type=float, default=None)
    p.add_argument("--xlim", nargs=2, type=float, default=None)
    a = p.parse_args()

    d = json.load(open(a.json))
    E = np.array(d["energy_eV"]) - a.ef
    key = next(k for k in ("sigma_xy_e2h", "shc", "sigma_shc_raw") if k in d)
    s = np.array(d[key])

    usetex = bool(shutil.which("latex") and shutil.which("dvipng"))
    style(usetex)
    print(f"[plot_hall] LaTeX path: {'usetex' if usetex else 'mathtext/cm'}")

    fig, ax = plt.subplots(figsize=(3.375, 2.35))
    me = max(1, len(E) // 22)
    if a.gap is not None:
        ax.axvspan(a.gap[0], a.gap[1], color="0.88", lw=0, zorder=0)       # trivial gap: plain gray
    if a.gap2 is not None:
        ax.axvspan(a.gap2[0], a.gap2[1], facecolor=OKABE_ITO[4], alpha=0.35,
                   hatch="///", edgecolor=OKABE_ITO[1], lw=0.0, zorder=0)   # topological gap: hatched
    ax.axhline(0.0, color="0.6", lw=0.6, zorder=1)
    ax.plot(E, s, color=OKABE_ITO[0], ls="-", marker="o", markevery=me,
            markerfacecolor="white", markeredgecolor=OKABE_ITO[0],
            markeredgewidth=1.0, zorder=3)

    if a.plateau and a.gap is not None:
        m = (E > a.gap[0]) & (E < a.gap[1])
        mu = s[m].mean()
        ax.axhline(mu, color=OKABE_ITO[1], ls="--", lw=1.0, zorder=2)
        ax.annotate(rf"$\langle\sigma\rangle = {mu:+.2f}\,e^2/h$",
                    xy=(0.5 * (a.gap[0] + a.gap[1]), mu),
                    xytext=(0.05, 0.62), textcoords="axes fraction",
                    color=OKABE_ITO[1], fontsize=8.5)

    ax.set_xlabel(r"$E - E_F\;[\mathrm{eV}]$" if a.ef else r"$E\;[\mathrm{eV}]$")
    ax.set_ylabel(a.ylabel)
    if a.xlim:
        ax.set_xlim(*a.xlim)
    else:
        ax.set_xlim(E.min(), E.max())
    if a.ymax is not None:
        ax.set_ylim(-a.ymax, a.ymax)
    fig.savefig(a.out, bbox_inches="tight", pad_inches=0.02)
    print(f"[plot_hall] wrote {a.out}")


if __name__ == "__main__":
    main()
