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


def _shade_gaps(ax, a):
    """Shade the trivial (gray) and topological (hatched) gap windows on an axis."""
    if a.gap is not None:
        ax.axvspan(a.gap[0], a.gap[1], color="0.88", lw=0, zorder=0)        # trivial gap: plain gray
    if a.gap2 is not None:
        ax.axvspan(a.gap2[0], a.gap2[1], facecolor=OKABE_ITO[4], alpha=0.35,
                   hatch="///", edgecolor=OKABE_ITO[1], lw=0.0, zorder=0)    # topological gap: hatched


def _draw_sigma(ax, E, s, a, with_xlabel):
    """The conductivity panel: sigma(E) with the gaps shaded and the plateau read off."""
    me = max(1, len(E) // 22)
    _shade_gaps(ax, a)
    ax.axhline(0.0, color="0.6", lw=0.6, zorder=1)
    ax.plot(E, s, color=OKABE_ITO[0], ls="-", marker="o", markevery=me,
            markerfacecolor="white", markeredgecolor=OKABE_ITO[0],
            markeredgewidth=1.0, zorder=3)
    if a.plateau and a.gap2 is not None:        # the quantized plateau sits in the topological gap
        m = (E > a.gap2[0]) & (E < a.gap2[1])
        mu = s[m].mean()
        ax.axhline(mu, color=OKABE_ITO[1], ls="--", lw=1.0, zorder=2)
        ax.annotate(rf"$\langle\sigma\rangle = {mu:+.2f}\,e^2/h$",
                    xy=(0.5 * (a.gap2[0] + a.gap2[1]), mu),
                    xytext=(0.05, 0.62), textcoords="axes fraction",
                    color=OKABE_ITO[1], fontsize=8.5)
    elif a.plateau and a.gap is not None:
        m = (E > a.gap[0]) & (E < a.gap[1])
        mu = s[m].mean()
        ax.axhline(mu, color=OKABE_ITO[1], ls="--", lw=1.0, zorder=2)
        ax.annotate(rf"$\langle\sigma\rangle = {mu:+.2f}\,e^2/h$",
                    xy=(0.5 * (a.gap[0] + a.gap[1]), mu),
                    xytext=(0.05, 0.62), textcoords="axes fraction",
                    color=OKABE_ITO[1], fontsize=8.5)
    ax.set_ylabel(a.ylabel)
    if a.ymax is not None:
        ax.set_ylim(-a.ymax, a.ymax)
    if with_xlabel:
        ax.set_xlabel(r"$E - E_F\;[\mathrm{eV}]$" if a.ef else r"$E\;[\mathrm{eV}]$")


def _draw_dos(ax, E, dos, a, with_xlabel):
    """The density-of-states panel: rho(E), same energy axis and gap shading."""
    _shade_gaps(ax, a)
    ax.plot(E, dos, color=OKABE_ITO[2], ls="-", zorder=3)
    ax.fill_between(E, dos, color=OKABE_ITO[2], alpha=0.18, lw=0, zorder=2)
    ax.set_ylabel(r"$\rho(E)\;[\mathrm{eV}^{-1}]$")
    ax.set_ylim(bottom=0.0)
    if with_xlabel:
        ax.set_xlabel(r"$E - E_F\;[\mathrm{eV}]$" if a.ef else r"$E\;[\mathrm{eV}]$")


def main():
    p = argparse.ArgumentParser()
    p.add_argument("json")
    p.add_argument("--out", default="hall.png")
    p.add_argument("--ef", type=float, default=0.0, help="shift x by -E_F (eV)")
    p.add_argument("--ylabel", default=r"$\sigma_{xy}\;[e^2/h]$")
    p.add_argument("--gap", nargs=2, type=float, default=None, metavar=("LO", "HI"),
                   help="shade a trivial gap window [LO,HI] (post-shift eV)")
    p.add_argument("--gap2", nargs=2, type=float, default=None, metavar=("LO", "HI"),
                   help="shade a topological gap window; the plateau is read here")
    p.add_argument("--plateau", action="store_true",
                   help="annotate the mean value in the topological (else trivial) gap")
    p.add_argument("--ymax", type=float, default=None)
    p.add_argument("--xlim", nargs=2, type=float, default=None)
    p.add_argument("--dos", default=None, metavar="DOS.json",
                   help="add a shared-axis density-of-states panel (hr_exactdiag.py dos output)")
    a = p.parse_args()

    d = json.load(open(a.json))
    E = np.array(d["energy_eV"]) - a.ef
    key = next(k for k in ("sigma_xy_e2h", "shc", "sigma_shc_raw") if k in d)
    s = np.array(d[key])

    usetex = bool(shutil.which("latex") and shutil.which("dvipng"))
    style(usetex)
    print(f"[plot_hall] LaTeX path: {'usetex' if usetex else 'mathtext/cm'}")

    if a.dos is not None:
        dd = json.load(open(a.dos))
        Ed = np.array(dd["energy_eV"]) - a.ef
        rho = np.array(dd["dos"])
        # Two stacked panels sharing the energy axis: (a) sigma, (b) DOS. Single-
        # column width; the conductivity gets the taller panel.
        fig, (ax1, ax2) = plt.subplots(2, 1, sharex=True, figsize=(3.375, 3.7),
                                       gridspec_kw={"height_ratios": [2, 1], "hspace": 0.08})
        _draw_sigma(ax1, E, s, a, with_xlabel=False)
        _draw_dos(ax2, Ed, rho, a, with_xlabel=True)
        ax1.text(0.035, 0.86, "(a)", transform=ax1.transAxes, fontweight="bold")
        ax2.text(0.035, 0.78, "(b)", transform=ax2.transAxes, fontweight="bold")
        axes = (ax1, ax2)
    else:
        fig, ax = plt.subplots(figsize=(3.375, 2.35))
        _draw_sigma(ax, E, s, a, with_xlabel=True)
        axes = (ax,)

    xlim = a.xlim if a.xlim else (E.min(), E.max())
    for ax in axes:
        ax.set_xlim(*xlim)

    fig.savefig(a.out, bbox_inches="tight", pad_inches=0.02)
    cap = a.out.rsplit(".", 1)[0] + ".caption.txt"
    with open(cap, "w") as f:
        f.write(_caption(a))
    print(f"[plot_hall] wrote {a.out} and {cap}")


def _caption(a):
    """A self-contained FIG. caption (plotting_style.md §7); units/encodings defined."""
    head = (a.ylabel or "sigma")
    for tok in ("$", "\\;", "\\,", "\\", "{", "}"):
        head = head.replace(tok, "")
    lines = [
        f"FIG. N. Intrinsic {head} versus energy"
        + (r" $E-E_F$" if a.ef else " $E$") + " (eV), from exact diagonalization.",
    ]
    if a.dos is not None:
        lines.append(
            "(a) the conductivity; solid blue with open circles. "
            "(b) the density of states rho(E) [eV^-1]; solid green, filled. "
            "Both panels share the energy axis.")
    else:
        lines.append("Solid blue with open circles.")
    if a.gap is not None:
        lines.append(f"Gray band: trivial gap [{a.gap[0]:+.2f},{a.gap[1]:+.2f}] eV "
                     "(rho->0, no plateau).")
    if a.gap2 is not None:
        lines.append(f"Hatched band: topological gap [{a.gap2[0]:+.2f},{a.gap2[1]:+.2f}] eV "
                     "(rho->0 and a flat, near-integer plateau).")
    if a.plateau:
        lines.append("Dashed orange: the mean conductivity across the marked gap.")
    lines.append("Exact-diagonalization reference (tools/hr_exactdiag.py); see the "
                 "JSON 'meta' for the k-mesh and broadening used.")
    return " ".join(lines) + "\n"


if __name__ == "__main__":
    main()
