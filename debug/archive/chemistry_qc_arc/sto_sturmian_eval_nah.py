"""Q4 of the STO -> shared-k Sturmian evaluator scoping: is the Gaussian-fit error material in
any Paper 58 NaH number? (debug only).

Same construction as tests/test_paper58_nah_ladder.py (Na{1s,2s,2p,3s} + H{1s}, atom-only
zetas, 12 electrons), evaluated with
  prod : the paper's shapes -- STO-6G 1s (H and Na), STO-3G 2s, 6-Gaussian LSQ 2p / 3s;
  L<n> : engine LSQ fits (fit_sto_shape) of ALL four shapes at n Gaussians, n = 6..12,
so the L-sequence converges onto the exact-STO value (H2/LiH: debug/sto_sturmian_eval_headtohead.py
shows MD over these fits converges monotonically onto the exact closed form, error ~ in h).

For each: E(Na), E(H), the ladder on the coarse 0.25-a0 grid, R_eq by the paper's 3-point
parabola through the grid minimum, D_e against the SAME-shape atoms, %FCI.
Writes debug/data/sto_sturmian_eval_nah.json
"""
from __future__ import annotations

import json
import os
import time

import numpy as np

from geovac import noci_engine as E

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "data", "sto_sturmian_eval_nah.json")
HA_EV = 27.211386245988
UP, DN = 0, 1
NA1S, NA2S, NA2PX, NA2PY, NA2PZ, NA3S, H1S = range(7)
ZETA_NA = {"1s": 10.63, "2s": 3.3, "2p": 3.44, "3s": 0.836}
GRID = np.arange(3.25, 4.51, 0.25)


def shapes(tag):
    if tag == "prod":
        sh = {"1s": E.STO6G_1S, "2s": E.STO3G_2S}
        for kind, (l, nr) in (("2p", (1, 2)), ("3s", (0, 3))):
            a, d, _ = E.fit_sto_shape(l, nr)
            sh[kind] = (a, d)
        return sh
    ng = int(tag[1:])
    sh = {}
    for kind, (l, nr) in (("1s", (0, 1)), ("2s", (0, 2)), ("2p", (1, 2)), ("3s", (0, 3))):
        a, d, _ = E.fit_sto_shape(l, nr, n_gauss=ng)
        sh[kind] = (a, d)
    return sh


def na_orbs(c, sh):
    z = ZETA_NA
    o = [E.BasisFn(c, (0, 0, 0), sh["1s"][0] * z["1s"] ** 2, sh["1s"][1]),
         E.BasisFn(c, (0, 0, 0), sh["2s"][0] * z["2s"] ** 2, sh["2s"][1])]
    o += [E.BasisFn(c, lmn, sh["2p"][0] * z["2p"] ** 2, sh["2p"][1])
          for lmn in ((1, 0, 0), (0, 1, 0), (0, 0, 1))]
    o.append(E.BasisFn(c, (0, 0, 0), sh["3s"][0] * z["3s"] ** 2, sh["3s"][1]))
    return o


def ladders():
    core = [(o, sp) for o in (NA1S, NA2S, NA2PX, NA2PY, NA2PZ) for sp in (UP, DN)]
    a = core + [(NA3S, UP), (H1S, DN)]
    b = core + [(H1S, UP), (NA3S, DN)]
    ih = core + [(H1S, UP), (H1S, DN)]
    ina = core + [(NA3S, UP), (NA3S, DN)]
    return {"cov2": [a, b], "noci3": [a, b, ih], "all4": [a, b, ih, ina]}


def point(R, sh):
    o0, hp = np.zeros(3), np.array([0.0, 0.0, R])
    orbs = na_orbs(o0, sh) + [E.BasisFn(hp, (0, 0, 0), sh["1s"][0], sh["1s"][1])]
    s, h, g = E.integral_set_md(orbs, [(o0, 11.0), (hp, 1.0)])
    vnn = 11.0 / R
    out = {}
    for lab, d in ladders().items():
        out[lab] = E.noci_ground_gensc(d, s, h, g)[0] + vnn
    ht, gt = E.transform_integrals(E.lowdin_orbitals(s), h, g)
    out["fci"] = E.fci_ground(ht, gt, 12) + vnn
    sid = np.eye(7)
    out["low3"] = E.noci_ground_gensc(ladders()["noci3"], sid, ht, gt)[0] + vnn
    out["S_3s_H"] = float(s[NA3S, H1S])
    return out


def parab(grid, y):
    i = int(np.argmin(y))
    if i == 0 or i == len(y) - 1:
        return float("nan"), float("nan")
    x0, x1, x2 = grid[i - 1], grid[i], grid[i + 1]
    y0, y1, y2 = y[i - 1], y[i], y[i + 1]
    den = (x0 - x1) * (x0 - x2) * (x1 - x2)
    a = (x2 * (y1 - y0) + x1 * (y0 - y2) + x0 * (y2 - y1)) / den
    b = (x2 * x2 * (y0 - y1) + x1 * x1 * (y2 - y0) + x0 * x0 * (y1 - y2)) / den
    c = y1 - a * x1 * x1 - b * x1
    r = -b / (2 * a)
    return float(r), float(a * r * r + b * r + c)


def main() -> None:
    res = {}
    for tag in ("prod", "L6", "L8", "L10", "L12"):
        t0 = time.time()
        sh = shapes(tag)
        o0 = np.zeros(3)
        s, h, g = E.integral_set_md(na_orbs(o0, sh), [(o0, 11.0)])
        ht, gt = E.transform_integrals(E.lowdin_orbitals(s), h, g)
        e_na = float(E.fci_ground(ht, gt, 11))
        chi = E.BasisFn(o0, (0, 0, 0), sh["1s"][0], sh["1s"][1])
        s1, h1, _ = E.integral_set_md([chi], [(o0, 1.0)])
        e_h = float(h1[0, 0] / s1[0, 0])
        pts = [point(float(R), sh) for R in GRID]
        e_diss = e_na + e_h
        row = {"E_Na": e_na, "E_H": e_h, "R3.5": pts[1]}
        dfci = None
        meth = {}
        for lab in ("cov2", "noci3", "all4", "fci", "low3"):
            y = np.array([p[lab] for p in pts])
            req, emin = parab(GRID, y)
            meth[lab] = {"R_eq": req, "D_e_eV": (e_diss - emin) * HA_EV, "E_min": emin}
        for lab in meth:
            meth[lab]["pct_FCI"] = 100 * meth[lab]["D_e_eV"] / meth["fci"]["D_e_eV"]
        # fixed-R (3.5) binding fractions, the test's form
        b_fci = e_diss - pts[1]["fci"]
        row["fixedR_pct"] = {lab: 100 * (e_diss - pts[1][lab]) / b_fci for lab in ("cov2", "noci3", "all4", "low3")}
        row["methods"] = meth
        row["sec"] = time.time() - t0
        res[tag] = row
        print(tag, f"E_Na={e_na:.6f} E_H={e_h:.7f}", {k: (round(v["R_eq"], 4), round(v["D_e_eV"], 4),
                                                       round(v["pct_FCI"], 2)) for k, v in meth.items()},
              {k: round(v, 2) for k, v in row["fixedR_pct"].items()}, f"{row['sec']:.0f}s", flush=True)
    with open(OUT, "w") as fh:
        json.dump(res, fh, indent=1)
    print("wrote", OUT)


if __name__ == "__main__":
    main()
