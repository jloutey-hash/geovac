"""QG1: the 3-D one-loop triangle with masses (gamma, gamma, 0) is elementary.
Closed form (solid angle, debug/ring_ext_QG_common.py) against direct 3-D integration of
int d^3 l prod_i (|l - x_i|^2 + m_i^2)^(-nu_i) at five random external triples, for
(nu) = (1,1,1) and (2,2,1) [the geminal loop], gamma = 1 and 3/2; plus the all-massive case
(gamma, gamma, gamma) used by the gamma -> 0 control.  Direct route: spherical coordinates about x3,
Gauss-Legendre in cos(theta) and r (mapped), trapezoid in phi; float64, independent of the closed form.
Gate: relative difference < 1e-8 at every triple."""
from __future__ import annotations
import json, random, sys
from pathlib import Path
import numpy as np
import mpmath as mp
sys.path.insert(0, str(Path(__file__).parent))
import ring_ext_QG_common as Cm

OUT = Path(__file__).parent / "data" / "ring_ext_QG_1_loop.json"


def direct(x1, x2, m, nu, Nang=160, Nr=240):
    """int d^3 l f, x3 = 0 at the origin of spherical coordinates."""
    x1, x2 = np.array(x1), np.array(x2)
    uu, wu = np.polynomial.legendre.leggauss(Nang)
    ph = np.linspace(0, 2 * np.pi, Nang, endpoint=False); wph = 2 * np.pi / Nang
    st = np.sqrt(1 - uu**2)
    n = np.stack([np.outer(st, np.cos(ph)), np.outer(st, np.sin(ph)), np.outer(uu, np.ones_like(ph))], -1)
    d1, d2 = n @ x1, n @ x2
    # radial: split [0, R1], [R1, inf) with r = R1 / s
    R = max(np.linalg.norm(x1), np.linalg.norm(x2)) + 2
    sr, wr = np.polynomial.legendre.leggauss(Nr)
    tot = 0.0
    for lo, hi in [(0, R / 3), (R / 3, 2 * R / 3), (2 * R / 3, R)]:
        r = (hi - lo) / 2 * sr + (hi + lo) / 2; w = (hi - lo) / 2 * wr
        for ri, wi in zip(r, w):
            f = (ri**2 - 2 * ri * d1 + x1 @ x1 + m[0]**2) ** -nu[0] * (ri**2 - 2 * ri * d2 + x2 @ x2 + m[1]**2) ** -nu[1] \
                * (ri**2 + m[2]**2) ** -nu[2] * ri**2
            tot += wi * np.sum(wu[:, None] * f) * wph
    # tail r = R/s, s in (0,1], dr = R/s^2 ds
    s = (sr + 1) / 2; ws = wr / 2
    for si, wi in zip(s, ws):
        ri = R / si
        f = (ri**2 - 2 * ri * d1 + x1 @ x1 + m[0]**2) ** -nu[0] * (ri**2 - 2 * ri * d2 + x2 @ x2 + m[1]**2) ** -nu[1] \
            * (ri**2 + m[2]**2) ** -nu[2] * ri**2
        tot += wi * R / si**2 * np.sum(wu[:, None] * f) * wph
    return tot


def closed(x1, x2, m, nu):
    x1 = [mp.mpf(v) for v in x1]; x2 = [mp.mpf(v) for v in x2]
    s12 = sum((a - b) ** 2 for a, b in zip(x1, x2)); s13 = sum(a * a for a in x1); s23 = sum(b * b for b in x2)
    if nu == (1, 1, 1):
        return Cm.I111(m[0], m[1], m[2], s12, s13, s23)
    return Cm.I221(m[0], m[1], m[2], s12, s13, s23)


def main():
    mp.mp.dps = 30
    rng = random.Random(20261017)
    rec = {}
    worst = 0
    for g in (1.0, 1.5):
        for masses, nu, lab in [((g, g, 0.0), (1, 1, 1), "111_gg0"), ((g, g, 0.0), (2, 2, 1), "221_gg0"),
                                ((g, g, g), (1, 1, 1), "111_ggg")]:
            rows = []
            for k in range(5):
                x1 = [rng.uniform(-2, 2) for _ in range(3)]; x2 = [rng.uniform(-2, 2) for _ in range(3)]
                a = direct(x1, x2, masses, nu); b = float(closed(x1, x2, masses, nu))
                rel = abs(a / b - 1); worst = max(worst, rel)
                rows.append({"x1": x1, "x2": x2, "direct": a, "closed": b, "rel": rel})
                print(f"gamma {g} {lab} triple {k}: direct {a:.15e} closed {b:.15e} rel {rel:.2e}", flush=True)
            rec[f"g{g}_{lab}"] = rows
    # massless limit of the closed form -> pi^3/(|x12||x13||x23|) (the direct integrator is not built
    # for the integrable 1/|l-x|^2 singularities of three massless lines, so this is a closed-form check)
    x1, x2 = [0.3, -1.1, 0.7], [1.2, 0.4, -0.5]
    s12 = sum((p - q) ** 2 for p, q in zip(x1, x2)); s13 = sum(p * p for p in x1); s23 = sum(q * q for q in x2)
    uniq = mp.pi ** 3 / mp.sqrt(mp.mpf(s12) * s13 * s23)
    lim = Cm.I111(mp.mpf("1e-12"), mp.mpf("1e-12"), mp.mpf("1e-12"), s12, s13, s23)
    print("massless limit closed/unique - 1 =", mp.nstr(lim / uniq - 1, 3))
    rec["massless_limit_closed"] = float(lim / uniq - 1)
    print("WORST rel:", worst, "GATE", "PASS" if worst < 1e-8 else "FAIL")
    OUT.write_text(json.dumps(rec, indent=1))


if __name__ == "__main__":
    main()
