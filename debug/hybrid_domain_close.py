"""Hybrid (AA|AB) domain close: the a_d >= a_c branches of the shell route.

Before: the shell-route hybrid evaluator NaN'd at the rate coincidence
a_c = a_d (Z_A/n_c = Z_B/n_d, e.g. Z_B = Z_A with n_c = n_d) and raised on
a_d > a_c (an Ei reached step 3). After: finite_power_exp has an exact b = 0
branch (ln / powers), e1_moment takes a negative rate (E_1 -> -Ei, principal
value) and c = 0, and log_moment_finite consumes the ln x the coincidence
leaves at step 3.

Checks:
  1. value vs two referees (hybrid_quadrature: (r_B, cos th_B) dblquad;
     and a nested quad in (r_A, r_B) triangle coordinates written here);
  2. function content of each new branch (atoms of the closed form);
  3. continuity: equal-rate value vs a_d = a_c(1 - 10^-k), k = 4, 6, 8;
  4. bit-identity on the a_d < a_c domain against the HEAD module.

Writes debug/data/hybrid_domain_close.json.
"""
from __future__ import annotations

import importlib.util
import json
import subprocess
import sys
import tempfile
import time
from fractions import Fraction
from pathlib import Path

import warnings
import numpy as np
warnings.filterwarnings('ignore')
import sympy as sp
from scipy import integrate

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from geovac.two_center_eri import (R_s, hybrid_closed_form, hybrid_quadrature,  # noqa: E402
                                   hybrid_terms, hybrid_F_radial, radial_poly,
                                   radial_norm, plm_signed, sph_norm_numeric,
                                   r_s)


def triangle_quadrature(ZA, oa, ob, oc, ZB, od, R: float) -> float:
    """Independent referee in (r_A, r_B) coordinates.

    d^3r = (r_A r_B / R) dr_A dr_B dphi on the triangle |r_B - R| <= r_A <= r_B + R;
    cos th_A = (r_A^2 + R^2 - r_B^2)/(2 R r_A), cos th_B = (r_A^2 - r_B^2 - R^2)/(2 R r_B).
    Shares only the Gaunt re-coupling (hybrid_terms) and F = conj(R_c) V_L.
    """
    coeffs_d, ad = radial_poly(ZB, od[0], od[1])
    Nd = radial_norm(ZB, od[0], od[1])
    Rd = sp.lambdify(r_s, sum(Nd * c * r_s ** k for k, c in coeffs_d.items())
                     * sp.exp(-ad * r_s), "numpy")
    total = 0.0
    for coeff, L, Lp, Mp, ld, md in hybrid_terms(ZA, oa, ob, oc, ZB, od):
        Ff = sp.lambdify(r_s, hybrid_F_radial(ZA, oa, ob, oc, L), "numpy")
        nA, nB = sph_norm_numeric(Lp, Mp), sph_norm_numeric(ld, md)

        def inner(rb):
            def f(ra):
                cA = np.clip((ra * ra + R * R - rb * rb) / (2 * R * ra), -1, 1)
                cB = np.clip((ra * ra - rb * rb - R * R) / (2 * R * rb), -1, 1)
                return float(np.real(
                    ra * rb / R * Ff(ra) * nA * plm_signed(Lp, Mp, np.array([cA]))[0]
                    * Rd(rb) * nB * plm_signed(ld, md, np.array([cB]))[0]))
            lo = abs(rb - R)
            if lo < 1e-300:
                lo = 1e-300
            v, _ = integrate.quad(f, lo, rb + R, epsabs=1e-15, epsrel=1e-12,
                                  limit=100)
            return v

        v1, _ = integrate.quad(inner, 0.0, R, epsabs=1e-14, epsrel=1e-12, limit=100)
        v2, _ = integrate.quad(inner, R, 80.0, epsabs=1e-14, epsrel=1e-12, limit=200)
        total += float(sp.re(coeff)) * 2 * np.pi * (v1 + v2)
    return total


def rate(Z, o):
    return Fraction(Z) / o[0]


POINTS = [
    # name, ZA, oa, ob, oc, ZB, od, R
    ("eq  (2p0 2p0|1s 1s_B) Z=3,3",   3, (2, 1, 0), (2, 1, 0), (1, 0, 0), 3, (1, 0, 0), "5/2"),
    ("eq  (2p0 2p0|1s 1s_B) Z=3,3",   3, (2, 1, 0), (2, 1, 0), (1, 0, 0), 3, (1, 0, 0), "7/5"),
    ("eq  (2p1 2p1|2p1 2p1_B) Z=2,2", 2, (2, 1, 1), (2, 1, 1), (2, 1, 1), 2, (2, 1, 1), "5/2"),
    ("eq  (2p0 2p0|2p0 1s_B) Z=2,1",  2, (2, 1, 0), (2, 1, 0), (2, 1, 0), 1, (1, 0, 0), "3"),
    ("eq  (3d0 3d0|1s 1s_B) Z=3,3",   3, (3, 2, 0), (3, 2, 0), (1, 0, 0), 3, (1, 0, 0), "5/2"),
    ("gt  (2p0 2p0|1s 1s_B) Z=3,4",   3, (2, 1, 0), (2, 1, 0), (1, 0, 0), 4, (1, 0, 0), "5/2"),
    ("gt  (2p0 2p0|1s 1s_B) Z=1,2",   1, (2, 1, 0), (2, 1, 0), (1, 0, 0), 2, (1, 0, 0), "5/2"),
    ("gt  (2p0 2p0|2p0 2p0_B) Z=1,3", 1, (2, 1, 0), (2, 1, 0), (2, 1, 0), 3, (2, 1, 0), "4"),
    ("gt  (2p1 2p1|2p1 2p1_B) Z=2,3", 2, (2, 1, 1), (2, 1, 1), (2, 1, 1), 3, (2, 1, 1), "7/5"),
    ("gt  (3d0 3d0|1s 1s_B) Z=2,3",   2, (3, 2, 0), (3, 2, 0), (1, 0, 0), 3, (1, 0, 0), "5/2"),
    # s-type one-centre pair (direct route) on the same domain
    ("eq  s-pair (1s 1s|2p0 2p0_B) Z=2,4", 2, (1, 0, 0), (1, 0, 0), (2, 1, 0), 4, (2, 1, 0), "5/2"),
    ("gt  s-pair (1s 1s|1s 1s_B) Z=1,3",   1, (1, 0, 0), (1, 0, 0), (1, 0, 0), 3, (1, 0, 0), "5/2"),
]


def head_module():
    """Import geovac/two_center_eri.py as of git HEAD, for the bit-identity check."""
    src = subprocess.run(["git", "show", "HEAD:geovac/two_center_eri.py"],
                         cwd=ROOT, capture_output=True, text=True, check=True).stdout
    tmp = Path(tempfile.mkdtemp()) / "tce_head.py"
    tmp.write_text(src, encoding="utf-8")
    spec = importlib.util.spec_from_file_location("tce_head", tmp)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def main() -> None:
    out = {"points": [], "continuity": [], "bit_identity": []}
    for name, ZA, oa, ob, oc, ZB, od, Rstr in POINTS:
        t0 = time.time()
        Rq = sp.Rational(Rstr)
        e = hybrid_closed_form(Fraction(ZA), oa, ob, oc, Fraction(ZB), od)
        val = sp.re(sp.N(e.subs(R_s, Rq), 40))
        funcs = sorted({type(f).__name__ for f in e.atoms(sp.Function)})
        has_logR = any(a.args[0] == R_s for a in e.atoms(sp.log))
        has_gamma = e.has(sp.EulerGamma)
        print('  closed form', round(time.time() - t0, 1), 's', flush=True)
        q1 = hybrid_quadrature(Fraction(ZA), oa, ob, oc, Fraction(ZB), od, float(Rq))
        q2 = triangle_quadrature(Fraction(ZA), oa, ob, oc, Fraction(ZB), od, float(Rq))
        v = float(val)
        row = dict(name=name, R=Rstr, a_c=str(rate(ZA, oc)), a_d=str(rate(ZB, od)),
                   value=str(sp.N(val, 20)), ref_dblquad=q1, ref_triangle=q2,
                   rel_dblquad=abs(v - q1) / abs(v), rel_triangle=abs(v - q2) / abs(v),
                   functions=funcs, log_of_R=has_logR, euler_gamma=has_gamma,
                   seconds=round(time.time() - t0, 1))
        out["points"].append(row)
        print(f"{name:36s} R={Rstr:4s} {v:+.15e}  dbl {row['rel_dblquad']:.1e}"
              f"  tri {row['rel_triangle']:.1e}  {funcs} lnR={has_logR} g={has_gamma}",
              flush=True)

    # continuity at the coincidence: approach a_d -> a_c from both sides
    e_eq = hybrid_closed_form(Fraction(3), (2, 1, 0), (2, 1, 0), (1, 0, 0),
                              Fraction(3), (1, 0, 0))
    v_eq = sp.N(e_eq.subs(R_s, sp.Rational(5, 2)), 40)
    for k in (4, 6, 8):
        for sgn in (-1, 1):
            ZB = Fraction(3) * (1 + Fraction(sgn, 10 ** k))
            e = hybrid_closed_form(Fraction(3), (2, 1, 0), (2, 1, 0), (1, 0, 0),
                                   ZB, (1, 0, 0))
            vv = sp.re(sp.N(e.subs(R_s, sp.Rational(5, 2)), 40))
            d = float(abs(vv - v_eq))
            out["continuity"].append(dict(eps=f"{sgn}e-{k}", diff=d))
            print(f"continuity eps={sgn}e-{k}: |v - v_eq| = {d:.3e}", flush=True)

    # bit-identity on the old domain a_d < a_c
    old = head_module()
    for q in ((3, (2, 1, 0), (2, 1, 0), (1, 0, 0), 1, (1, 0, 0)),
              (3, (2, 1, 1), (2, 1, 1), (2, 1, 1), 2, (2, 1, 1)),
              (3, (3, 2, 0), (3, 2, 0), (1, 0, 0), 1, (1, 0, 0)),
              (3, (1, 0, 0), (1, 0, 0), (3, 2, 1), 2, (2, 1, 1))):
        ZA, oa, ob, oc, ZB, od = q
        a = hybrid_closed_form(Fraction(ZA), oa, ob, oc, Fraction(ZB), od)
        b = old.hybrid_closed_form(Fraction(ZA), oa, ob, oc, Fraction(ZB), od)
        same = sp.srepr(a) == sp.srepr(b)
        out["bit_identity"].append(dict(quartet=str(q), identical=same))
        print(f"bit-identical {q}: {same}", flush=True)

    path = ROOT / "debug" / "data" / "hybrid_domain_close.json"
    path.write_text(json.dumps(out, indent=1), encoding="utf-8")
    print("wrote", path)


if __name__ == "__main__":
    main()
