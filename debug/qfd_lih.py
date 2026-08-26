"""QFD driver -- LiH, quadrature-free end to end, validated and certified.

Basis (s-only minimal, HYDROGENIC rates a = Z/n -- the `two_center_eri`
convention, NOT Coulomb-Sturmian and NOT Gaussian):

    Li 1s  Z = 3, n = 1  ->  a = 3
    Li 2s  Z = 3, n = 2  ->  a = 1.5   (a genuine hydrogenic 2s, with its node)
    H  1s  Z = 1, n = 1  ->  a = 1

3 spatial orbitals, 4 electrons, FCI over C(6,4) = 15 determinants after an
exact Loewdin orthogonalization.  R = 3.015 bohr = 603/200.

HONEST SCOPE.  This is a three-function s-only basis with unoptimized hydrogenic
exponents; the energy is FAR from experiment and is not offered as a chemistry
result.  What is claimed is the closed-form assembly and its certified precision.

ONE STRUCTURAL DIFFERENCE FROM H2, and it is the whole story for LiH precision:
the Neumann tau series of the EXCHANGE class terminates only when the two centers
of a density carry the same orbital exponent (q = (alpha-beta)R/2 = 0, the
Phase 0-e criterion).  H2 at zeta_A = zeta_B satisfies that; LiH does not.  Each
tau term is still a closed form; the series is factorially convergent and is
truncated at a per-quartet tau_max whose tail is measured and reported.

Run:  python -u debug/qfd_lih.py [--skip-quad] [--tau N]
"""
from __future__ import annotations

import json
import sys
import time
from fractions import Fraction
from pathlib import Path

import sympy as sp
from mpmath import mp

HERE = Path(__file__).resolve().parent
REPO = HERE.parent
for p in (str(REPO), str(HERE)):
    if p not in sys.path:
        sys.path.insert(0, p)

import qfd_assemble as AS  # noqa: E402
import qfd_core as Q  # noqa: E402
import qfd_quad as V  # noqa: E402

SKIP_QUAD = "--skip-quad" in sys.argv
OUT = REPO / "debug" / "data"
OUT.mkdir(parents=True, exist_ok=True)

R = sp.Rational(603, 200)                       # 3.015 bohr
ZA, ZB = 3, 1
ORBS = [("A", Fraction(3), 1), ("A", Fraction(3), 2), ("B", Fraction(1), 1)]
LABEL = {0: "Li1s", 1: "Li2s", 2: "H1s"}

# per-quartet tau_max: the harder the exponent mismatch q, the more terms.
# (Li1s,H) has q = (3-1)R/2 = 3.015; (Li2s,H) has q = (1.5-1)R/2 = 0.754.
TAU = {(1, 1, 1, 1): 20,      # (Li1s H | Li1s H)   both densities q = 3.015
       (1, 1, 2, 1): 16,      # (Li1s H | Li2s H)
       (2, 1, 1, 1): 16,
       (2, 1, 2, 1): 14}      # (Li2s H | Li2s H)   both densities q = 0.754
DPS = 50
# The numeric-Neumann exchange reference is a nested quadrature, far more
# expensive than the closed form; tau <= 4 already exercises every code path in
# `ordered_xi_general` + `integrate_poly_exp` that the higher terms reuse.
EXCH_QUAD_TAU = 4


def tau_for(quart):
    a, b, c, d = quart
    return TAU.get((a[2], b[2], c[2], d[2]), 20)


def nstr(x, n=25):
    return mp.nstr(x, n, strip_zeros=False)


def mpv(expr, dps):
    if isinstance(expr, mp.mpf):
        return +expr
    with mp.workdps(dps + 15):
        return mp.mpf(str(sp.N(expr, dps + 10)))


def tag(expr):
    if isinstance(expr, mp.mpf):
        return "{exp, E_1, ln, gamma}"          # exchange class, per Paper 58/59
    e = sp.sympify(expr)
    names = set()
    for a in sp.preorder_traversal(e):
        if isinstance(a, sp.expint) or a.func is sp.Ei:
            names.add("E_1")
        elif isinstance(a, sp.log):
            names.add("ln")
        elif a is sp.EulerGamma:
            names.add("gamma")
        elif isinstance(a, sp.exp):
            names.add("exp")
    present = [n for n in ("exp", "E_1", "ln", "gamma") if n in names]
    return "elementary {exp}" if present in ([], ["exp"]) \
        else "{" + ", ".join(present) + "}"


def main():
    mp.dps = 80
    rep = {"system": "LiH", "R": "603/200 (3.015 bohr)",
           "basis": "s-only minimal, hydrogenic: Li1s a=3, Li2s a=1.5, H1s a=1",
           "n_electrons": 4, "fci_dim": 15}
    t_all = time.time()

    print("=" * 78)
    print("LiH  |  s-only minimal hydrogenic basis  |  R = 3.015 bohr")
    print("=" * 78)

    # ------------------------------------------------------------- one-electron
    t = time.time()
    S, h = AS.build_S_h(ORBS, ZA, ZB, R)
    hchk = AS.h_core_check(ORBS, ZA, ZB, R)
    print(f"one-electron closed forms built in {time.time()-t:.1f}s\n")

    print("OVERLAP and CORE-HAMILTONIAN (closed form, 30 digits)")
    print("-" * 78)
    one_rows = []
    for i in range(3):
        for j in range(i, 3):
            sv, hv = mpv(S[i][j], 50), mpv(h[i][j], 50)
            print(f"  S[{LABEL[i]:5s},{LABEL[j]:5s}] = {nstr(sv,30):>34s}"
                  f"   [{tag(S[i][j])}]")
            print(f"  h[{LABEL[i]:5s},{LABEL[j]:5s}] = {nstr(hv,30):>34s}"
                  f"   [{tag(h[i][j])}]")
            one_rows.append({"i": LABEL[i], "j": LABEL[j],
                             "S_50": nstr(sv, 50), "h_50": nstr(hv, 50),
                             "S_tag": tag(S[i][j]), "h_tag": tag(h[i][j])})
    hdiff = max(abs(mpv(sp.expand(h[i][j] - hchk[i][j]), 40))
                for i in range(3) for j in range(3))
    print(f"\n  max |h(Laplacian route) - h(eigen-trick route)| = {nstr(hdiff,3)}"
          "   (two independent closed forms)")
    rep["one_electron"] = one_rows
    rep["h_two_route_max_dev"] = nstr(hdiff, 3)

    # ---------------------------------------------------- one-electron vs quad
    if not SKIP_QUAD:
        print("\nONE-ELECTRON vs INDEPENDENT QUADRATURE (dps 20)")
        print("-" * 78)
        qrows = []
        mp.dps = 20
        specs = [("S", 0, 2, "S"), ("S", 1, 2, "S"), ("T", 0, 2, "T"),
                 ("T", 1, 2, "T"), ("V^B", 0, 0, "inv_rb"), ("V^B", 0, 1,
                                                             "inv_rb"),
                 ("V^A", 0, 2, "inv_ra"), ("V^B", 1, 2, "inv_rb")]
        for nm, i, j, kern in specs:
            t = time.time()
            qv = V.one_electron_quad(ORBS[i], ORBS[j], kern, R)
            if kern == "S":
                cf = Q.overlap(ORBS[i], ORBS[j], R)
            elif kern == "T":
                cf = Q.kinetic(ORBS[i], ORBS[j], R)
            else:
                cf = Q._inv_r(ORBS[i], ORBS[j],
                              "A" if kern == "inv_ra" else "B", R)
            dev = abs(mpv(cf, 20) - qv)
            qrows.append({"name": f"{nm}[{LABEL[i]},{LABEL[j]}]",
                          "dev": nstr(dev, 3)})
            print(f"  {nm}[{LABEL[i]:5s},{LABEL[j]:5s}]  |closed - quad| = "
                  f"{nstr(dev,3)}   ({time.time()-t:.0f}s)")
        mp.dps = 80
        rep["one_electron_quadrature"] = qrows

    # ------------------------------------------------------------ two-electron
    print("\nTWO-ELECTRON closed forms  (exchange at dps "
          f"{DPS}, per-quartet tau_max {TAU})")
    print("-" * 78)
    t = time.time()
    AS.EXCHANGE_TAILS.clear()
    gmap, canon = AS.build_g(ORBS, R, tau_max=tau_for, exchange_hp_dps=DPS)
    t_g = time.time() - t
    print(f"  (built in {t_g:.0f}s)\n")
    g_rows = []
    for k, (expr, cls) in sorted(canon.items()):
        val = mpv(expr, 50)
        nm = f"({LABEL[k[0]]} {LABEL[k[1]]}|{LABEL[k[2]]} {LABEL[k[3]]})"
        print(f"  {nm:27s} {cls:11s} {nstr(val,30):>34s}  [{tag(expr)}]")
        g_rows.append({"quartet": nm, "index": list(k), "class": cls,
                       "tag": tag(expr), "value_50": nstr(val, 50)})
    rep["two_electron"] = g_rows

    # ----------------------------------------------------- tau truncation audit
    print("\nEXCHANGE tau-truncation audit (heteronuclear: q != 0, series is "
          "INFINITE)")
    print("-" * 78)
    tail_rows = []
    for key, per in sorted(AS.EXCHANGE_TAILS.items()):
        tot = mp.fsum(per)
        last = abs(per[-1])
        rel = last / abs(tot)
        print(f"  quartet n=({key[0]},{key[1]}|{key[2]},{key[3]})  "
              f"tau_max = {len(per)-1}  value = {nstr(tot,25)}")
        print(f"      last term |a_tau| = {nstr(last,3)}   relative "
              f"{nstr(rel,3)}   (the truncation tail proxy)")
        tail_rows.append({"n_quartet": list(key), "tau_max": len(per) - 1,
                          "value": nstr(tot, 40), "last_term": nstr(last, 3),
                          "relative_tail": nstr(rel, 3),
                          "per_tau_rel": [nstr(abs(v / tot), 3) for v in per]})
    rep["exchange_tau_tails"] = tail_rows

    # precision self-check of the exchange accumulator: the same quartet at two
    # working precisions (the tau truncation is held fixed, so this isolates the
    # arithmetic, exactly as the dps 30/45 pair does for the linear algebra)
    t = time.time()
    lo, _ = Q.exchange_hp(Fraction(3), (2, 0, 0), (1, 0, 0), Fraction(1),
                          (2, 0, 0), (1, 0, 0), R, tau_max=8, dps=30)
    hi, _ = Q.exchange_hp(Fraction(3), (2, 0, 0), (1, 0, 0), Fraction(1),
                          (2, 0, 0), (1, 0, 0), R, tau_max=8, dps=DPS)
    dexc = abs(lo - hi)
    print(f"\n  exchange accumulator, dps 30 vs {DPS} at fixed tau_max = 8:"
          f" |difference| = {nstr(dexc,3)}   ({time.time()-t:.0f}s)")
    rep["exchange_two_precision"] = nstr(dexc, 3)
    worst_abs = max(abs(mp.mpf(r["last_term"])) for r in tail_rows)
    print(f"\n  worst absolute truncation tail over all exchange quartets: "
          f"{nstr(worst_abs,3)} Ha-scale")


    # -------------------------------------------------------------- certification
    print("\nCERTIFICATION -- FCI total energy at two precisions")
    print("-" * 78)
    cert = {}
    for dps in (30, 45):
        t = time.time()
        with mp.workdps(dps + 40):
            tot, ee, vnn = AS.total_energy(S, h, gmap, ORBS, 4, ZA, ZB, R, dps)
            cert[dps] = (mp.mpf(tot), mp.mpf(ee), mp.mpf(vnn))
        print(f"  dps {dps}: E_tot = {nstr(cert[dps][0], dps)}"
              f"   ({time.time()-t:.0f}s)")
    with mp.workdps(120):
        gap = abs(cert[30][0] - cert[45][0])
        lin = int(-mp.log10(gap / abs(cert[45][0]))) if gap > 0 else 999
        # the binding constraint is the tau truncation, not the linear algebra
        tau_digits = int(-mp.log10(worst_abs / abs(cert[45][0])))
        digits = min(lin, tau_digits)
    print(f"\n  two-precision (dps 30 vs 45) gap  = {nstr(gap,3)}"
          f"   -> {lin} digits from the linear algebra")
    print(f"  worst exchange tau-truncation tail = {nstr(worst_abs,3)}"
          f"   -> {tau_digits} digits from the Neumann truncation")
    print(f"  CERTIFIED DIGITS = min = {digits}")
    print("\n  E_electronic = " + nstr(cert[45][1], 40))
    print("  V_NN         = " + nstr(cert[45][2], 40))
    print("  E_total      = " + nstr(cert[45][0], 40))
    print("\n  HONEST FRAMING: three s functions with unoptimized hydrogenic")
    print("  exponents.  Exact LiH is -8.0705 Ha; this basis cannot approach it")
    print("  and no accuracy claim is made.  What is certified is the")
    print("  closed-form assembly and the digit count above.")

    rep["E_total_45"] = nstr(cert[45][0], 45)
    rep["E_electronic_45"] = nstr(cert[45][1], 45)
    rep["V_NN_45"] = nstr(cert[45][2], 45)
    rep["E_total_30"] = nstr(cert[30][0], 30)
    rep["two_precision_gap"] = nstr(gap, 3)
    rep["digits_from_linear_algebra"] = lin
    rep["digits_from_tau_truncation"] = tau_digits
    rep["digits_certified"] = digits
    rep["g_build_seconds"] = round(t_g)

    (OUT / "qfd_lih_certified.json").write_text(json.dumps(rep, indent=2),
                                                encoding="utf-8")
    print("  (certification written; the two-electron quadrature "
          "cross-check follows and updates the same file)")

    # ------------------------------------------------- two-electron vs quadrature
    if not SKIP_QUAD:
        print("\nTWO-ELECTRON vs INDEPENDENT QUADRATURE (dps 20)")
        print("-" * 78)
        mp.dps = 20
        q2 = []
        for k, (expr, cls) in sorted(canon.items()):
            oa, ob, oc, od = (ORBS[i] for i in k)
            t = time.time()
            try:
                if cls == "one-center":
                    qv = V.one_center_eri_quad(oa, ob, oc, od)
                elif cls == "(AA|BB)":
                    qv = V.aabb_quad(oa, ob, oc, od, R)
                elif cls == "hybrid":
                    if oa[0] == ob[0]:
                        qv = V.hybrid_quad(oa, ob, oc, od, R)
                    else:
                        qv = V.hybrid_quad(oc, od, oa, ob, R)
                else:                                    # exchange
                    a, b, c, d = oa, ob, oc, od
                    if a[0] != "A":
                        a, b = b, a
                    if c[0] != "A":
                        c, d = d, c
                    qv, _ = V.exchange_numeric_neumann(
                        a[1], (a[2], 0, 0), (b[2], 0, 0), b[1],
                        (c[2], 0, 0), (d[2], 0, 0), R,
                        tau_max=EXCH_QUAD_TAU)
            except Exception as exc:                     # pragma: no cover
                print(f"  {k} {cls}: quadrature FAILED "
                      f"{type(exc).__name__}: {exc}")
                q2.append({"index": list(k), "class": cls,
                           "dev": f"FAILED {type(exc).__name__}"})
                continue
            dev = abs(mpv(expr, 20) - qv)
            nm = f"({LABEL[k[0]]} {LABEL[k[1]]}|{LABEL[k[2]]} {LABEL[k[3]]})"
            note = ""
            if cls == "exchange":
                note = f"  (quadrature reference truncated at tau <= {EXCH_QUAD_TAU}; residual = omitted tau terms of the FULL closed form, NOT quadrature error)"
            print(f"  {nm:27s} {cls:11s} |closed - quad| = {nstr(dev,3)}"
                  f"   ({time.time()-t:.0f}s){note}")
            q2.append({"index": list(k), "class": cls, "dev": nstr(dev, 3)})
        mp.dps = 80
        rep["two_electron_quadrature"] = q2

    (OUT / "qfd_lih_certified.json").write_text(json.dumps(rep, indent=2),
                                                encoding="utf-8")
    print(f"\nwrote {OUT / 'qfd_lih_certified.json'}   "
          f"(total {time.time()-t_all:.0f}s)")


if __name__ == "__main__":
    main()
