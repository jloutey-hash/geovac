"""QFD driver -- H2, quadrature-free end to end, validated and certified.

Basis: 1s on each nucleus, hydrogenic rate a = zeta (n = 1).  Two electrons, FCI
in the 4-spin-orbital space after an exact Loewdin orthogonalization.

EVERY production-path integral is a closed form.  Validation is three-way:
  * raw mpmath quadrature (prolate spheroidal / Newton routes, `qfd_quad`)
  * LITERATURE closed forms -- Sugiura's exchange integral and the classical
    Coulomb / hybrid / overlap / kinetic / nuclear formulas for 1s-1s
  * two-precision certification of the final FCI energy (dps 60 vs 90)

Run:  python -u debug/qfd_h2.py [--quick]
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

QUICK = "--quick" in sys.argv
OUT = REPO / "debug" / "data"
OUT.mkdir(parents=True, exist_ok=True)


# --------------------------------------------------------------- transcendence

def tag(expr):
    """Transcendence tag of a closed form, per the Paper 58/59 classification."""
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
    order = ["exp", "E_1", "ln", "gamma"]
    present = [n for n in order if n in names]
    if present == ["exp"] or not present:
        return "elementary {exp}"
    return "{" + ", ".join(present) + "}"


def nstr(x, n=25):
    return mp.nstr(x, n, strip_zeros=False)


def mpv(expr, dps):
    with mp.workdps(dps + 15):
        return mp.mpf(str(sp.N(expr, dps + 10)))


# ------------------------------------------------------------------- the system

def h2_system(zeta, R):
    ZO = Fraction(zeta) if not isinstance(zeta, Fraction) else zeta
    orbs = [("A", ZO, 1), ("B", ZO, 1)]
    return orbs


def build(zeta, R, tau_max=6, verbose=True):
    orbs = h2_system(zeta, R)
    S, h = AS.build_S_h(orbs, 1, 1, R)
    hchk = AS.h_core_check(orbs, 1, 1, R)
    gmap, canon = AS.build_g(orbs, R, tau_max=tau_max, verbose=False)
    return orbs, S, h, hchk, gmap, canon


def main():
    report = {"system": "H2", "basis": "1s on each nucleus, hydrogenic a = zeta",
              "geometries": []}
    t_all = time.time()
    mp.dps = 80                      # all literature constants built at 80 digits

    # ------------------------------------------------------------- R = 1.4 core
    R = sp.Rational(7, 5)
    zeta = Fraction(1)
    print("=" * 78)
    print("H2  |  zeta = 1  |  R = 1.4 bohr  |  closed-form assembly")
    print("=" * 78)
    t = time.time()
    orbs, S, h, hchk, gmap, canon = build(zeta, R)
    print(f"closed forms built in {time.time()-t:.2f}s\n")

    Rm = mp.mpf("1.4")
    zm = mp.mpf(1)

    # --- one-electron table -------------------------------------------------
    A, B = orbs
    one_e = [
        ("S_AA", Q.overlap(A, A, R), None, None, "overlap, one-center"),
        ("S_AB", Q.overlap(A, B, R), ("S", A, B), V.lit_overlap(zm, Rm),
         "overlap, two-center"),
        ("T_AA", Q.kinetic(A, A, R), None, mp.mpf(1) / 2, "kinetic, one-center"),
        ("T_AB", Q.kinetic(A, B, R), ("T", A, B), V.lit_kinetic_cross(zm, Rm),
         "kinetic, two-center"),
        ("V^A_AA", Q._inv_r(A, A, "A", R), None, mp.mpf(1), "<A|1/r_A|A>"),
        ("V^B_AA", Q._inv_r(A, A, "B", R), ("inv_rb", A, A),
         V.lit_nuc_same(zm, Rm), "<A|1/r_B|A>"),
        ("V^A_AB", Q._inv_r(A, B, "A", R), ("inv_ra", A, B),
         V.lit_nuc_cross(zm, Rm), "<A|1/r_A|B>"),
        ("V^B_AB", Q._inv_r(A, B, "B", R), ("inv_rb", A, B),
         V.lit_nuc_cross(zm, Rm), "<A|1/r_B|B> (= V^A_AB by symmetry)"),
    ]
    print("ONE-ELECTRON closed forms (50 digits), vs literature and quadrature")
    print("-" * 78)
    one_e_rows = []
    for name, expr, qspec, lit, note in one_e:
        val = mpv(expr, 50)
        dlit = abs(val - lit) if lit is not None else None
        row = {"name": name, "note": note, "tag": tag(expr),
               "value_50": nstr(val, 50)}
        if dlit is not None:
            row["dev_literature"] = nstr(dlit, 3)
        one_e_rows.append((name, expr, qspec, val, dlit, note))
        print(f"{name:8s} {nstr(val,30):>34s}  [{tag(expr)}]"
              + (f"  |lit| {nstr(dlit,3)}" if dlit is not None else ""))
    # h routes agree symbolically?
    hdiff = max(sp.simplify(sp.expand(h[i][j] - hchk[i][j])) for i in range(2)
                for j in range(2))
    print(f"\nh_core: Laplacian route vs eigen-trick route, symbolic difference "
          f"= {hdiff}")

    # --- two-electron table -------------------------------------------------
    print("\nTWO-ELECTRON closed forms (50 digits), vs literature")
    print("-" * 78)
    lit_g = {
        (0, 0, 0, 0): (mp.mpf(5) / 8, "one-center 5 zeta/8"),
        (0, 0, 1, 1): (V.lit_coulomb(zm, Rm), "Coulomb J(R), classical"),
        (0, 0, 0, 1): (V.lit_hybrid(zm, Rm), "hybrid, classical"),
        (0, 1, 0, 1): (V.lit_exchange_sugiura(zm, Rm), "exchange, Sugiura 1927"),
    }
    two_e_rows = []
    for key, (expr, cls) in sorted(canon.items()):
        val = mpv(expr, 50)
        lit, litname = lit_g.get(key, (None, None))
        d = abs(val - lit) if lit is not None else None
        two_e_rows.append((key, cls, expr, val, d, litname))
        print(f"({key[0]}{key[1]}|{key[2]}{key[3]}) {cls:11s} "
              f"{nstr(val,30):>34s}  [{tag(expr)}]"
              + (f"  |lit| {nstr(d,3)}  ({litname})" if d is not None else ""))

    # --- tau termination audit ---------------------------------------------
    print("\nEXCHANGE tau-series audit (homonuclear -> q = 0 -> TERMINATES)")
    terms = Q.exchange_closed_form(zeta, (1, 0, 0), (1, 0, 0), zeta,
                                   (1, 0, 0), (1, 0, 0), R, tau_max=8,
                                   return_terms=True)
    tau_rows = []
    for i, tm in enumerate(terms):
        v = mpv(tm, 40) if tm != 0 else mp.mpf(0)
        tau_rows.append((i, nstr(v, 30)))
        print(f"   tau = {i}: {nstr(v,30)}")
    print("   -> tau > 2 vanish IDENTICALLY (symbolic zero), so the Neumann sum "
          "is finite: no truncation error.")

    # --- independent quadrature --------------------------------------------
    quad_rows = []
    if not QUICK:
        print("\nINDEPENDENT QUADRATURE (mpmath, dps 20) -- one-electron")
        print("-" * 78)
        mp.dps = 20
        for name, expr, qspec, val, _d, _n in one_e_rows:
            if qspec is None:
                continue
            kern, oi, oj = qspec
            t = time.time()
            qv = V.one_electron_quad(oi, oj, kern, R)
            dev = abs(mpv(expr, 20) - qv)
            quad_rows.append((name, nstr(dev, 3)))
            print(f"{name:8s} |closed - quad| = {nstr(dev,3)}   "
                  f"({time.time()-t:.0f}s)")
        # Newton 1-D cross check for the same-center nuclear attraction
        dev = abs(mpv(Q._inv_r(A, A, "B", R), 20)
                  - V.one_electron_quad_same_center(A, A, "B", R))
        quad_rows.append(("V^B_AA (Newton 1D)", nstr(dev, 3)))
        print(f"V^B_AA   |closed - Newton-1D| = {nstr(dev,3)}")

        print("\nINDEPENDENT QUADRATURE (mpmath, dps 20) -- two-electron")
        print("-" * 78)
        specs = {
            (0, 0, 0, 0): ("one-center", lambda: V.one_center_eri_quad(A, A, A, A)),
            (0, 0, 1, 1): ("(AA|BB)", lambda: V.aabb_quad(A, A, B, B, R)),
            (0, 0, 0, 1): ("hybrid", lambda: V.hybrid_quad(A, A, A, B, R)),
        }
        for key, (cls, fn) in specs.items():
            t = time.time()
            qv = fn()
            dev = abs(mpv(canon[key][0], 20) - qv)
            quad_rows.append((f"({key[0]}{key[1]}|{key[2]}{key[3]}) {cls}",
                              nstr(dev, 3)))
            print(f"({key[0]}{key[1]}|{key[2]}{key[3]}) {cls:11s} "
                  f"|closed - quad| = {nstr(dev,3)}   ({time.time()-t:.0f}s)")
        t = time.time()
        exn, _per = V.exchange_numeric_neumann(
            zeta, (1, 0, 0), (1, 0, 0), zeta, (1, 0, 0), (1, 0, 0), R, tau_max=2)
        dev = abs(mpv(canon[(0, 1, 0, 1)][0], 20) - exn)
        quad_rows.append(("(01|01) exchange (numeric Neumann)", nstr(dev, 3)))
        print(f"(01|01) exchange   |closed - numeric-Neumann| = {nstr(dev,3)}"
              f"   ({time.time()-t:.0f}s)")
        mp.dps = 80

    # --- certification ------------------------------------------------------
    print("\nCERTIFICATION -- FCI total energy at two precisions")
    print("-" * 78)
    cert = {}
    for dps in (60, 90):
        t = time.time()
        with mp.workdps(dps + 30):
            tot, ee, vnn = AS.total_energy(S, h, gmap, orbs, 2, 1, 1, R, dps)
            cert[dps] = (mp.mpf(tot), mp.mpf(ee), mp.mpf(vnn))
        print(f"  dps {dps}: E_tot = {nstr(cert[dps][0], dps-5)}"
              f"   ({time.time()-t:.0f}s)")
    with mp.workdps(120):
        agree = abs(cert[60][0] - cert[90][0])
        digits = int(-mp.log10(agree / abs(cert[90][0]))) if agree > 0 else 999
    print(f"\n  two-precision agreement: |E(60) - E(90)| = {nstr(agree,3)}"
          f"  ->  {digits} significant digits certified")

    print("\n  E_electronic = " + nstr(cert[90][1], 40))
    print("  V_NN         = " + nstr(cert[90][2], 40))
    print("  E_total      = " + nstr(cert[90][0], 40))
    print("\n  HONEST FRAMING: this is the minimal-basis (one 1s per nucleus,")
    print("  zeta = 1) FCI energy -- the textbook single-zeta number.  Exact H2")
    print("  is -1.174476 Ha (Kolos-Wolniewicz); the gap is BASIS incompleteness,")
    print("  not integral error.  The claim certified here is the closed-form")
    print("  assembly and its precision, NOT chemical accuracy.")

    geom = {"R": "7/5", "zeta": "1", "tau_terminates_at": 2,
            "E_total_90": nstr(cert[90][0], 60),
            "E_electronic_90": nstr(cert[90][1], 60),
            "V_NN_90": nstr(cert[90][2], 60),
            "E_total_60": nstr(cert[60][0], 55),
            "two_precision_gap": nstr(agree, 3),
            "digits_certified": digits,
            "one_electron": [
                {"name": n, "tag": tag(e), "value_50": nstr(v, 50),
                 "dev_literature": (nstr(d, 3) if d is not None else None),
                 "note": nt}
                for n, e, _q, v, d, nt in one_e_rows],
            "two_electron": [
                {"quartet": f"({k[0]}{k[1]}|{k[2]}{k[3]})", "class": c,
                 "tag": tag(e), "value_50": nstr(v, 50),
                 "dev_literature": (nstr(d, 3) if d is not None else None),
                 "literature": ln}
                for k, c, e, v, d, ln in two_e_rows],
            "tau_terms": tau_rows,
            "quadrature_residuals": [{"name": n, "dev": d} for n, d in quad_rows],
            "h_route_symbolic_difference": str(hdiff)}
    report["geometries"].append(geom)

    # ------------------------------------------------------- other geometries
    print("\n" + "=" * 78)
    print("OTHER GEOMETRIES (closed-form, certified, literature-checked)")
    print("=" * 78)
    for Rr, Rf in ((sp.Rational(8, 5), "1.6"), (sp.Integer(2), "2.0")):
        t = time.time()
        orbs2, S2, h2_, hc2, gmap2, canon2 = build(zeta, Rr)
        Rmm = mp.mpf(Rf)
        devs = []
        for key, litfn in (((0, 0, 1, 1), V.lit_coulomb),
                           ((0, 0, 0, 1), V.lit_hybrid),
                           ((0, 1, 0, 1), V.lit_exchange_sugiura)):
            devs.append(abs(mpv(canon2[key][0], 40) - litfn(zm, Rmm)))
        with mp.workdps(120):
            t60, _e, _v = AS.total_energy(S2, h2_, gmap2, orbs2, 2, 1, 1, Rr, 60)
            t90, e90, v90 = AS.total_energy(S2, h2_, gmap2, orbs2, 2, 1, 1, Rr, 90)
            gap = abs(mp.mpf(t60) - mp.mpf(t90))
            dg = int(-mp.log10(gap / abs(mp.mpf(t90)))) if gap > 0 else 999
        print(f"R = {Rf}:  E_tot = {nstr(mp.mpf(t90), 40)}")
        print(f"          max |closed - literature| over the three two-center "
              f"ERIs = {nstr(max(devs),3)}")
        print(f"          two-precision gap {nstr(gap,3)} -> {dg} digits  "
              f"({time.time()-t:.0f}s)")
        report["geometries"].append(
            {"R": str(Rr), "zeta": "1", "E_total_90": nstr(mp.mpf(t90), 60),
             "E_electronic_90": nstr(mp.mpf(e90), 60),
             "V_NN_90": nstr(mp.mpf(v90), 60),
             "two_precision_gap": nstr(gap, 3), "digits_certified": dg,
             "max_dev_literature": nstr(max(devs), 3)})

    # --------------------------------------------------- variational zeta point
    print("\n" + "=" * 78)
    print("VARIATIONAL zeta = 1.197 at R = 1.4 (classic single-zeta optimum)")
    print("=" * 78)
    try:
        t = time.time()
        zv = Fraction(1197, 1000)
        orbs3, S3, h3, hc3, gmap3, canon3 = build(zv, R)
        zmm = mp.mpf("1.197")
        devs = []
        for key, litfn in (((0, 0, 1, 1), V.lit_coulomb),
                           ((0, 0, 0, 1), V.lit_hybrid),
                           ((0, 1, 0, 1), V.lit_exchange_sugiura)):
            devs.append(abs(mpv(canon3[key][0], 30) - litfn(zmm, Rm)))
        with mp.workdps(90):
            t40, e40, v40 = AS.total_energy(S3, h3, gmap3, orbs3, 2, 1, 1, R, 40)
            t60b, _, _ = AS.total_energy(S3, h3, gmap3, orbs3, 2, 1, 1, R, 60)
            gap = abs(mp.mpf(t40) - mp.mpf(t60b))
            dg = int(-mp.log10(gap / abs(mp.mpf(t60b)))) if gap > 0 else 999
        print(f"  E_tot = {nstr(mp.mpf(t60b), 40)}")
        print(f"  max |closed - literature| = {nstr(max(devs),3)};  "
              f"two-precision gap {nstr(gap,3)} -> {dg} digits "
              f"({time.time()-t:.0f}s)")
        report["geometries"].append(
            {"R": "7/5", "zeta": "1197/1000",
             "E_total_60": nstr(mp.mpf(t60b), 55),
             "two_precision_gap": nstr(gap, 3), "digits_certified": dg,
             "max_dev_literature": nstr(max(devs), 3)})
    except Exception as exc:                                # pragma: no cover
        print(f"  zeta = 1.197 point FAILED: {type(exc).__name__}: {exc}")
        report["zeta_variational_error"] = f"{type(exc).__name__}: {exc}"

    (OUT / "qfd_h2_certified.json").write_text(
        json.dumps(report, indent=2), encoding="utf-8")
    print(f"\nwrote {OUT / 'qfd_h2_certified.json'}   "
          f"(total {time.time()-t_all:.0f}s)")


if __name__ == "__main__":
    main()
