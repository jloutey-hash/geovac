"""QFD table extension -- certified minimal diatomics beyond H2 and LiH.

Extends the quadrature-free diatomic (QFD) build (geovac/qfd_core.py +
geovac/qfd_assemble.py, Paper 58, CHANGELOG v4.105.0) from the two systems it
was demonstrated on (H2, LiH) to a small TABLE of s-only minimal diatomics that
an s-only basis can honestly hold:

    H2+     1 electron   1s_A / 1s_B, zeta = 1        R = 2.0 and 1.4 bohr
    HeH+    2 electrons  He 1s (Z=2) / H 1s (Z=1)     R = 1.46 bohr
    He2^2+  2 electrons  1s / 1s, Z = 2 both          R = 1.3  bohr
    BeH+    4 electrons  Be 1s,2s (Z=4) / H 1s        R = 2.5  bohr

Basis convention: HYDROGENIC, decay rate a = Z_orbital / n (the
`geovac.two_center_eri` convention -- NOT Coulomb-Sturmian, NOT Gaussian; see
memory/polyatomic_state_of_play.md).  All functions are l = 0.

HONEST SCOPE.  These are minimal s-only bases with unoptimized hydrogenic
exponents.  The energies are far from exact and NO accuracy claim is made; what
is certified is the closed-form assembly and its digit count.  Where a textbook
minimal-basis number exists (H2+ LCAO) it is quoted as a sanity anchor, not as
a validation of chemistry.

THE ONE STRUCTURAL AXIS.  The Neumann tau series of the EXCHANGE class
TERMINATES exactly when both centres of a density carry the same orbital
exponent (q = (alpha-beta)R/2 = 0, the Phase 0-e criterion).

    H2+, He2^2+   homonuclear at equal Z_orbital -> q = 0 -> FINITE sum,
                  no truncation anywhere, energy exact to working precision.
    HeH+, BeH+    heteronuclear -> infinite sum, truncated at a per-quartet
                  tau_max whose tail is MEASURED and BOUNDED here, with the
                  amplification factor recomputed for each system (never copied
                  from LiH).

Run:
    python -u debug/qfd_table_ext.py                 # everything
    python -u debug/qfd_table_ext.py --only he2_2plus
    python -u debug/qfd_table_ext.py --skip-quad
"""
from __future__ import annotations

import argparse
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

from geovac import qfd_assemble as AS       # noqa: E402
from geovac import qfd_core as Q            # noqa: E402
import qfd_quad as V                        # noqa: E402

OUT = REPO / "debug" / "data"
OUT.mkdir(parents=True, exist_ok=True)
JSON_PATH = OUT / "qfd_table_ext_certified.json"

DPS_LO, DPS_HI = 40, 60          # the two-precision certification pair
EXCH_DPS = 50                    # working precision of the exchange accumulator
EXCH_QUAD_TAU = 3                # tau depth of the numeric-Neumann cross-check


# ---------------------------------------------------------------- system table

def sys_h2plus(rlabel, R):
    return dict(
        key="h2plus_R" + rlabel, name="H2+",
        R=R, rlabel=rlabel, ZA=1, ZB=1, n_elec=1,
        orbs=[("A", Fraction(1), 1), ("B", Fraction(1), 1)],
        labels=["H_A 1s", "H_B 1s"],
        homonuclear=True, tau={}, tau_default=2, tau_verify=8,
        basis="1s on each proton, hydrogenic a = 1",
        note="one electron: the 2x2 generalized secular problem (h, S); the "
             "two-electron tensor is assembled but cannot contribute.")


def sys_heh():
    return dict(
        key="hehplus", name="HeH+",
        R=sp.Rational(73, 50), rlabel="1.46", ZA=2, ZB=1, n_elec=2,
        orbs=[("A", Fraction(2), 1), ("B", Fraction(1), 1)],
        labels=["He 1s", "H 1s"],
        homonuclear=False,
        tau={(1, 1, 1, 1): 16}, tau_default=16, tau_verify=None,
        basis="He 1s a=2, H 1s a=1 (hydrogenic)",
        note="heteronuclear: the exchange tau series does not terminate.")


def sys_he2():
    return dict(
        key="he2_2plus", name="He2^2+",
        R=sp.Rational(13, 10), rlabel="1.3", ZA=2, ZB=2, n_elec=2,
        orbs=[("A", Fraction(2), 1), ("B", Fraction(2), 1)],
        labels=["He_A 1s", "He_B 1s"],
        homonuclear=True, tau={}, tau_default=2, tau_verify=8,
        basis="1s on each He, hydrogenic a = 2",
        note="homonuclear at equal exponent: the exchange tau series TERMINATES.")


def sys_behplus():
    return dict(
        key="behplus", name="BeH+",
        R=sp.Rational(5, 2), rlabel="2.5", ZA=4, ZB=1, n_elec=4,
        orbs=[("A", Fraction(4), 1), ("A", Fraction(4), 2),
              ("B", Fraction(1), 1)],
        labels=["Be 1s", "Be 2s", "H 1s"],
        homonuclear=False,
        # q = (a_A - a_B) R / 2 sets the decay: (Be1s,H) q = 3.75,
        # (Be2s,H) q = 1.25.  Larger q -> slower Neumann convergence.
        tau={(1, 1, 1, 1): 22, (1, 1, 2, 1): 18, (2, 1, 1, 1): 18,
             (2, 1, 2, 1): 16},
        tau_default=22, tau_verify=None,
        basis="Be 1s a=4, Be 2s a=2, H 1s a=1 (hydrogenic)",
        note="heteronuclear, 4 electrons, FCI dimension C(6,4) = 15.")


ALL_SYSTEMS = [sys_h2plus("2.0", sp.Integer(2)),
               sys_h2plus("1.4", sp.Rational(7, 5)),
               sys_heh(), sys_he2(), sys_behplus()]


# --------------------------------------------------------------------- helpers

def nstr(x, n=25):
    return mp.nstr(mp.mpf(x), n, strip_zeros=False)


def mpv(expr, dps=50):
    if isinstance(expr, mp.mpf):
        return +expr
    with mp.workdps(dps + 15):
        return mp.mpf(str(sp.N(expr, dps + 10)))


def tag(expr):
    """Transcendence seeds actually present in the built closed form."""
    if isinstance(expr, mp.mpf):
        return "{exp, E_1, ln, gamma}"          # exchange class, Paper 58/59
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
    present = [nm for nm in ("exp", "E_1", "ln", "gamma") if nm in names]
    if present in ([], ["exp"]):
        return "elementary {exp}"
    return "{" + ", ".join(present) + "}"


def agree_digits(a, b):
    a, b = mp.mpf(a), mp.mpf(b)
    if a == b:
        return 10 ** 6
    scale = abs(b) if b != 0 else mp.mpf(1)
    return int(mp.floor(-mp.log10(abs(a - b) / scale)))


def tau_tail_bound(per):
    """Geometric bound on sum_{tau > tau_max} a_tau from the measured terms.

    The observed |a_tau| decrease monotonically in the tail (factorial
    convergence -- the per-tau table is printed so this is checked, not merely
    asserted), so with r the last observed ratio,

        |tail| <= |a_taumax| * r / (1 - r).

    Returns (bound, |a_taumax|, r).  If the last ratio is >= 1 the series has
    not entered its decaying regime and the last term itself is returned as the
    (weak) bound, flagged by r >= 1.
    """
    nz = [(i, abs(v)) for i, v in enumerate(per) if v != 0]
    if not nz:
        return mp.mpf(0), mp.mpf(0), mp.mpf(0)
    if len(nz) < 2:
        return nz[-1][1], nz[-1][1], mp.mpf(0)
    (i_last, a_last), (i_prev, a_prev) = nz[-1], nz[-2]
    r = (a_last / a_prev) ** (mp.mpf(1) / (i_last - i_prev))
    if r >= 1:
        return a_last, a_last, r
    return a_last * r / (1 - r), a_last, r


# ------------------------------------------------------------------ one system

def run_system(spec, skip_quad=False):
    key, name, R = spec["key"], spec["name"], spec["R"]
    orbs, ZA, ZB = spec["orbs"], spec["ZA"], spec["ZB"]
    n, N = len(orbs), spec["n_elec"]
    LAB = spec["labels"]
    rep = {"key": key, "system": name, "R": str(sp.nsimplify(R)),
           "R_bohr": nstr(mpv(R, 30), 12), "basis": spec["basis"],
           "n_electrons": N, "n_orbitals": n,
           "homonuclear": spec["homonuclear"], "note": spec["note"]}
    t_sys = time.time()

    print("=" * 78)
    print(name + "   |   " + spec["basis"] + "   |   R = " + rep["R_bohr"]
          + " bohr   |   " + str(N) + " electron(s)")
    print("=" * 78)

    # -------------------------------------------------------- one-electron
    t = time.time()
    S, h = AS.build_S_h(orbs, ZA, ZB, R)
    hchk = AS.h_core_check(orbs, ZA, ZB, R)
    hdiff_zero = all(sp.simplify(h[i][j] - hchk[i][j]) == 0
                     for i in range(n) for j in range(n))
    print("one-electron closed forms built in %.1fs; two independent h routes "
          "agree symbolically: %s" % (time.time() - t, hdiff_zero))
    rep["h_two_routes_symbolically_identical"] = bool(hdiff_zero)

    rows = []
    print("\nOVERLAP / CORE HAMILTONIAN (closed form, 30 digits)")
    print("-" * 78)
    for i in range(n):
        for j in range(i, n):
            sv, hv = mpv(S[i][j]), mpv(h[i][j])
            print("  S[%-7s,%-7s] = %34s  [%s]"
                  % (LAB[i], LAB[j], nstr(sv, 30), tag(S[i][j])))
            print("  h[%-7s,%-7s] = %34s  [%s]"
                  % (LAB[i], LAB[j], nstr(hv, 30), tag(h[i][j])))
            rows.append({"i": LAB[i], "j": LAB[j], "S_50": nstr(sv, 50),
                         "h_50": nstr(hv, 50), "S_tag": tag(S[i][j]),
                         "h_tag": tag(h[i][j])})
    rep["one_electron"] = rows

    # ------------------------------------------- one-electron vs quadrature
    if not skip_quad:
        print("\nONE-ELECTRON vs INDEPENDENT QUADRATURE (raw prolate "
              "spheroidal, dps 20)")
        print("-" * 78)
        qrows = []
        mp.dps = 20
        cross = [(i, j) for i in range(n) for j in range(n)
                 if orbs[i][0] != orbs[j][0]]
        same = [(i, j) for i in range(n) for j in range(i, n)
                if orbs[i][0] == orbs[j][0]]
        specs = []
        for (i, j) in cross[:2]:
            specs += [("S", i, j, "S"), ("T", i, j, "T"),
                      ("V^A", i, j, "inv_ra")]
        for (i, j) in same[:2]:
            other = "B" if orbs[i][0] == "A" else "A"
            specs.append(("V^" + other, i, j,
                          "inv_ra" if other == "A" else "inv_rb"))
        for nm, i, j, kern in specs:
            if kern == "S":
                cf = Q.overlap(orbs[i], orbs[j], R)
            elif kern == "T":
                cf = Q.kinetic(orbs[i], orbs[j], R)
            else:
                cf = Q._inv_r(orbs[i], orbs[j],
                              "A" if kern == "inv_ra" else "B", R)
            qv = V.one_electron_quad(orbs[i], orbs[j], kern, R)
            dev = abs(mpv(cf, 20) - qv)
            print("  %-4s[%-7s,%-7s]  |closed - quad| = %s"
                  % (nm, LAB[i], LAB[j], nstr(dev, 3)))
            qrows.append({"name": "%s[%s,%s]" % (nm, LAB[i], LAB[j]),
                          "abs_dev": nstr(dev, 3)})
        mp.dps = 80
        rep["one_electron_quadrature"] = qrows

        # literature closed forms, available for the zeta = 1 H2+ pair
        if name == "H2+":
            mp.dps = 40
            Rm = mp.mpf(str(sp.N(R, 45)))
            lits = [("S_AB", V.lit_overlap(mp.mpf(1), Rm),
                     Q.overlap(orbs[0], orbs[1], R)),
                    ("T_AB", V.lit_kinetic_cross(mp.mpf(1), Rm),
                     Q.kinetic(orbs[0], orbs[1], R)),
                    ("V^B_AA", V.lit_nuc_same(mp.mpf(1), Rm),
                     Q._inv_r(orbs[0], orbs[0], "B", R)),
                    ("V^A_AB", V.lit_nuc_cross(mp.mpf(1), Rm),
                     Q._inv_r(orbs[0], orbs[1], "A", R))]
            print("\n  vs LITERATURE closed forms (classical 1s two-centre "
                  "formulas, dps 40)")
            lrows = []
            for nm, lv, cf in lits:
                dev = abs(mpv(cf, 40) - lv)
                print("    %-8s |closed - literature| = %s" % (nm, nstr(dev, 3)))
                lrows.append({"name": nm, "abs_dev": nstr(dev, 3)})
            mp.dps = 80
            rep["literature_checks"] = lrows

    # -------------------------------------------------------- two-electron
    use_hp = not spec["homonuclear"]
    tau_map = spec["tau"]
    tdef = spec["tau_default"]

    def tau_for(quart):
        a, b, c, d = quart
        return tau_map.get((a[2], b[2], c[2], d[2]), tdef)

    print("\nTWO-ELECTRON closed forms  (%s; tau_max %s)"
          % ("numeric tau accumulation at dps %d" % EXCH_DPS if use_hp
             else "fully symbolic", tau_map if tau_map else tdef))
    print("-" * 78)
    AS.EXCHANGE_TAILS.clear()
    t = time.time()
    gmap, canon = AS.build_g(orbs, R, tau_max=tau_for,
                             exchange_hp_dps=EXCH_DPS if use_hp else None)
    t_g = time.time() - t
    print("  (built in %.0fs)\n" % t_g)
    g_rows = []
    for k, (expr, cls) in sorted(canon.items()):
        val = mpv(expr)
        nm = "(%s %s|%s %s)" % (LAB[k[0]], LAB[k[1]], LAB[k[2]], LAB[k[3]])
        print("  %-33s %-11s %34s  [%s]" % (nm, cls, nstr(val, 30), tag(expr)))
        g_rows.append({"quartet": nm, "index": list(k), "class": cls,
                       "tag": tag(expr), "value_50": nstr(val, 50)})
    rep["two_electron"] = g_rows
    rep["g_build_seconds"] = round(t_g)

    # ------------------------------------ exchange tau: terminate, or bound it
    exch_keys = [k for k, (_e, c) in canon.items() if c == "exchange"]
    tail_rows, total_bound = [], mp.mpf(0)
    if spec["homonuclear"]:
        print("\nEXCHANGE tau series -- TERMINATION CHECK (q = 0, homonuclear "
              "at equal exponent)")
        print("-" * 78)
        term_rows = []
        for k in sorted(exch_keys):
            a, b, c, d = (orbs[i] for i in k)
            if a[0] != "A":
                a, b = b, a
            if c[0] != "A":
                c, d = d, c
            terms = Q.exchange_closed_form(
                a[1], (a[2], 0, 0), (b[2], 0, 0), b[1], (c[2], 0, 0),
                (d[2], 0, 0), R, tau_max=spec["tau_verify"], return_terms=True)
            simp = [sp.simplify(x) for x in terms]
            nonzero = [i for i, x in enumerate(simp) if x != 0]
            last_nz = max(nonzero) if nonzero else -1
            ok = all(simp[i] == 0 for i in range(last_nz + 1, len(simp)))
            nm = "(%s %s|%s %s)" % (LAB[k[0]], LAB[k[1]], LAB[k[2]], LAB[k[3]])
            print("  %-33s last nonzero tau = %d; every tau in [%d, %d] is a "
                  "SYMBOLIC ZERO: %s"
                  % (nm, last_nz, last_nz + 1, spec["tau_verify"], ok))
            print("      per-tau: "
                  + ", ".join(nstr(mpv(x, 25), 4) for x in terms))
            term_rows.append({"quartet": nm, "last_nonzero_tau": last_nz,
                              "verified_zero_through_tau": spec["tau_verify"],
                              "terminates": bool(ok),
                              "per_tau": [nstr(mpv(x, 25), 6) for x in terms]})
        rep["exchange_termination"] = term_rows
        rep["tau_truncation_error"] = "0 (finite sum -- no truncation)"
    else:
        print("\nEXCHANGE tau-truncation audit (heteronuclear: q != 0, the "
              "series is INFINITE)")
        print("-" * 78)
        for kk, per in sorted(AS.EXCHANGE_TAILS.items()):
            tot = mp.fsum(per)
            bound, a_last, r = tau_tail_bound(per)
            total_bound += bound
            half = len(per) // 2
            monotone = all(abs(per[i + 1]) <= abs(per[i])
                           for i in range(half, len(per) - 1))
            print("  n-quartet (%d,%d|%d,%d)  tau_max = %d   value = %s"
                  % (kk[0], kk[1], kk[2], kk[3], len(per) - 1, nstr(tot, 25)))
            print("      |a_taumax| = %s   last ratio r = %s   tail <= %s   "
                  "(tail terms monotone: %s)"
                  % (nstr(a_last, 3), nstr(r, 3), nstr(bound, 3), monotone))
            tail_rows.append({"n_quartet": list(kk), "tau_max": len(per) - 1,
                              "value_40": nstr(tot, 40),
                              "last_term_abs": nstr(a_last, 3),
                              "last_ratio": nstr(r, 3),
                              "tail_bound": nstr(bound, 3),
                              "tail_terms_monotone": bool(monotone),
                              "per_tau_rel": [nstr(abs(v / tot), 3) if tot != 0
                                              else "0" for v in per]})
        rep["exchange_tau_tails"] = tail_rows

    # ------------------------------------------- two-electron vs quadrature
    if not skip_quad:
        print("\nTWO-ELECTRON vs INDEPENDENT QUADRATURE (dps 20)")
        print("-" * 78)
        mp.dps = 20
        q2 = []
        for k, (expr, cls) in sorted(canon.items()):
            oa, ob, oc, od = (orbs[i] for i in k)
            nm = "(%s %s|%s %s)" % (LAB[k[0]], LAB[k[1]], LAB[k[2]], LAB[k[3]])
            t = time.time()
            try:
                if cls == "exchange":
                    a, b, c, d = oa, ob, oc, od
                    if a[0] != "A":
                        a, b = b, a
                    if c[0] != "A":
                        c, d = d, c
                    qv, _ = V.exchange_numeric_neumann(
                        a[1], (a[2], 0, 0), (b[2], 0, 0), b[1],
                        (c[2], 0, 0), (d[2], 0, 0), R, tau_max=EXCH_QUAD_TAU)
                    # tau-MATCHED comparison: rebuild the closed form at the
                    # same truncation, else the deviation is dominated by the
                    # tau terms the quadrature reference simply omits.
                    if use_hp:
                        cfv, _ = Q.exchange_hp(
                            a[1], (a[2], 0, 0), (b[2], 0, 0), b[1],
                            (c[2], 0, 0), (d[2], 0, 0), R,
                            tau_max=EXCH_QUAD_TAU, dps=25)
                        cfv = mp.mpf(str(cfv))
                    else:
                        cfv = mpv(Q.exchange_closed_form(
                            a[1], (a[2], 0, 0), (b[2], 0, 0), b[1],
                            (c[2], 0, 0), (d[2], 0, 0), R,
                            tau_max=EXCH_QUAD_TAU), 20)
                    dev = abs(cfv - qv)
                    print("  %-33s %-11s |closed - quad| = %s   (%.0fs)  "
                          "(tau-matched at tau <= %d)"
                          % (nm, cls, nstr(dev, 3), time.time() - t,
                             EXCH_QUAD_TAU))
                    q2.append({"index": list(k), "class": cls,
                               "abs_dev": nstr(dev, 3),
                               "note": "tau-matched at tau <= %d"
                                       % EXCH_QUAD_TAU})
                    continue
                if cls == "one-center":
                    qv = V.one_center_eri_quad(oa, ob, oc, od)
                elif cls == "(AA|BB)":
                    qv = V.aabb_quad(oa, ob, oc, od, R)
                else:
                    qv = (V.hybrid_quad(oa, ob, oc, od, R) if oa[0] == ob[0]
                          else V.hybrid_quad(oc, od, oa, ob, R))
            except Exception as exc:                      # pragma: no cover
                print("  %-33s %-11s quadrature FAILED %s: %s"
                      % (nm, cls, type(exc).__name__, exc))
                q2.append({"index": list(k), "class": cls,
                           "abs_dev": "FAILED %s: %s"
                                      % (type(exc).__name__, exc)})
                continue
            dev = abs(mpv(expr, 20) - qv)
            print("  %-33s %-11s |closed - quad| = %s   (%.0fs)"
                  % (nm, cls, nstr(dev, 3), time.time() - t))
            q2.append({"index": list(k), "class": cls, "abs_dev": nstr(dev, 3)})
        mp.dps = 80
        rep["two_electron_quadrature"] = q2

    # -------------------------------------------------------- certification
    print("\nCERTIFICATION -- total energy at dps %d and %d" % (DPS_LO, DPS_HI))
    print("-" * 78)
    cert = {}
    for dps in (DPS_LO, DPS_HI):
        t = time.time()
        with mp.workdps(dps + 40):
            tot, ee, vnn = AS.total_energy(S, h, gmap, orbs, N, ZA, ZB, R, dps)
            cert[dps] = (mp.mpf(tot), mp.mpf(ee), mp.mpf(vnn))
        print("  dps %d: E_tot = %s   (%.0fs)"
              % (dps, nstr(cert[dps][0], min(dps, 45)), time.time() - t))

    with mp.workdps(150):
        gap = abs(cert[DPS_LO][0] - cert[DPS_HI][0])
        lin = min(agree_digits(cert[DPS_LO][0], cert[DPS_HI][0]), DPS_LO)
    rep["two_precision_gap"] = nstr(gap, 3) if gap > 0 else "0 (bit-identical)"
    rep["digits_from_linear_algebra"] = lin

    # amplification factor, computed FOR THIS SYSTEM (never copied)
    with mp.workdps(80):
        Sm, _hm, _gm = AS.evaluate(S, h, gmap, n, 60)
        lam = min(mp.eigsy(Sm, eigvals_only=True))
        Xnorm4 = 1 / lam ** 2
        rdm_weight = N * (N - 1) // 2
        raw_amp = rdm_weight * Xnorm4
    rep["lambda_min_S"] = nstr(lam, 12)
    rep["amplification_raw"] = nstr(raw_amp, 6)

    if spec["homonuclear"] or N < 2:
        tau_dig = 10 ** 6
        e_bound = mp.mpf(0)
        print("\n  exchange tau series terminates (or there is no two-electron "
              "term at all): truncation error is exactly 0")
    else:
        amp = 1
        while amp < 4 * raw_amp:              # >= 4x margin, next power of two
            amp *= 2
        e_bound = amp * total_bound
        tau_dig = (int(mp.floor(-mp.log10(e_bound / abs(cert[DPS_HI][0]))))
                   if e_bound > 0 else 10 ** 6)
        print("\n  sum of per-quartet tail bounds = %s" % nstr(total_bound, 3))
        print("  amplification: 2-RDM weight N(N-1)/2 = %d times ||S^-1/2||^4 "
              "= lambda_min^-2 = %s  ->  %s, rounded up to %d"
              % (rdm_weight, nstr(Xnorm4, 6), nstr(raw_amp, 6), amp))
        print("  energy-level tau-tail bound   = %s Ha   ->  %d digits"
              % (nstr(e_bound, 3), tau_dig))
        rep["amplification_used"] = amp
        rep["energy_tau_tail_bound_Ha"] = nstr(e_bound, 3)
    rep["digits_from_tau_truncation"] = (None if tau_dig == 10 ** 6
                                         else tau_dig)

    # exchange accumulator arithmetic self-check (heteronuclear only)
    exch_dig = 10 ** 6
    if use_hp and exch_keys:
        k = sorted(exch_keys)[0]
        a, b, c, d = (orbs[i] for i in k)
        if a[0] != "A":
            a, b = b, a
        if c[0] != "A":
            c, d = d, c
        tm = 8
        lo, _ = Q.exchange_hp(a[1], (a[2], 0, 0), (b[2], 0, 0), b[1],
                              (c[2], 0, 0), (d[2], 0, 0), R, tau_max=tm, dps=30)
        hi, _ = Q.exchange_hp(a[1], (a[2], 0, 0), (b[2], 0, 0), b[1],
                              (c[2], 0, 0), (d[2], 0, 0), R, tau_max=tm,
                              dps=EXCH_DPS)
        with mp.workdps(80):
            dexc = abs(lo - hi)
            exch_dig = agree_digits(lo, hi)
        print("  exchange accumulator, dps 30 vs %d at fixed tau = %d: "
              "|difference| = %s  (%d digits)"
              % (EXCH_DPS, tm, nstr(dexc, 3), exch_dig))
        rep["exchange_two_precision"] = nstr(dexc, 3)
        rep["digits_from_exchange_arithmetic"] = exch_dig

    digits = min(lin, tau_dig, exch_dig)
    rep["digits_certified"] = digits
    print("\n  CERTIFIED DIGITS = min(linear algebra %d, tau tail %s, exchange "
          "arithmetic %s) = %d"
          % (lin, "inf" if tau_dig == 10 ** 6 else str(tau_dig),
             "inf" if exch_dig == 10 ** 6 else str(exch_dig), digits))

    rep["E_total"] = nstr(cert[DPS_HI][0], digits)
    rep["E_electronic"] = nstr(cert[DPS_HI][1], digits)
    rep["V_NN"] = nstr(cert[DPS_HI][2], digits)
    rep["E_total_full"] = nstr(cert[DPS_HI][0], DPS_HI)
    print("\n  E_electronic = " + rep["E_electronic"])
    print("  V_NN         = " + rep["V_NN"])
    print("  E_total      = " + rep["E_total"])

    # ------------------------------------------------- system-specific anchors
    if name == "H2+":
        with mp.workdps(60):
            Sm2, hm2, gm2 = AS.evaluate(S, h, gmap, n, 50)
            X = AS.lowdin(Sm2, n)
            ht, _gt = AS.transform(X, hm2, gm2, n)
            gen = min(mp.eigsy(ht, eigvals_only=True))
            lcao = (hm2[0, 0] + hm2[0, 1]) / (1 + Sm2[0, 1])
            d1 = min(agree_digits(gen, cert[DPS_HI][1]), 40)
            d2 = min(agree_digits(lcao, cert[DPS_HI][1]), 40)
        print("\n  ANCHOR: the 1-electron FCI equals the lowest generalized "
              "eigenvalue of (h, S) to %d digits" % d1)
        print("  ANCHOR: and equals the closed LCAO sigma_g form "
              "(h_AA + h_AB)/(1 + S_AB) to %d digits" % d2)
        print("  Textbook context: the unscreened (zeta = 1) LCAO H2+ curve has "
              "its minimum at\n    R = 2.49 bohr with E = -0.5648 Ha; the values "
              "here are that same minimal-basis\n    curve evaluated at the "
              "stated R, not at its own minimum.")
        rep["anchor_generalized_eigenvalue_digits"] = d1
        rep["anchor_lcao_formula_digits"] = d2
        rep["anchor_note"] = (
            "the one-electron FCI reproduces both the lowest generalized "
            "eigenvalue of (h, S) and the closed LCAO sigma_g expression "
            "(h_AA + h_AB)/(1 + S_AB).  Textbook minimal-basis context: the "
            "zeta = 1 LCAO H2+ curve minimizes at R = 2.49 bohr, E = -0.5648 Ha.")

    rep["seconds"] = round(time.time() - t_sys)
    print("\n  (%ds for this system)\n" % rep["seconds"])
    return rep


# ------------------------------------------------------------------------ main

def main(argv=None):
    ap = argparse.ArgumentParser(description="QFD diatomic table extension")
    ap.add_argument("--only", default="", help="comma-separated system keys")
    ap.add_argument("--skip-quad", action="store_true")
    ap.add_argument("--out", default=str(JSON_PATH),
                    help="output JSON (use a distinct file when running "
                         "systems in parallel processes, then --merge)")
    ap.add_argument("--merge", default="",
                    help="comma-separated JSON shards to merge into --out "
                         "before/instead of running")
    args = ap.parse_args(argv)

    out_path = Path(args.out)
    mp.dps = 80
    keys = [k.strip() for k in args.only.split(",") if k.strip()]
    todo = [s for s in ALL_SYSTEMS if not keys or s["key"] in keys]
    assert todo, "no system matched %s" % keys

    store = {}
    if out_path.exists():
        store = json.loads(out_path.read_text(encoding="utf-8"))
    store.setdefault("systems", {})
    for shard in [s.strip() for s in args.merge.split(",") if s.strip()]:
        sh = json.loads(Path(shard).read_text(encoding="utf-8"))
        store["systems"].update(sh.get("systems", {}))
        print("merged %d system(s) from %s" % (len(sh.get("systems", {})),
                                               shard))
    store["meta"] = {"dps_pair": [DPS_LO, DPS_HI], "exchange_dps": EXCH_DPS,
                     "basis_convention": "hydrogenic, rate a = Z_orbital / n; "
                                         "s-type (l = 0) only",
                     "generator": "debug/qfd_table_ext.py"}

    if args.merge and not args.only:
        todo = []
    for spec in todo:
        rep = run_system(spec, skip_quad=args.skip_quad)
        store["systems"][spec["key"]] = rep
        out_path.write_text(json.dumps(store, indent=2), encoding="utf-8")
        print("wrote %s" % out_path)
    if not todo:
        out_path.write_text(json.dumps(store, indent=2), encoding="utf-8")
        print("wrote %s" % out_path)

    print("\n" + "=" * 78)
    print("SUMMARY")
    print("=" * 78)
    print("%-12s %7s %2s %7s  E_total" % ("system", "R", "N", "digits"))
    for _k, r in store["systems"].items():
        print("%-12s %7s %2d %7d  %s"
              % (r["system"], r["R_bohr"], r["n_electrons"],
                 r["digits_certified"], r["E_total"]))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
