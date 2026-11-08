"""Poly-3: does the genuine-integral route actually give better polyatomics?

Paper 58 established the methodological break -- genuine two-centre orbital
products instead of atomic-sector surrogates -- and demonstrated it on NaH, a
DIATOMIC. It has never been run on a polyatomic. A quick scan (2026-08-14) put
BeH2 at R_eq +10.1% against the composed builder's 11.7%: a wash, but from two
different minimal bases, so the comparison was confounded.

This does it properly. A matched-basis A/B is not available even in principle:
the composed builder's block architecture cannot be run on an arbitrary orbital
set, because it carries DUPLICATE heavy-atom functions at different Z_eff per
block -- an arrangement that is only self-consistent because its cross-block
ERIs vanish identically. So instead of forcing a confounded A/B, this runs the
more informative experiment:

    THE BASIS LADDER. Grow the basis and watch where the geometry error goes.

That directly tests the load-bearing thesis (v4.73.0 A/B/C, plus the 2026-08-14
H2 two-axis split): polyatomic accuracy is limited by the ORBITAL BASIS, not by
integral fidelity and not by the missing cross-centre physics.

METHOD SPLIT, and why.
  Ladder  -> RHF. Correlation and basis are different axes; mixing them
             confounds exactly the question being asked. RHF geometry is
             well-defined, converges with basis, and is cheap enough to reach
             polarization functions (M = 19). A fifth rung adding heavy-atom 3d
             (M = 25) was dropped on cost -- ~64 s/point even at n_gauss = 6 --
             and is added back only if the L0-L3 trend is ambiguous.
  Anchor  -> FCI at the minimal basis, which is what the composed builder's
             11.7% / 19.4% figures are, so the comparison is like-for-like.
  (fci_ground here is DENSE over C(2M, N) determinants, so it is minimal-basis
  only: H2O at M = 9 would already be 43,758 dets and ~15 TB of dense H.)

PRE-REGISTERED PREDICTIONS (written before running):
  P1  BeH2 R_eq error decreases monotonically along the ladder and lands well
      inside the composed builder's 11.7%. A plateau near 10% would falsify the
      thesis -- something other than basis would be limiting the genuine route.
  P2  H2O behaves the same way; "bent" costs nothing extra. All three-centre
      integrals are present (McMurchie-Davidson is centre-agnostic), so
      bentness should not be a variable at all.
  P3  The water BOND ANGLE lands near experiment. This is the qualitative one:
      the composed builder structurally CANNOT produce an angle (MolecularSpec
      carries no nuclei; Paper 58 Obs. no_angle), so any sane angle is a
      capability the prior methodology did not have.

PRE-REGISTERED GATES:
  GO        both R_eq errors fall below 5% by the top rung AND the water angle
            lands in 100-115 deg
  PARTIAL   improvement real but stalls above 5%
  STOP      errors flat along the ladder -> basis is NOT the limiter, and the
            exact-vs-accurate framing needs revisiting for polyatomics

D_e is reported for BeH2 only, as a secondary observable, against fragments in
the SAME basis at each rung -- so basis-set superposition error is present and
UNCORRECTED. For H2O it is not reported at all: the O ground state is an
open-shell 3P, and closed-shell RHF is not a valid description of it. Stated,
not hidden.

Run from repo root:  python debug/poly3_genuine_polyatomic.py
"""

from __future__ import annotations

import json
import sys
import time
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[1]
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))

from geovac import noci_engine as GE  # noqa: E402

BOHR = 0.529177210903
HARTREE_EV = 27.211386245988

# Experimental references. BeH2 r_e(Be-H) = 1.3264 Ang (linear D_inf_h);
# H2O r_e(O-H) = 0.95784 Ang, theta_e = 104.508 deg.
EXP = {
    "BeH2": {"r": 1.3264 / BOHR, "angle": 180.0},
    "H2O": {"r": 0.95784 / BOHR, "angle": 104.508},
}
COMPOSED_REF = {"BeH2": 11.7, "H2O": 19.4}   # Paper 17, CLAUDE.md best-results table

P3 = ((1, 0, 0), (0, 1, 0), (0, 0, 1))
D6 = ((2, 0, 0), (0, 2, 0), (0, 0, 2), (1, 1, 0), (1, 0, 1), (0, 1, 1))
N_ELEC = {"BeH2": 6, "H2O": 10}
Z_HEAVY = {"BeH2": 4.0, "H2O": 8.0}
HEAVY = {"BeH2": "Be", "H2O": "O"}


def shapes(n_gauss: int = 6):
    """Slater radial shapes r^{n_r-1} e^{-r} fitted by n_gauss Gaussians.

    n_gauss = 6, NOT the 10 that the Phase 0-h finding recommends. That rule was
    written for validating a closed form against a reference at 1e-6, where the
    absolute fit error is the whole measurement. Here the observable is a
    geometry -- a DIFFERENCE over a smooth family -- so what matters is whether
    the fit error varies with geometry, not its absolute size. Cost forced the
    question: the build is O(n_gauss^4), measured 146 s vs 21 s per point at
    M = 19, and the first attempt at n_gauss = 10 ran 5.5 h without finishing.

    The choice is CONTROLLED, not assumed -- `ngauss_control()` re-runs the L0
    geometry at both 6 and 10 and reports the drift. Read that before trusting
    any number in the ladder.
    """
    out = {}
    for kind, (l, n_r) in (("1s", (0, 1)), ("2s", (0, 2)), ("2p", (1, 2)),
                           ("3s", (0, 3)), ("3p", (1, 3)), ("3d", (2, 3))):
        a, d, q = GE.fit_sto_shape(l, n_r, n_gauss=n_gauss)
        out[kind] = (a, d)
        out[kind + "_q"] = q
    return out


# Each rung is the FULL basis at that level (not a delta). Single-zeta
# Slater-rule exponents on the heavy atom; H standard.
LADDER = {
    "BeH2": [
        ("L0 minimal", {"Be": [("1s", 3.68), ("2s", 0.96), ("2p", 0.96)],
                        "H": [("1s", 1.00)]}),
        ("L1 +H 2s", {"Be": [("1s", 3.68), ("2s", 0.96), ("2p", 0.96)],
                      "H": [("1s", 1.00), ("2s", 0.80)]}),
        ("L2 +Be 3s3p", {"Be": [("1s", 3.68), ("2s", 0.96), ("2p", 0.96),
                                ("3s", 0.60), ("3p", 0.60)],
                         "H": [("1s", 1.00), ("2s", 0.80)]}),
        ("L3 +H 2p", {"Be": [("1s", 3.68), ("2s", 0.96), ("2p", 0.96),
                             ("3s", 0.60), ("3p", 0.60)],
                      "H": [("1s", 1.00), ("2s", 0.80), ("2p", 1.00)]}),
    ],
    "H2O": [
        ("L0 minimal", {"O": [("1s", 7.66), ("2s", 2.25), ("2p", 2.23)],
                        "H": [("1s", 1.00)]}),
        ("L1 +H 2s", {"O": [("1s", 7.66), ("2s", 2.25), ("2p", 2.23)],
                      "H": [("1s", 1.00), ("2s", 0.80)]}),
        ("L2 +O 3s3p", {"O": [("1s", 7.66), ("2s", 2.25), ("2p", 2.23),
                              ("3s", 1.20), ("3p", 1.20)],
                        "H": [("1s", 1.00), ("2s", 0.80)]}),
        ("L3 +H 2p", {"O": [("1s", 7.66), ("2s", 2.25), ("2p", 2.23),
                            ("3s", 1.20), ("3p", 1.20)],
                      "H": [("1s", 1.00), ("2s", 0.80), ("2p", 1.00)]}),
    ],
}


def _add_center(orbs, center, spec, sh):
    for kind, z in spec:
        if kind in ("2p", "3p"):
            for lmn in P3:
                orbs.append(GE.sto_shape_basis(center, kind, z, sh, lmn))
        elif kind == "3d":
            for lmn in D6:
                orbs.append(GE.sto_shape_basis(center, kind, z, sh, lmn))
        else:
            orbs.append(GE.sto_shape_basis(center, kind, z, sh, (0, 0, 0)))


def geometry(r, angle_deg):
    half = np.radians(angle_deg) / 2.0
    A = np.array([0.0, 0.0, 0.0])
    H1 = np.array([r * np.sin(half), 0.0, r * np.cos(half)])
    H2 = np.array([-r * np.sin(half), 0.0, r * np.cos(half)])
    return A, H1, H2


def build(system, rung, r, angle_deg, sh):
    A, H1, H2 = geometry(r, angle_deg)
    orbs = []
    _add_center(orbs, A, rung[HEAVY[system]], sh)
    _add_center(orbs, H1, rung["H"], sh)
    _add_center(orbs, H2, rung["H"], sh)
    nuc = [(A, Z_HEAVY[system]), (H1, 1.0), (H2, 1.0)]
    return orbs, nuc


def vnn(nuc):
    e = 0.0
    for i in range(len(nuc)):
        for j in range(i + 1, len(nuc)):
            e += nuc[i][1] * nuc[j][1] / np.linalg.norm(nuc[i][0] - nuc[j][0])
    return e


def rhf(S, h, g, n_elec, max_iter=200, tol=1e-10, damp=0.35):
    """Closed-shell RHF over a non-orthogonal basis. Returns electronic energy.

    Loewdin-orthogonalizes, core-Hamiltonian guess, linear density damping.
    g is chemist (pq|rs), so J_pq = sum_rs D_rs (pq|rs), K_pq = sum_rs D_rs (pr|qs).
    """
    assert n_elec % 2 == 0, "closed-shell RHF needs an even electron count"
    n_occ = n_elec // 2
    X = GE.lowdin_orbitals(S)
    D = np.zeros_like(h)
    e_old = 0.0
    for _ in range(max_iter):
        J = np.einsum("rs,pqrs->pq", D, g, optimize=True)
        K = np.einsum("rs,prqs->pq", D, g, optimize=True)
        F = h + J - 0.5 * K
        e_new = 0.5 * np.sum(D * (h + F))
        Fo = X.T @ F @ X
        _, C = np.linalg.eigh(Fo)
        C = X @ C
        Dn = 2.0 * (C[:, :n_occ] @ C[:, :n_occ].T)
        if abs(e_new - e_old) < tol and np.abs(Dn - D).max() < 1e-9:
            return e_new
        D = damp * D + (1.0 - damp) * Dn
        e_old = e_new
    return e_new                                  # non-converged; caller checks


def one_electron_exact(S, h):
    """Exact 1-electron ground state (H atom): lowest eigenvalue of h in the
    orthogonalized basis. RHF is not defined for an odd electron count."""
    X = GE.lowdin_orbitals(S)
    return float(np.linalg.eigvalsh(X.T @ h @ X)[0])


def total(system, rung, r, angle_deg, sh, method="rhf"):
    orbs, nuc = build(system, rung, r, angle_deg, sh)
    S, h, g = GE.integral_set_md(orbs, nuc)
    if method == "rhf":
        e = rhf(S, h, g, N_ELEC[system])
    else:
        X = GE.lowdin_orbitals(S)
        ht, gt = GE.transform_integrals(X, h, g)
        e = GE.fci_ground(ht, gt, N_ELEC[system])
    return e + vnn(nuc), len(orbs)


def beh2_fragments(rung, sh):
    """Be + 2H in the SAME basis (BSSE present, uncorrected)."""
    O = np.array([0.0, 0.0, 0.0])
    orbs = []
    _add_center(orbs, O, rung["Be"], sh)
    S, h, g = GE.integral_set_md(orbs, [(O, 4.0)])
    e_be = rhf(S, h, g, 4)
    orbs_h = []
    _add_center(orbs_h, O, rung["H"], sh)
    Sh, hh, _ = GE.integral_set_md(orbs_h, [(O, 1.0)])
    e_h = one_electron_exact(Sh, hh)
    return e_be + 2 * e_h


def _parab(xs, ys):
    c = np.polyfit(xs, ys, 2)
    if c[0] <= 0:
        return None, None
    x = -c[1] / (2 * c[0])
    return float(x), float(np.polyval(c, x))


def _scan(f, centre, width, n):
    xs = np.linspace(centre - width, centre + width, n)
    ys = [f(x) for x in xs]
    k = int(np.argmin(ys))
    if 0 < k < len(xs) - 1:
        xq, yq = _parab(xs[k - 1:k + 2], ys[k - 1:k + 2])
        if xq is not None:
            return xq, yq
    return float(xs[k]), float(ys[k])


def relax(system, rung, sh, method, r0, a0):
    """Alternating r / angle relaxation, coarse then fine."""
    r_eq, a_eq = r0, a0
    e = None
    for width_r, width_a, n in ((0.40, 20.0, 5), (0.10, 5.0, 5)):
        r_eq, e = _scan(lambda x: total(system, rung, x, a_eq, sh, method)[0],
                        r_eq, width_r, n)
        if system == "H2O":
            a_eq, e = _scan(lambda x: total(system, rung, r_eq, x, sh, method)[0],
                            a_eq, width_a, n)
    return r_eq, a_eq, e


def run(system, sh, results):
    ref = EXP[system]
    print(f"\n{'=' * 78}")
    print(f"{system}  --  genuine cross-centre integrals, 3-centre block INCLUDED")
    print(f"  experiment: r_e = {ref['r']:.4f} bohr"
          + (f", angle = {ref['angle']:.2f} deg" if system == "H2O" else "  (linear)"))
    print(f"  composed builder (Paper 17): R_eq error {COMPOSED_REF[system]}%"
          + ("  [and cannot represent an angle at all]" if system == "H2O" else ""))
    print(f"\n  {'rung':<14}{'M':>4}{'r_eq/bohr':>12}{'err':>9}"
          + (f"{'angle':>9}{'err':>8}" if system == "H2O" else "")
          + f"{'E_RHF':>15}{'s':>7}")
    print("  " + "-" * (74 if system == "H2O" else 58))

    for name, rung in LADDER[system]:
        t0 = time.time()
        _, M = total(system, rung, ref["r"] * 1.1, ref["angle"], sh)
        r_eq, a_eq, e = relax(system, rung, sh, "rhf", ref["r"] * 1.1, ref["angle"])
        err_r = 100 * (r_eq - ref["r"]) / ref["r"]
        row = {"rung": name, "M": M, "method": "rhf", "r_eq": r_eq,
               "err_r_pct": err_r, "angle": a_eq,
               "err_angle_deg": a_eq - ref["angle"], "E": e,
               "secs": round(time.time() - t0, 1)}
        if system == "BeH2":
            row["De_eV"] = (beh2_fragments(rung, sh) - e) * HARTREE_EV
        results.setdefault(system, []).append(row)
        ang = f"{a_eq:>9.2f}{a_eq - ref['angle']:>+8.2f}" if system == "H2O" else ""
        print(f"  {name:<14}{M:>4}{r_eq:>12.4f}{err_r:>+8.2f}%{ang}"
              f"{e:>15.6f}{row['secs']:>7.1f}")

    # FCI anchor at the minimal basis -- like-for-like with the composed figure
    name, rung = LADDER[system][0]
    t0 = time.time()
    r_eq, a_eq, e = relax(system, rung, sh, "fci", ref["r"] * 1.1, ref["angle"])
    err_r = 100 * (r_eq - ref["r"]) / ref["r"]
    results[system].append({"rung": "L0 minimal (FCI anchor)", "method": "fci",
                            "r_eq": r_eq, "err_r_pct": err_r, "angle": a_eq,
                            "err_angle_deg": a_eq - ref["angle"], "E": e,
                            "secs": round(time.time() - t0, 1)})
    ang = f"{a_eq:>9.2f}{a_eq - ref['angle']:>+8.2f}" if system == "H2O" else ""
    print(f"  {'L0 FCI anchor':<14}{'':>4}{r_eq:>12.4f}{err_r:>+8.2f}%{ang}"
          f"{e:>15.6f}{round(time.time() - t0, 1):>7.1f}")


def ngauss_control(results):
    """Is n_gauss = 6 good enough for a GEOMETRY? Measured, not assumed.

    Re-runs the L0 relaxation at n_gauss 6 and 10. If R_eq and the angle agree
    closely, the cheap fit is validated for this observable and the whole ladder
    is trustworthy at n_gauss = 6. If they do not, every number below is suspect
    and the ladder has to be re-costed instead of re-interpreted.
    """
    print(f"\n{'=' * 78}\nCONTROL -- is n_gauss = 6 enough for a geometry?")
    print(f"  {'system':<7}{'n_gauss':>9}{'r_eq/bohr':>12}{'angle':>10}{'s':>8}")
    print("  " + "-" * 48)
    ctrl = {}
    for system in ("BeH2", "H2O"):
        name, rung = LADDER[system][0]
        for ng in (6, 10):
            t0 = time.time()
            sh = shapes(ng)
            r_eq, a_eq, _ = relax(system, rung, sh, "rhf",
                                  EXP[system]["r"] * 1.1, EXP[system]["angle"])
            ctrl.setdefault(system, {})[ng] = {"r_eq": r_eq, "angle": a_eq}
            print(f"  {system:<7}{ng:>9}{r_eq:>12.4f}{a_eq:>10.2f}"
                  f"{time.time() - t0:>8.1f}")
        d_r = abs(ctrl[system][6]["r_eq"] - ctrl[system][10]["r_eq"])
        d_a = abs(ctrl[system][6]["angle"] - ctrl[system][10]["angle"])
        pct = 100 * d_r / EXP[system]["r"]
        ok = "PASS" if pct < 0.5 else "FAIL"
        print(f"  {system} drift: {d_r:.5f} bohr ({pct:.3f}% of r_e), "
              f"{d_a:.3f} deg  -> {ok}")
        ctrl[system]["drift_pct_of_re"] = pct
        ctrl[system]["drift_deg"] = d_a
        ctrl[system]["verdict"] = ok
    results["ngauss_control"] = ctrl
    return all(ctrl[s]["verdict"] == "PASS" for s in ("BeH2", "H2O"))


def main():
    print("Poly-3 -- the genuine-integral route on polyatomics, basis ladder")
    results = {}
    ok = ngauss_control(results)
    if not ok:
        print("\n  !! n_gauss control FAILED -- ladder numbers below are not "
              "trustworthy at n_gauss=6; re-cost before interpreting.")

    sh = shapes(n_gauss=6)
    print("\n  <fit|STO> at n_gauss=6: "
          + ", ".join(f"{k}:{sh[k + '_q']:.7f}"
                      for k in ("1s", "2s", "2p", "3s", "3p", "3d")))
    for system in ("BeH2", "H2O"):
        run(system, sh, results)

    out = REPO / "debug" / "data" / "poly3_genuine_polyatomic.json"
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(results, indent=2))

    print(f"\n{'=' * 78}\nGATE")
    for system in ("BeH2", "H2O"):
        rhf_rows = [r for r in results[system] if r["method"] == "rhf"]
        first, last = rhf_rows[0]["err_r_pct"], rhf_rows[-1]["err_r_pct"]
        verdict = ("GO" if abs(last) < 5.0 else
                   "PARTIAL" if abs(last) < abs(first) - 0.5 else "STOP")
        print(f"  {system:<5} r_eq {first:+.2f}% -> {last:+.2f}%   "
              f"(composed {COMPOSED_REF[system]}%)   {verdict}")
    wa = [r for r in results["H2O"] if r["method"] == "rhf"][-1]["angle"]
    print(f"  H2O angle at top rung: {wa:.2f} deg   "
          f"(gate 100-115) {'PASS' if 100 <= wa <= 115 else 'FAIL'}")
    print(f"\nwrote {out}")


if __name__ == "__main__":
    main()
