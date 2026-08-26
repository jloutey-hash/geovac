"""NaH basis-extension probe: does the in-basis FCI ceiling rise?

QUESTION (PI direction 2026-08-22). Paper 58's NaH result attributes its D_e
shortfall to minimal-basis incompleteness, naming three missing ingredients:
polarization functions, diffuse functions, and BSSE correction. None has been
tried. This driver adds the first two and measures whether the ceiling moves.

THE MEASUREMENT. NOCI-3 is a compact *approximation* to in-basis FCI, and it
already recovers 91.1% of it. So NOCI is not the binding constraint -- the
BASIS is. The number to watch is therefore the in-basis FCI D_e:

    stored baseline (M=7):  FCI D_e = 0.043189 Ha = 1.175 eV = 59.9% of exp
                            NOCI-3   = 0.039357 Ha = 1.071 eV = 54.6% of exp
    experiment:                        0.072072 Ha = 1.961 eV

If the FCI ceiling climbs toward 1.961 eV as functions are added, the Paper 58
attribution is confirmed and basis size is the live accuracy lever. If it
plateaus well short, the deficit is something else and the attribution needs
revisiting.

DESIGN NOTES.

* Na zetas are NOT re-optimized. They were fitted on the isolated Na atom
  (1s 10.63, 2s 3.3, 2p 3.44, 3s 0.836) and adding an H-centred function does
  not change the Na atom. Reusing them preserves the "fragment prediction, not
  molecular fit" property of the original probe.

* The diffuse-H exponent is NOT variationally determined. H- is barely bound
  (EA = 0.0277 Ha) and in a small basis its variational minimum runs away to
  zeta -> 0, so an atom-only optimization is not available. Instead the driver
  SCANS a grid of diffuse exponents and reports the whole curve. Consequence,
  disclosed: the best-of-grid D_e is an upper envelope over that one parameter,
  NOT a parameter-free prediction like the M=7 baseline. Read the curve, not
  the maximum.

* Atom references are recomputed in each extended basis, so a function that
  helps the free atom as much as the molecule shows no spurious gain. There is
  still no counterpoise correction, so every D_e here (baseline included) is
  BSSE-inflated. BSSE correction would push D_e DOWN -- it is a disclosure item,
  not a route to closing an underbinding gap.

* Dense FCI (geovac.noci_engine.fci_ground) caps out around dim ~2000, i.e.
  M=8 for 12 electrons. Beyond that the driver freezes the Na 1s/2s/2p core
  (5 orbitals, 10 electrons), leaving 2 active electrons. Two gates validate
  this before it is used (n_core=0 identity; frozen-vs-all-electron D_e at M=7).

Exploratory. No paper claim, no CHANGELOG entry, no test.
"""

from __future__ import annotations

import json
import os
import sys
import time
from itertools import combinations
from typing import Dict, List, Sequence, Tuple

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from geovac.noci_engine import (
    HARTREE_TO_EV,
    STO3G_2S,
    STO6G_1S,
    BasisFn,
    det_pair_gensc,
    fci_ground,
    fit_sto_shape,
    integral_set_md,
    lowdin_orbitals,
    noci_ground_gensc,
    transform_integrals,
)

UP, DN = 0, 1
R_EQ_EXP = 3.566            # a0  (1.887 A, NIST CCCBDB)
D_E_EXP_EV = 1.961          # eV  (15,815 cm^-1, Huang et al. JCP 133 044301)
D_E_EXP_HA = D_E_EXP_EV / HARTREE_TO_EV

Z_NA = {"1s": 10.63, "2s": 3.3, "2p": 3.44, "3s": 0.836}   # stored optimum
Z_H = 1.0

# stored M=7 baseline, for the reproduction gate
BASE_FCI_DE = 0.04318852868348699
BASE_NOCI3_DE = 0.039357282832042984
BASE_FCI_REQ = 3.5945536011363233


# ------------------------------------------------------------------ shapes

def build_shapes() -> Dict[str, Tuple[np.ndarray, np.ndarray]]:
    """STO shapes as Gaussian expansions. 1s/2s hardcoded (N1/N2 lineage);
    the rest least-squares fitted, quality printed."""
    shapes = {"1s": STO6G_1S, "2s": STO3G_2S}
    for kind, (l, n_r) in (("2p", (1, 2)), ("3s", (0, 3)), ("3p", (1, 3))):
        a, dco, q = fit_sto_shape(l, n_r)
        print(f"[shapes] fitted {kind:3s}: <fit|STO> = {q:.6f}")
        shapes[kind] = (a, dco)
    return shapes


# ------------------------------------------------------------------ bases
#
# A "basis spec" is a list of (label, center_tag, shape_key, cartesian_lmn,
# zeta). center_tag is 'Na' or 'H'. Building it this way keeps the atom
# references and the molecule on literally the same function list.

def base_spec() -> List[tuple]:
    return [
        ("Na1s",  "Na", "1s", (0, 0, 0), Z_NA["1s"]),
        ("Na2s",  "Na", "2s", (0, 0, 0), Z_NA["2s"]),
        ("Na2px", "Na", "2p", (1, 0, 0), Z_NA["2p"]),
        ("Na2py", "Na", "2p", (0, 1, 0), Z_NA["2p"]),
        ("Na2pz", "Na", "2p", (0, 0, 1), Z_NA["2p"]),
        ("Na3s",  "Na", "3s", (0, 0, 0), Z_NA["3s"]),
        ("H1s",   "H",  "1s", (0, 0, 0), Z_H),
    ]


def ladder_specs(z_diff: float) -> Dict[str, List[tuple]]:
    """Basis ladder. Each rung adds the next physically-motivated shell."""
    b0 = base_spec()
    b1 = b0 + [("Hdiff", "H", "1s", (0, 0, 0), z_diff)]
    b2 = b1 + [("H2pz", "H", "2p", (0, 0, 1), Z_H)]
    b3 = b2 + [("Na3pz", "Na", "3p", (0, 0, 1), Z_NA["3s"])]
    b4 = b2 + [("H2px", "H", "2p", (1, 0, 0), Z_H),
               ("H2py", "H", "2p", (0, 1, 0), Z_H)]
    b5 = b4 + [("Na3pz", "Na", "3p", (0, 0, 1), Z_NA["3s"]),
               ("Na3px", "Na", "3p", (1, 0, 0), Z_NA["3s"]),
               ("Na3py", "Na", "3p", (0, 1, 0), Z_NA["3s"])]
    return {
        "B0 base (M=7)": b0,
        "B1 +H diffuse s (M=8)": b1,
        "B2 +H 2pz (M=9)": b2,
        "B3 +Na 3pz (M=10)": b3,
        "B4 +H 2px,2py (M=11)": b4,
        "B5 +Na 3p full (M=14)": b5,
    }


def build_orbs(spec, shapes, pos_na, pos_h) -> List[BasisFn]:
    centers = {"Na": pos_na, "H": pos_h}
    out = []
    for _lab, tag, key, lmn, zeta in spec:
        a, d = shapes[key]
        out.append(BasisFn(centers[tag], lmn, a * zeta ** 2, d))
    return out


def spec_subset(spec, tag) -> List[tuple]:
    return [row for row in spec if row[1] == tag]


# ------------------------------------------------------- frozen-core helper

def frozen_core(h: np.ndarray, g: np.ndarray, n_core: int):
    """Return (E_core, h_eff, g_active) for n_core doubly-occupied orbitals.

    Orthonormal spatial basis, g in chemist notation g[p,q,r,s] = (pq|rs).
    n_core = 0 must be the identity map (gate G-FC1).
    """
    if n_core == 0:
        return 0.0, h, g
    c = slice(0, n_core)
    a = slice(n_core, h.shape[0])

    e_core = 2.0 * np.trace(h[c, c])
    gc = g[c, c, c, c]
    e_core += 2.0 * np.einsum("iijj->", gc) - np.einsum("ijji->", gc)

    h_eff = h[a, a].copy()
    h_eff += 2.0 * np.einsum("pqii->pq", g[a, a, c, c])
    h_eff -= np.einsum("piiq->pq", g[a, c, c, a])
    return float(e_core), h_eff, g[a, a, a, a]


def fci_energy(h, g, n_elec, n_core: int = 0) -> float:
    e_core, h_eff, g_act = frozen_core(h, g, n_core)
    return e_core + fci_ground(h_eff, g_act, n_elec - 2 * n_core)


# ------------------------------------------------------------------ pieces

def atom_energy(spec, shapes, tag: str, z_nuc: float, n_elec: int,
                n_core: int) -> float:
    """Isolated-atom in-basis FCI, using only that atom's own functions."""
    orig = np.array([0.0, 0.0, 0.0])
    sub = spec_subset(spec, tag)
    orbs = build_orbs(sub, shapes, orig, orig)
    s, h, g = integral_set_md(orbs, [(orig, z_nuc)])
    x = lowdin_orbitals(s)
    ht, gt = transform_integrals(x, h, g)
    return fci_energy(ht, gt, n_elec, n_core)


def atom_energy_cp(spec, shapes, tag: str, r: float, n_core: int) -> float:
    """Counterpoise atom energy: the atom in the FULL molecular basis.

    All basis functions sit at their molecular positions, but only THIS atom's
    nucleus is present (the partner's charge is set to zero -- "ghost"
    functions). This is the Boys-Bernardi counterpoise reference.

    Why it is needed here. Without it, the molecule can borrow the partner's
    basis functions to lower its own energy while the isolated-atom reference
    cannot, so D_e is inflated by basis-set superposition error. BSSE GROWS
    WITH BASIS SIZE -- which is precisely the axis this driver varies -- so an
    uncorrected D_e ladder cannot distinguish real basis improvement from
    increased borrowing. Every D_e in the uncorrected table is contaminated
    this way, the baseline included.
    """
    pos_na = np.array([0.0, 0.0, 0.0])
    pos_h = np.array([0.0, 0.0, r])
    orbs = build_orbs(spec, shapes, pos_na, pos_h)
    if tag == "Na":
        nuclei, n_elec = [(pos_na, 11.0)], 11
    else:
        nuclei, n_elec = [(pos_h, 1.0)], 1
        n_core = 0                      # H has no core to freeze
    s, h, g = integral_set_md(orbs, nuclei)
    x = lowdin_orbitals(s)
    ht, gt = transform_integrals(x, h, g)
    return fci_energy(ht, gt, n_elec, n_core)


def molecule_energies(spec, shapes, r: float, n_core: int,
                      dets: Sequence = None):
    """Return (E_fci, E_noci3, s_min) at separation r, both incl. V_nn.

    s_min is the smallest eigenvalue of the AO overlap. It is returned because
    near-linear-dependence in the basis makes the Loewdin S^{-1/2} transform
    ill-conditioned and the resulting "energies" non-variational -- which is a
    way to manufacture a spuriously good D_e by adding a function that nearly
    duplicates one already present. See the LINDEP gate.
    """
    pos_na = np.array([0.0, 0.0, 0.0])
    pos_h = np.array([0.0, 0.0, r])
    orbs = build_orbs(spec, shapes, pos_na, pos_h)
    s, h, g = integral_set_md(orbs, [(pos_na, 11.0), (pos_h, 1.0)])
    vnn = 11.0 / r
    s_min = float(np.linalg.eigvalsh(s).min())

    x = lowdin_orbitals(s)
    ht, gt = transform_integrals(x, h, g)
    e_fci = fci_energy(ht, gt, 12, n_core) + vnn

    e_noci = None
    if dets is not None:
        e_noci = noci_ground_gensc(dets, s, h, g)[0] + vnn
    return e_fci, e_noci, s_min


def noci3_dets(spec) -> List[List[Tuple[int, int]]]:
    """The published 3-determinant ladder: 2 covalent + Na+H-.

    Indices are looked up by label so the pattern survives basis extension;
    the ADDED orbitals stay empty in these three determinants, which is the
    point -- it isolates the effect of orbital-shape quality on a fixed
    determinant pattern.
    """
    lab = [row[0] for row in spec]
    core_lab = ["Na1s", "Na2s", "Na2px", "Na2py", "Na2pz"]
    core = [(lab.index(o), sp) for o in core_lab for sp in (UP, DN)]
    na3s, h1s = lab.index("Na3s"), lab.index("H1s")
    return [
        core + [(na3s, UP), (h1s, DN)],
        core + [(h1s, UP), (na3s, DN)],
        core + [(h1s, UP), (h1s, DN)],
    ]


# ------------------------------------------------------------------ well

R_GRID = [2.6, 2.9, 3.1, 3.3, 3.45, 3.566, 3.7, 3.85, 4.0, 4.3,
          4.7, 5.2, 6.0, 7.5, 10.0]


def well(rs: List[float], es: List[float], e_diss: float):
    """Minimum by parabolic refinement through the three lowest points."""
    es_a = np.array(es)
    k = int(np.argmin(es_a))
    if k in (0, len(es) - 1):
        return {"R_eq": rs[k], "E_min": es[k],
                "D_e_Ha": e_diss - es[k], "edge": True}
    x = np.array(rs[k - 1:k + 2])
    y = es_a[k - 1:k + 2]
    c = np.polyfit(x, y, 2)
    r_eq = -c[1] / (2 * c[0])
    e_min = np.polyval(c, r_eq)
    return {"R_eq": float(r_eq), "E_min": float(e_min),
            "D_e_Ha": float(e_diss - e_min), "edge": False}


LINDEP_TOL = 1e-6      # standard threshold on the smallest AO overlap eigenvalue


def scan(spec, shapes, n_core: int, with_noci: bool = True):
    """Scan the PES.

    NOTE ON THE NOCI COLUMN. `noci_ground_gensc` works on the raw AO integrals
    with explicit ALL-ELECTRON determinants, so it cannot be mixed with a
    frozen-core dissociation reference -- doing so subtracts a 2-electron
    reference from a 12-electron energy. NOCI is therefore only computed on the
    all-electron path (n_core = 0). The first version of this driver got this
    wrong and printed NOCI D_e values ABOVE the FCI ceiling, which is
    variationally impossible and was the tell.
    """
    if with_noci and n_core != 0:
        raise ValueError("NOCI is all-electron; cannot pair with frozen core")

    e_na = atom_energy(spec, shapes, "Na", 11.0, 11, n_core)
    e_h = atom_energy(spec, shapes, "H", 1.0, 1, 0)
    e_diss = e_na + e_h

    dets = noci3_dets(spec) if with_noci else None
    fci, noci, smins = [], [], []
    for r in R_GRID:
        ef, en, sm = molecule_energies(spec, shapes, r, n_core, dets)
        fci.append(ef)
        smins.append(sm)
        if en is not None:
            noci.append(en)

    out = {"E_Na": e_na, "E_H": e_h, "E_diss": e_diss,
           "s_min": min(smins), "lindep": min(smins) < LINDEP_TOL,
           "fci": well(R_GRID, fci, e_diss)}
    if noci:
        out["noci3"] = well(R_GRID, noci, e_diss)
    return out


def scan_cp(spec, shapes, n_core: int):
    """Counterpoise-corrected PES.

    The interaction energy is formed pointwise against ghost-basis atom
    references at the SAME R:

        dE_CP(R) = E_mol(R) - E_Na(ghost H, R) - E_H(ghost Na, R)

    so the borrowing available to the molecule is also available to each
    reference and cancels. D_e = -min_R dE_CP(R). Also returns the raw
    (uncorrected) curve so the BSSE magnitude is visible.
    """
    e_na_iso = atom_energy(spec, shapes, "Na", 11.0, 11, n_core)
    e_h_iso = atom_energy(spec, shapes, "H", 1.0, 1, 0)
    e_diss_iso = e_na_iso + e_h_iso

    raw, cp, smins = [], [], []
    for r in R_GRID:
        e_mol, _, sm = molecule_energies(spec, shapes, r, n_core, None)
        e_na = atom_energy_cp(spec, shapes, "Na", r, n_core)
        e_h = atom_energy_cp(spec, shapes, "H", r, n_core)
        raw.append(e_mol)
        cp.append(e_mol - e_na - e_h + e_diss_iso)   # shifted to common zero
        smins.append(sm)

    return {"E_diss": e_diss_iso, "s_min": min(smins),
            "raw": well(R_GRID, raw, e_diss_iso),
            "cp": well(R_GRID, cp, e_diss_iso),
            "bsse_at_min_Ha": float(np.min(np.array(raw))
                                    - np.min(np.array(cp)))}


def pct(d_ha: float) -> float:
    return 100.0 * d_ha / D_E_EXP_HA


# ------------------------------------------------------------------ main

def main() -> None:
    print("=== NaH basis extension: does the in-basis FCI ceiling rise? ===")
    print(f"    experiment: R_eq = {R_EQ_EXP} a0, D_e = {D_E_EXP_EV} eV "
          f"({D_E_EXP_HA:.6f} Ha)\n")
    shapes = build_shapes()
    results = {}

    # ---------------- Gate G-B0: reproduce the stored baseline -------------
    print("\n[G-B0] reproducing the stored M=7 baseline (all-electron)")
    t0 = time.time()
    b0 = base_spec()
    r0 = scan(b0, shapes, n_core=0, with_noci=True)
    d_fci = abs(r0["fci"]["D_e_Ha"] - BASE_FCI_DE)
    d_noci = abs(r0["noci3"]["D_e_Ha"] - BASE_NOCI3_DE)
    print(f"       FCI   D_e = {r0['fci']['D_e_Ha']:.6f} Ha  "
          f"(stored {BASE_FCI_DE:.6f}, delta {d_fci:.2e})")
    print(f"       NOCI3 D_e = {r0['noci3']['D_e_Ha']:.6f} Ha  "
          f"(stored {BASE_NOCI3_DE:.6f}, delta {d_noci:.2e})")
    print(f"       R_eq(FCI) = {r0['fci']['R_eq']:.4f} a0 "
          f"(stored {BASE_FCI_REQ:.4f})")
    # R grid differs from the original probe, so agreement is to grid
    # resolution, not to machine precision.
    gate_b0 = d_fci < 2e-3 and d_noci < 2e-3
    print(f"       gate (both within 2 mHa): "
          f"{'PASS' if gate_b0 else 'FAIL'}   [{time.time()-t0:.0f} s]")
    results["B0_allelectron"] = r0

    # ---------------- Gate G-FC1: frozen_core(0) is the identity ----------
    print("\n[G-FC1] frozen_core(n_core=0) identity check")
    pos_na, pos_h = np.array([0.0, 0.0, 0.0]), np.array([0.0, 0.0, 3.566])
    orbs = build_orbs(b0, shapes, pos_na, pos_h)
    s, h, g = integral_set_md(orbs, [(pos_na, 11.0), (pos_h, 1.0)])
    x = lowdin_orbitals(s)
    ht, gt = transform_integrals(x, h, g)
    e_all = fci_energy(ht, gt, 12, 0)
    ec, he, ge = frozen_core(ht, gt, 0)
    ident = (ec == 0.0 and he is ht and ge is gt)
    print(f"       identity map: {'PASS' if ident else 'FAIL'}")

    # ---------------- Gate G-FC2: frozen core preserves D_e ---------------
    print("\n[G-FC2] frozen Na 1s/2s/2p (n_core=5) vs all-electron, M=7")
    t0 = time.time()
    r0fc = scan(b0, shapes, n_core=5, with_noci=False)
    dd = r0fc["fci"]["D_e_Ha"] - r0["fci"]["D_e_Ha"]
    print(f"       all-electron FCI D_e = {r0['fci']['D_e_Ha']:.6f} Ha")
    print(f"       frozen-core  FCI D_e = {r0fc['fci']['D_e_Ha']:.6f} Ha")
    print(f"       difference           = {dd:+.6f} Ha "
          f"({100*abs(dd)/r0['fci']['D_e_Ha']:.2f}% of D_e)")
    gate_fc = abs(dd) < 5e-3
    print(f"       gate (|delta D_e| < 5 mHa): "
          f"{'PASS' if gate_fc else 'FAIL'}   [{time.time()-t0:.0f} s]")
    results["B0_frozencore"] = r0fc

    if not (gate_b0 and ident and gate_fc):
        print("\n*** GATES FAILED -- ladder not run. ***")
        return

    # ---------------- diffuse-exponent scan on the B1 rung ----------------
    print("\n[scan] diffuse-H exponent (B1, M=8, all-electron)")
    print("       a diffuse zeta close to the H 1s zeta (1.0) nearly duplicates")
    print("       it; s_min is the linear-dependence monitor. Rows flagged")
    print("       LINDEP have an ill-conditioned Loewdin transform and their")
    print("       D_e is NOT variationally meaningful.")
    print(f"       {'zeta':>6} {'FCI D_e (Ha)':>14} {'eV':>8} {'% exp':>7} "
          f"{'R_eq':>7} {'s_min':>10}  flag")
    diff_rows = []
    best = None
    for zd in (0.10, 0.15, 0.20, 0.25, 0.30, 0.40, 0.55, 0.75):
        spec = ladder_specs(zd)["B1 +H diffuse s (M=8)"]
        rr = scan(spec, shapes, n_core=0, with_noci=False)
        de = rr["fci"]["D_e_Ha"]
        flag = "LINDEP" if rr["lindep"] else ""
        print(f"       {zd:6.2f} {de:14.6f} {de*HARTREE_TO_EV:8.3f} "
              f"{pct(de):7.1f} {rr['fci']['R_eq']:7.3f} {rr['s_min']:10.2e}"
              f"  {flag}")
        diff_rows.append({"zeta": zd, "D_e_Ha": de, "s_min": rr["s_min"],
                          "lindep": rr["lindep"], "R_eq": rr["fci"]["R_eq"]})
        if not rr["lindep"] and (best is None or de > best[1]):
            best = (zd, de)
    z_diff = best[0]
    print(f"       -> best CLEAN-conditioned zeta_diffuse = {z_diff}")
    print("       (still an upper envelope over this one parameter, not a")
    print("        parameter-free prediction like the M=7 baseline)")
    results["diffuse_scan"] = diff_rows

    # ---------------- the ladder ------------------------------------------
    print(f"\n[ladder] FCI ceiling, frozen core (n_core=5), "
          f"zeta_diffuse = {z_diff}")
    print(f"  {'basis':<24} {'M':>3} {'FCI D_e':>9} {'eV':>7} {'% exp':>6} "
          f"{'R_eq':>7} {'s_min':>10}  flag")
    specs = ladder_specs(z_diff)
    for name, spec in specs.items():
        t0 = time.time()
        rr = scan(spec, shapes, n_core=5, with_noci=False)
        de = rr["fci"]["D_e_Ha"]
        flag = "LINDEP" if rr["lindep"] else ("edge" if rr["fci"]["edge"] else "")
        print(f"  {name:<24} {len(spec):>3} {de:9.6f} "
              f"{de*HARTREE_TO_EV:7.3f} {pct(de):6.1f} "
              f"{rr['fci']['R_eq']:7.3f} {rr['s_min']:10.2e}  {flag}"
              f"   [{time.time()-t0:.0f} s]")
        results[name] = rr

    # NOCI compactness cross-check, all-electron only (see scan() docstring).
    print("\n[noci] does 3-determinant compactness survive a basis addition?")
    print("       (all-electron; the 3 dets leave added orbitals EMPTY, so this")
    print("        isolates orbital-shape quality at fixed determinant count)")
    for name in ("B0 base (M=7)", "B1 +H diffuse s (M=8)"):
        rr = scan(specs[name], shapes, n_core=0, with_noci=True)
        dn, df = rr["noci3"]["D_e_Ha"], rr["fci"]["D_e_Ha"]
        print(f"  {name:<24} NOCI3 {dn*HARTREE_TO_EV:6.3f} eV / "
              f"FCI {df*HARTREE_TO_EV:6.3f} eV = {100*dn/df:5.1f}% "
              f"{'  <-- ABOVE FCI, IMPOSSIBLE' if dn > df + 1e-9 else ''}")
        results[name + " [AE+noci]"] = rr

    # ---------------- verdict ---------------------------------------------
    print("\n=== where we land ===")
    base_ev = results["B0 base (M=7)"]["fci"]["D_e_Ha"] * HARTREE_TO_EV
    top = max(results[n]["fci"]["D_e_Ha"] for n in specs) * HARTREE_TO_EV
    print(f"  FCI ceiling, minimal basis : {base_ev:.3f} eV "
          f"({100*base_ev/D_E_EXP_EV:.1f}% of experiment)")
    print(f"  FCI ceiling, best rung     : {top:.3f} eV "
          f"({100*top/D_E_EXP_EV:.1f}% of experiment)")
    print(f"  experiment                 : {D_E_EXP_EV:.3f} eV")
    print(f"  gap closed by basis        : "
          f"{100*(top-base_ev)/(D_E_EXP_EV-base_ev):.1f}% of the deficit")

    out = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                       "data", "noci_nah_basis_extension.json")
    with open(out, "w") as fh:
        json.dump({"exp": {"R_eq_a0": R_EQ_EXP, "D_e_eV": D_E_EXP_EV},
                   "z_na": Z_NA, "z_diffuse": z_diff,
                   "R_grid": R_GRID, "results": results}, fh, indent=2)
    print(f"\n  written: {out}")


def main_cp() -> None:
    """Counterpoise-corrected ladder -- the only version whose D_e trend is
    interpretable, since BSSE grows with basis size."""
    print("=== NaH basis extension, COUNTERPOISE-CORRECTED ===")
    print(f"    experiment: R_eq = {R_EQ_EXP} a0, D_e = {D_E_EXP_EV} eV\n")
    print("    raw  = atom references in atom-only bases (BSSE-inflated,")
    print("           this is what the stored M=7 baseline reports)")
    print("    CP   = Boys-Bernardi ghost-basis references at the same R\n")
    shapes = build_shapes()

    z_diff = float(sys.argv[2]) if len(sys.argv) > 2 else 0.30
    specs = ladder_specs(z_diff)
    print(f"[ladder] frozen core (n_core=5), zeta_diffuse = {z_diff}\n")
    print(f"  {'basis':<24} {'M':>3} | {'raw eV':>7} {'%exp':>5} {'R_eq':>6} "
          f"| {'CP eV':>7} {'%exp':>5} {'R_eq':>6} | {'BSSE eV':>8}")
    results = {}
    for name, spec in specs.items():
        t0 = time.time()
        rr = scan_cp(spec, shapes, n_core=5)
        raw, cp = rr["raw"], rr["cp"]
        print(f"  {name:<24} {len(spec):>3} | "
              f"{raw['D_e_Ha']*HARTREE_TO_EV:7.3f} {pct(raw['D_e_Ha']):5.1f} "
              f"{raw['R_eq']:6.3f} | "
              f"{cp['D_e_Ha']*HARTREE_TO_EV:7.3f} {pct(cp['D_e_Ha']):5.1f} "
              f"{cp['R_eq']:6.3f} | "
              f"{rr['bsse_at_min_Ha']*HARTREE_TO_EV:8.3f}"
              f"   [{time.time()-t0:.0f} s]")
        results[name] = rr

    print("\n=== where we land (counterpoise-corrected) ===")
    names = list(specs)
    lo = results[names[0]]["cp"]["D_e_Ha"] * HARTREE_TO_EV
    hi = max(results[n]["cp"]["D_e_Ha"] for n in names) * HARTREE_TO_EV
    print(f"  CP D_e, minimal basis : {lo:.3f} eV ({100*lo/D_E_EXP_EV:.1f}%)")
    print(f"  CP D_e, best rung     : {hi:.3f} eV ({100*hi/D_E_EXP_EV:.1f}%)")
    print(f"  experiment            : {D_E_EXP_EV:.3f} eV")
    if D_E_EXP_EV - lo > 1e-9:
        print(f"  deficit closed        : "
              f"{100*(hi-lo)/(D_E_EXP_EV-lo):.1f}%")

    out = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                       "data", "noci_nah_basis_extension_cp.json")
    with open(out, "w") as fh:
        json.dump({"exp": {"R_eq_a0": R_EQ_EXP, "D_e_eV": D_E_EXP_EV},
                   "z_diffuse": z_diff, "R_grid": R_GRID,
                   "results": results}, fh, indent=2)
    print(f"\n  written: {out}")


if __name__ == "__main__":
    if len(sys.argv) > 1 and sys.argv[1] == "--cp":
        main_cp()
    else:
        main()
