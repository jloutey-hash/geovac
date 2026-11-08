"""Flesh out the atomic xTC angular-sparsity result across the whole monoatomic
library (Z = 1..56).

The tracked engine ``geovac.xtc_angular_sparsity.fast_angular`` returns EXACT
integer block counts (angular support is radial-independent), and the fill-in is
a function ONLY of the reference density's angular multipoles.  So the monoatomic
sweep reduces to a sweep over the open-shell angular character of each atom's
ground-state valence configuration.

Hypothesis under test (predict-before-compute, per the W3 protocol):
  xTC fill-in = 0  <=>  the reference density is SPHERICAL
                   <=>  (Unsoeld) the open subshell is closed, half-filled
                        high-spin, or s-type.
  Non-spherical p  -> refmult {(0,0),(2,0)}   -> the paper's 0/4/96 pattern.
  Non-spherical d  -> refmult {(0,0),(2,0),(4,0)} -> NEW (transition metals).

Prediction written BEFORE running (see PREDICTION dict).  No f-shell atoms are in
scope (Z<=56 stops at Ba; La is Z=57), so only s/p/d open shells occur.
"""
from __future__ import annotations

import json
from typing import Dict, List, Tuple

from geovac.xtc_angular_sparsity import fast_angular, reference_density_multipoles

LM = Tuple[int, int]

# --- pre-registered prediction (written before the run) --------------------
PREDICTION = {
    "spherical_classes_zero_fillin": True,   # closed / half-filled-hs / s-open -> 0
    "p_open_pattern_lmax_1_2_3": [0, 4, 96],  # matches paper tab:xtc_fillin
    "d_open_has_L4_multipole": True,          # refmult contains (4,0)
    "m_filling_independent": True,            # support depends only on (l, spherical?)
    "core_support_invariant": True,           # spherical core adds only (0,0)
}


# ---------------------------------------------------------------------------
# high-spin (Hund) open-subshell angular occupation
# ---------------------------------------------------------------------------
def open_subshell(l: int, k: int) -> List[LM]:
    """(l,m) spatial-orbital occupation for k electrons in an l-subshell, high spin.

    Fill each m once (spin up) in the order 0, +1, -1, +2, -2, ... then pair.
    Only the SET of occupied spatial orbitals (with multiplicity) matters for the
    diagonal density multipoles, and the specific m-order does not affect angular
    SUPPORT -- verified explicitly below.
    """
    ms = [0]
    for a in range(1, l + 1):
        ms += [a, -a]           # 0, +1, -1, +2, -2, ...
    occ: List[LM] = []
    # first pass: one electron per m (singly occupied)
    for i in range(min(k, 2 * l + 1)):
        occ.append((l, ms[i]))
    # second pass: pair up the rest
    for i in range(k - (2 * l + 1)):
        occ.append((l, ms[i]))
    return occ


def is_spherical(l: int, k: int) -> bool:
    """Unsoeld: closed (k=0 or 4l+2), half-filled (k=2l+1), or s (l=0) -> spherical."""
    return l == 0 or k == 0 or k == (2 * l + 1) or k == (4 * l + 2)


# ---------------------------------------------------------------------------
# ground-state open (valence) subshell for Z = 1..56
# (l, k) of the shell that determines the angular character; None -> closed/s only
# Uses standard ground-state configurations; the s-shells are always spherical.
# ---------------------------------------------------------------------------
def ground_open(Z: int) -> Tuple[int, int, str]:
    """Return (l, k, label) of the angular-determining open subshell.

    For angular support only the NON-s open subshell matters (s density is
    spherical).  Atoms whose only open shell is s (H, alkali, alkaline-earth,
    noble) return (0, k, ...) which is_spherical() -> True.
    """
    # (symbol, angular-determining subshell l, occupation k)
    TABLE: Dict[int, Tuple[str, int, int]] = {
        1: ("H", 0, 1), 2: ("He", 0, 2),
        3: ("Li", 0, 1), 4: ("Be", 0, 2),
        5: ("B", 1, 1), 6: ("C", 1, 2), 7: ("N", 1, 3), 8: ("O", 1, 4),
        9: ("F", 1, 5), 10: ("Ne", 1, 6),
        11: ("Na", 0, 1), 12: ("Mg", 0, 2),
        13: ("Al", 1, 1), 14: ("Si", 1, 2), 15: ("P", 1, 3), 16: ("S", 1, 4),
        17: ("Cl", 1, 5), 18: ("Ar", 1, 6),
        19: ("K", 0, 1), 20: ("Ca", 0, 2),
        21: ("Sc", 2, 1), 22: ("Ti", 2, 2), 23: ("V", 2, 3), 24: ("Cr", 2, 5),
        25: ("Mn", 2, 5), 26: ("Fe", 2, 6), 27: ("Co", 2, 7), 28: ("Ni", 2, 8),
        29: ("Cu", 2, 10), 30: ("Zn", 2, 10),
        31: ("Ga", 1, 1), 32: ("Ge", 1, 2), 33: ("As", 1, 3), 34: ("Se", 1, 4),
        35: ("Br", 1, 5), 36: ("Kr", 1, 6),
        37: ("Rb", 0, 1), 38: ("Sr", 0, 2),
        39: ("Y", 2, 1), 40: ("Zr", 2, 2), 41: ("Nb", 2, 4), 42: ("Mo", 2, 5),
        43: ("Tc", 2, 5), 44: ("Ru", 2, 7), 45: ("Rh", 2, 8), 46: ("Pd", 2, 10),
        47: ("Ag", 2, 10), 48: ("Cd", 2, 10),
        49: ("In", 1, 1), 50: ("Sn", 1, 2), 51: ("Sb", 1, 3), 52: ("Te", 1, 4),
        53: ("I", 1, 5), 54: ("Xe", 1, 6),
        55: ("Cs", 0, 1), 56: ("Ba", 0, 2),
    }
    sym, l, k = TABLE[Z]
    return l, k, sym


def classify(l: int, k: int) -> str:
    if is_spherical(l, k):
        return "SPHERICAL"
    return {1: "P-OPEN", 2: "D-OPEN"}[l]


# ---------------------------------------------------------------------------
# run
# ---------------------------------------------------------------------------
def run_class(ref_ang: List[LM], label: str) -> Dict:
    rows = {}
    for lmax in (1, 2, 3):
        r = fast_angular(lmax, ref_ang, Lmax=2 * lmax)
        rows[lmax] = dict(
            n_coulomb=r["n_coulomb"], n_l3=r["n_l3"], n_fill=r["n_fill_in"],
            l3_density=round(r["l3_density"], 4), coul_density=round(r["coul_density"], 4),
        )
    refmult = reference_density_multipoles(ref_ang, 6)
    return dict(label=label, refmult=sorted(refmult.keys()), rows=rows)


def main() -> None:
    print("=" * 78)
    print("PREDICTION (pre-registered):", json.dumps(PREDICTION))
    print("=" * 78)

    # --- representative occupations for the three classes -------------------
    reps = {
        "SPHERICAL(s2)":  open_subshell(0, 2),      # e.g. Be/Mg/Ca valence
        "SPHERICAL(p3hs)": open_subshell(1, 3),     # half-filled p (N, P) -- test Unsoeld
        "SPHERICAL(d5hs)": open_subshell(2, 5),     # half-filled d (Mn, Cr) -- test Unsoeld
        "P-OPEN(p2)":     open_subshell(1, 2),      # C/Si
        "P-OPEN(p4)":     open_subshell(1, 4),      # O/S -- same class as p2?
        "D-OPEN(d2)":     open_subshell(2, 2),      # Ti/Zr  -- NEW
        "D-OPEN(d1)":     open_subshell(2, 1),      # Sc/Y   -- NEW
        "D-OPEN(d8)":     open_subshell(2, 8),      # Ni/Rh  -- NEW
    }

    print("\n--- CLASS TABLE (exact integer fill-in; basis = s+p / s+p+d / s+p+d+f) ---")
    print(f"{'class(rep occ)':<18} {'refmult(L,M)':<24} "
          f"{'fill@sp':>8} {'fill@spd':>9} {'fill@spdf':>10} {'dens@spdf':>10}")
    results = {}
    for label, occ in reps.items():
        rc = run_class(occ, label)
        results[label] = rc
        rm = ",".join(f"({a},{b})" for a, b in rc["refmult"])
        f1 = rc["rows"][1]["n_fill"]; f2 = rc["rows"][2]["n_fill"]; f3 = rc["rows"][3]["n_fill"]
        dens = rc["rows"][3]["l3_density"]
        print(f"{label:<18} {rm:<24} {f1:>8} {f2:>9} {f3:>10} {dens:>10}")

    # --- validations --------------------------------------------------------
    print("\n--- VALIDATIONS ---")
    # (1) m-filling independence: p2 as (0,+1) vs (0,0)-paired
    from geovac.xtc_angular_sparsity import REF_ANG
    a = fast_angular(2, REF_ANG["C_2p2"], Lmax=4)["n_fill_in"]
    b = fast_angular(2, REF_ANG["C_2p0sq"], Lmax=4)["n_fill_in"]
    print(f"  m-filling independence (C 2p^2 two fillings): {a} == {b} -> {a == b}")
    # (2) half-filled -> spherical -> 0 fill-in
    p3 = fast_angular(3, open_subshell(1, 3), Lmax=6)["n_fill_in"]
    d5 = fast_angular(3, open_subshell(2, 5), Lmax=6)["n_fill_in"]
    print(f"  half-filled p^3 fill-in @spdf: {p3} (predict 0);  d^5: {d5} (predict 0)")
    # (3) core support-invariance: p2 alone vs p2 + spherical core
    core = [(0, 0)] * 6 + [(1, -1), (1, -1), (1, 0), (1, 0), (1, 1), (1, 1)]  # [Ne]-like
    bare = fast_angular(2, open_subshell(1, 2), Lmax=4)["n_fill_in"]
    cored = fast_angular(2, core + open_subshell(1, 2), Lmax=4)["n_fill_in"]
    print(f"  core support-invariance (p^2 bare vs +[Ne] core): {bare} == {cored} -> {bare == cored}")
    # (4) p2 vs p4 same class
    p2 = fast_angular(3, open_subshell(1, 2), Lmax=6)["n_fill_in"]
    p4 = fast_angular(3, open_subshell(1, 4), Lmax=6)["n_fill_in"]
    print(f"  p^2 vs p^4 fill-in @spdf: {p2} == {p4} -> {p2 == p4}")

    # --- per-atom periodic sweep -------------------------------------------
    print("\n--- PER-ATOM SWEEP (Z = 1..56), fill-in @ s+p+d+f ---")
    fill_by_class = {"SPHERICAL": 0, "P-OPEN": p2, "D-OPEN": results["D-OPEN(d2)"]["rows"][3]["n_fill"]}
    per_atom = {}
    counts = {"SPHERICAL": [], "P-OPEN": [], "D-OPEN": []}
    for Z in range(1, 57):
        l, k, sym = ground_open(Z)
        cls = classify(l, k)
        per_atom[Z] = dict(sym=sym, l=l, k=k, cls=cls, fill_spdf=fill_by_class[cls])
        counts[cls].append(sym)
    for cls in ("SPHERICAL", "P-OPEN", "D-OPEN"):
        syms = counts[cls]
        print(f"  {cls:<10} (n={len(syms):>2}, fill@spdf={fill_by_class[cls]:>3}): {' '.join(syms)}")

    out = dict(prediction=PREDICTION, classes=results, per_atom=per_atom,
               fill_by_class=fill_by_class,
               class_counts={c: len(s) for c, s in counts.items()})
    with open("debug/data/xtc_monoatomic_sweep.json", "w") as fh:
        json.dump(out, fh, indent=2)
    print("\nwrote debug/data/xtc_monoatomic_sweep.json")


if __name__ == "__main__":
    main()
