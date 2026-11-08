"""Axis check: is the TC three-body non-collapse the SAME obstruction as >=3-projection
wildness, or a shape match?  (2026-08-25, /aha follow-on 1.)

THE /aha PROPOSAL.  Both look like "pairwise reduces, triple doesn't":
  - TC (v5.0.2): the plane-wave two-body collapse is ABELIAN (two correlator momenta at a
    vertex add to ONE resultant); GeoVac's Gaunt/6j analog is NON-ABELIAN (Y.Y = sum_Lambda,
    a CG multiplet), so the three-body operator does not collapse.
  - Projections (v5.1.0): two center-projections are tame (Halmos), three are *-wild and
    the molecular triple is irreducible M_6.

WHY IT FAILS, in two independent ways.

(a) MISREAD NUMBERS.  The TC memo's "rank 4 -> 9 -> 16" runs over the columns
    l=1 | l=2 (L_corr=2) | l=2 (L_corr=3) -- it is growth in ANGULAR MOMENTUM / basis, NOT in
    the number of coupled legs.  And the non-collapse is already present at the SMALLEST
    basis: 87.5% of external pairs have shared-vertex rank >= 2 at l=1.  The TC obstruction
    fires at TWO correlator legs (Y.Y is already a multiplet).  There is no 2->3 threshold
    on that axis at all.

(b) A THEOREM SEPARATES THEM.  The TC obstruction is compact-group representation theory:
    SO(3)/SU(2) coupling.  By Peter-Weyl every compact group's representation algebra is
    TYPE I -- tame -- however non-abelian it is.  So the TC obstruction can never BE wildness.
    Three orthogonal projections are not a compact-group representation (arbitrary subspaces
    escape Peter-Weyl), which is exactly how they reach the wild regime.

THIS DRIVER measures (b): it separates "non-abelian" from "irreducible/wild" by computing
commutant dimensions.  A non-abelian compact-group algebra stays REDUCIBLE (commutant dim =
sum of squared multiplicities > 1); three molecular center-projections are IRREDUCIBLE
(commutant dim 1).  Non-abelian therefore does not imply wild, and the two axes are
separated by a theorem rather than by a measurement.
"""
from __future__ import annotations
import importlib.util
import json
import sys

import numpy as np

sys.path.insert(0, "debug")
TOL = 1e-8


def _load(name, path):
    s = importlib.util.spec_from_file_location(name, path)
    m = importlib.util.module_from_spec(s)
    s.loader.exec_module(m)
    return m


LIN = _load("lin", "debug/beh2_ci_exact_landscape.py")


# ---------------------------------------------------------------- commutant
def commutant_dim(As):
    """dim {X : [A_i, X] = 0 for all i}."""
    n = As[0].shape[0]
    rows = [np.kron(A, np.eye(n)) - np.kron(np.eye(n), A.T) for A in As]
    s = np.linalg.svd(np.vstack(rows), compute_uv=False)
    return int((s < TOL * max(1.0, s[0])).sum()) + (n * n - len(s))


# ---------------------------------------------------------------- su(2) generators
def su2_generators(js):
    """Block-diagonal J_x, J_y, J_z for the direct sum of spin-j irreps in `js`."""
    blocks = {"x": [], "y": [], "z": []}
    for j in js:
        d = int(round(2 * j)) + 1
        m = np.array([j - k for k in range(d)])
        Jz = np.diag(m).astype(complex)
        off = np.sqrt(j * (j + 1) - m[1:] * (m[1:] + 1))
        Jp = np.zeros((d, d), complex)
        for k in range(d - 1):
            Jp[k, k + 1] = off[k]
        Jm = Jp.conj().T
        blocks["x"].append((Jp + Jm) / 2)
        blocks["y"].append((Jp - Jm) / (2j))
        blocks["z"].append(Jz)
    from scipy.linalg import block_diag
    return [block_diag(*blocks[a]) for a in ("x", "y", "z")]


def nonabelian(As):
    """max ||[A_i, A_j]|| -- confirms the generators genuinely do not commute."""
    return max(float(np.linalg.norm(As[i] @ As[j] - As[j] @ As[i], 2))
               for i in range(len(As)) for j in range(i + 1, len(As)))


# ---------------------------------------------------------------- molecular projections
def beh2_projectors(d1, d2, n_centers=3):
    S1 = LIN._beh(d1)
    S2 = LIN.PAR @ LIN._beh(d2) @ LIN.PAR
    SH = LIN.PAR @ LIN._hh(d1 + d2) @ LIN.PAR
    I = np.eye(LIN.NS)
    G = np.block([[I, S1, S2], [S1.T, I, SH], [S2.T, SH.T, I]])
    X = np.linalg.cholesky(G).T
    Ps = [X[:, 2 * k:2 * k + 2] @ np.linalg.pinv(X[:, 2 * k:2 * k + 2]) for k in range(3)]
    return Ps[:n_centers]


# ---------------------------------------------------------------- run
if __name__ == "__main__":
    out = {}

    print("=== (b) Non-abelian does NOT imply wild: compact-group algebras stay reducible ===")
    print(f"{'object':44s} {'dim':>4} {'max||[A,B]||':>13} {'commutant':>10}  verdict")
    su2_cases = [
        ("SU(2) on l=0 + l=1 + l=2  (Gaunt-style menu)", [0, 1, 2]),
        ("SU(2) on l=1 + l=1        (multiplicity 2)", [1, 1]),
        ("SU(2) on l=1 + l=2 + l=2", [1, 2, 2]),
        ("SU(2) on l=1 alone        (irreducible)", [1]),
    ]
    for label, js in su2_cases:
        As = su2_generators(js)
        c = commutant_dim(As)
        na = nonabelian(As)
        n = As[0].shape[0]
        verdict = "IRREDUCIBLE" if c == 1 else f"REDUCIBLE ({c} blocks-worth)"
        print(f"{label:44s} {n:4d} {na:13.4f} {c:10d}  {verdict}")
        out[label] = dict(dim=n, noncommutativity=na, commutant=c)

    print()
    print("=== the contrast: molecular center-projections at the same non-commutativity ===")
    for n_c, tag in [(2, "TWO centers   (Halmos tame)"), (3, "THREE centers (*-wild point)")]:
        Ps = beh2_projectors(2.6, 2.6, n_c)
        c = commutant_dim(Ps)
        na = nonabelian(Ps) if n_c > 1 else 0.0
        verdict = "IRREDUCIBLE" if c == 1 else f"REDUCIBLE ({c})"
        print(f"{tag:44s} {Ps[0].shape[0]:4d} {na:13.4f} {c:10d}  {verdict}")
        out[tag] = dict(dim=int(Ps[0].shape[0]), noncommutativity=na, commutant=c)

    print()
    print("READING: every SU(2) case is strongly non-abelian yet REDUCIBLE whenever the rep")
    print("carries more than one irrep -- Peter-Weyl type I (tame).  The three-center")
    print("projection algebra, at comparable non-commutativity, is IRREDUCIBLE.  Non-abelian")
    print("and wild are independent properties; the TC and composition obstructions live on")
    print("different axes and are separated by a theorem, not by a measurement.")

    with open("debug/data/aha_axis_check_tame_vs_wild.json", "w") as fh:
        json.dump(out, fh, indent=2)
    print("\nwrote debug/data/aha_axis_check_tame_vs_wild.json")
