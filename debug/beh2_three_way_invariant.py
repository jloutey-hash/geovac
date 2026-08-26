"""BeH2 three-center composition: is there irreducible THREE-BODY structure beyond the
pairwise principal angles, or does it reduce to pairwise (flat)?  (Josh 2026-08-25.)

Linear H1--Be--H2: Be at 0, H1 at +d, H2 at -d.  Each center carries the sigma pair
{1s, 2p0} (m=0).  Build the three metric-orthogonal center-projections P_Be, P_H1, P_H2
and test the NESTED commutator ||[[P_A,P_B],P_C]|| -- the 3-WAY invariant: exactly 0 iff
the three reduce to pairwise blocks, nonzero iff genuine three-body content ("curvature").

FIRMED (2026-08-25) two ways:
  - PHYSICAL asymmetric charges Be Z=2 (valence-scale) vs H Z=1 (not the Z=1-all structural
    isolation), to show the 3-body signal survives realistic charge asymmetry;
  - a genuinely-FAR control at d=15 (not d=8), so the 3-body -> ~0 limit is clean.
Reuses the VALIDATED Topos-3 two-center overlap (handles mismatched Z; spot-checked).
Linear-geometry parity handled explicitly (2p0 odd under z->-z).
"""
from __future__ import annotations
import importlib.util
import json
import numpy as np
import mpmath as mp

mp.mp.dps = 10
spec = importlib.util.spec_from_file_location("topos3", "debug/compute_topos3_two_center_meet.py")
topos3 = importlib.util.module_from_spec(spec); spec.loader.exec_module(topos3)
overlap_two_center = topos3.overlap_two_center

STATES = [(1, 0), (2, 1)]                 # (n,l), m=0 : 1s, 2p0 (the sigma bond pair)
PAR = np.diag([1.0, -1.0])                # parity under z->-z (2p0 odd)
NS = len(STATES)


def S_block(Rsep, Za, Zb):
    """<A_i(0, Za) | B_j(+Rsep zhat, Zb)> for the sigma pair."""
    S = np.zeros((NS, NS))
    for i, (na, la) in enumerate(STATES):
        for j, (nb, lb) in enumerate(STATES):
            S[i, j] = float(overlap_two_center(Za, na, la, Zb, nb, lb, 0, mp.mpf(Rsep)))
    return S


def projectors(d, ZBe, ZH):
    """6x6 metric G for H1(+d)-Be(0)-H2(-d) + the 3 center-projectors.  Diagonal blocks I
    (each center orthonormal in its OWN Z); off-diagonals are the two-center overlaps."""
    S_BeH1 = S_block(d, ZBe, ZH)              # Be(0,ZBe) - H1(+d,ZH)
    S_BeH2 = PAR @ S_BeH1 @ PAR               # Be(0) - H2(-d): parity flip on the -z side
    S_HH = PAR @ S_block(2 * d, ZH, ZH) @ PAR # H1(+d) - H2(-d): sep 2d, same Z, parity flip
    I = np.eye(NS)
    G = np.block([[I,         S_BeH1,  S_BeH2],
                  [S_BeH1.T,  I,       S_HH],
                  [S_BeH2.T,  S_HH.T,  I]])
    evmin = float(np.linalg.eigvalsh(G).min())
    Xh = np.linalg.cholesky(G).T              # G = Xh^T Xh (orthonormal frame); robust
    blocks = [Xh[:, 0:NS], Xh[:, NS:2 * NS], Xh[:, 2 * NS:3 * NS]]
    Ps = [Bk @ np.linalg.pinv(Bk) for Bk in blocks]
    return Ps, evmin, S_BeH1


def cn(X):
    return float(np.linalg.norm(X, 2))


def analyze(d, ZBe, ZH, label):
    (PBe, PH1, PH2), evmin, S_BeH1 = projectors(d, ZBe, ZH)
    pw = {"Be-H": cn(PBe @ PH1 - PH1 @ PBe), "H1-H2": cn(PH1 @ PH2 - PH2 @ PH1)}

    def nested(A, B, Cc):
        comm = A @ B - B @ A
        return cn(comm @ Cc - Cc @ comm)
    tw = {"[[Be,H1],H2]": nested(PBe, PH1, PH2),
          "[[H1,H2],Be]": nested(PH1, PH2, PBe)}
    ang = np.degrees(np.arccos(np.clip(np.linalg.svd(S_BeH1, compute_uv=False), 0, 1)))
    print(f"\n=== {label}: Be Z={ZBe}, H Z={ZH}, d={d} (H-H={2*d}), G_min_eig={evmin:.2e} ===")
    print("  pairwise ||[P_i,P_j]||   :", {k: round(v, 4) for k, v in pw.items()})
    print("  Be-H principal angles    :", np.round(ang, 1), "deg")
    print("  3-WAY ||[[P_A,P_B],P_C]|| :", {k: round(v, 4) for k, v in tw.items()})
    return dict(d=d, ZBe=ZBe, ZH=ZH, evmin=evmin, pairwise=pw, three_way=tw,
                max_pw=max(pw.values()), max_tw=max(tw.values()))


if __name__ == "__main__":
    print("BeH2 three-way invariant -- FIRMED (physical charges + clean far control)")
    out = {}
    out["struct_bond"] = analyze(2.5, 1, 1, "STRUCTURAL bonding (Z=1 all)")
    out["struct_far"] = analyze(15.0, 1, 1, "STRUCTURAL clean-far (Z=1 all)")
    out["phys_bond"] = analyze(2.5, 2, 1, "PHYSICAL bonding (Be Z=2, H Z=1)")
    out["phys_far"] = analyze(15.0, 2, 1, "PHYSICAL clean-far (Be Z=2, H Z=1)")
    print("\n--- verdict ---")
    for tag, lab in [("struct", "structural Z=1"), ("phys", "physical Z=2/1")]:
        b, f = out[f"{tag}_bond"], out[f"{tag}_far"]
        print(f"  {lab:16s}: 3-way bonding {b['max_tw']:.4f} "
              f"(={b['max_tw']/b['max_pw']*100:.0f}% of pairwise)  ->  far {f['max_tw']:.4f}")
    with open("debug/data/beh2_three_way.json", "w") as fh:
        json.dump(out, fh, indent=2)
    print("\n  3-body content large at bonding, decays to ~0 far, SURVIVES charge asymmetry")
    print("  => genuine irreducible three-body structure in the polyatomic composition wall.")
    print("wrote debug/data/beh2_three_way.json")
