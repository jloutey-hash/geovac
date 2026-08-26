"""Validation ladder for the N-electron correlated-CI engine (step 5, revised).

GATE 1  gamma = 0 EXACTLY.  Then f == 1, f' == 0, so F == n_pairs (a number) and
        |G> = n_pairs |Phi_0>.  Therefore, ELEMENT BY ELEMENT,
            S_FG[I] == n_pairs * S[I, Phi_0],    H_FG[I] == n_pairs * H[I, Phi_0]
            S_GG    == n_pairs^2 * S[Phi_0,Phi_0], H_GG == n_pairs^2 * H[Phi_0,Phi_0]
        This checks the F-carrying assembly against the plain assembly at the level
        of individual matrix elements, not just an energy.

        (The earlier "gamma -> 0" energy control was WRONG and is retired: as gamma->0
        the geminal's residual vanishes in NORM, but the surviving direction is
        proportional to sum_ij r_ij -- a genuine Hylleraas-type correlating function
        whose energy contribution does NOT vanish.  That is why it "failed" at 1e-4.)

GATE 2  N = 2, no geminal, vs debug/ctf12_r12ci_he.py plain FCI.
GATE 2b N = 3, no geminal, vs geovac.transcorrelated_sturmian plain FCI (independent
        N-electron determinant code).
GATE 3  N = 2, WITH geminal, vs the He engine -- using RAW Sturmians so that orbital 0
        is exactly R_1 and |G> = f * R_1(r1) R_1(r2) is the SAME FUNCTION the reference
        engine uses.  (With Loewdin orbitals the geminal sits on a different reference
        function, so the two calculations span different spaces and must not agree.)
"""
import importlib.util
import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
os.chdir(ROOT)
sys.path.insert(0, HERE)
sys.path.insert(0, ROOT)

import r12ci_ne_engine as E  # noqa: E402
from geovac import transcorrelated_sturmian as TC  # noqa: E402

_spec = importlib.util.spec_from_file_location(
    "ctf12_r12ci_he", os.path.join(HERE, "ctf12_r12ci_he.py"))
ctf = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(ctf)

NG, NX = 600, 128
out = {"gates": []}

print("=" * 80)
print("GATE 1 -- gamma = 0 exactly: F == n_pairs, so the F-blocks must reduce")
print("=" * 80)
print(f"{'N':>2}{'ns':>4}{'npairs':>8}{'max|S_FG dev|':>15}{'max|H_FG dev|':>15}"
      f"{'S_GG dev':>12}{'H_GG dev':>12}{'ok':>6}")
g1_ok = True
for (N, ns, ms2, Z, k, Ng) in [(2, 3, 0, 2.0, 1.7, 300), (3, 3, 1, 3.0, 2.4, 200)]:
    S1, H1, dets = E.build(ns, N, Z, k, 0.0, ms2, Ng=Ng, with_geminal=True,
                           Lmax=10, nx=160)
    nd = len(dets)
    npair = N * (N - 1) // 2
    i0 = 0                                   # ref_det defaults to dets[0]
    sdev = np.abs(S1[:nd, nd] - npair * S1[:nd, i0]).max()
    hdev = np.abs(H1[:nd, nd] - npair * H1[:nd, i0]).max()
    sgg = abs(S1[nd, nd] - npair ** 2 * S1[i0, i0])
    hgg = abs(H1[nd, nd] - npair ** 2 * H1[i0, i0])
    scale = max(abs(H1[i0, i0]), 1.0)
    ok = max(sdev, sgg) < 1e-9 and max(hdev, hgg) < 1e-8 * scale
    g1_ok &= ok
    print(f"{N:>2}{ns:>4}{npair:>8}{sdev:>15.2e}{hdev:>15.2e}{sgg:>12.2e}"
          f"{hgg:>12.2e}{'OK' if ok else 'FAIL':>6}")
    out["gates"].append(dict(gate=1, N=N, ns=ns, npairs=npair, s_fg_dev=sdev,
                             h_fg_dev=hdev, s_gg_dev=sgg, h_gg_dev=hgg, ok=bool(ok)))
print(f"GATE 1: {'PASS' if g1_ok else 'FAIL'}")

print()
print("=" * 80)
print("GATE 2 -- N=2, no geminal, vs the independent He engine (plain FCI)")
print("=" * 80)
print(f"{'ns':>3}{'this engine':>16}{'ctf12 engine':>16}{'diff':>12}{'ok':>6}")
g2_ok = True
for ns in (2, 3, 4):
    S0, H0, _ = E.build(ns, 2, 2.0, 1.7, 0.7, 0, Ng=NG, with_geminal=False, nx=NX)
    e0, _ = E.solve(S0, H0)
    eref = ctf.he_energy(ns, 1.7, 0.7, n_gem=0, Ng=NG, nx=NX)[0]
    ok = abs(e0 - eref) < 2e-6
    g2_ok &= ok
    print(f"{ns:>3}{e0:>16.9f}{eref:>16.9f}{e0 - eref:>12.2e}{'OK' if ok else 'FAIL':>6}")
    out["gates"].append(dict(gate=2, ns=ns, mine=e0, ref=eref, diff=e0 - eref, ok=bool(ok)))
print(f"GATE 2: {'PASS' if g2_ok else 'FAIL'}")

print()
print("=" * 80)
print("GATE 2b -- N=3 (Li), no geminal, vs geovac.transcorrelated_sturmian FCI")
print("=" * 80)
print(f"{'ns':>3}{'k':>7}{'this engine':>16}{'geovac module':>16}{'diff':>12}{'ok':>6}")
g2b_ok = True
for ns, k in [(2, 2.4), (3, 2.4), (3, 2.0)]:
    S0, H0, _ = E.build(ns, 3, 3.0, k, 0.7, 1, Ng=800, with_geminal=False, nx=NX)
    e0, _ = E.solve(S0, H0)
    syst = TC.build_atomic_system(ns, k, 0.7, 3, 3, Ng=800, nx=NX, with_L3=False)
    eref, _ = TC.ground(TC.build_fci_matrix(syst, syst.plain()), hermitian=True)
    ok = abs(e0 - eref) < 2e-6
    g2b_ok &= ok
    print(f"{ns:>3}{k:>7.2f}{e0:>16.9f}{eref:>16.9f}{e0 - eref:>12.2e}"
          f"{'OK' if ok else 'FAIL':>6}")
    out["gates"].append(dict(gate="2b", ns=ns, k=k, mine=e0, ref=eref,
                            diff=e0 - eref, ok=bool(ok)))
print(f"GATE 2b: {'PASS' if g2b_ok else 'FAIL'}")

print()
print("=" * 80)
print("GATE 3 -- N=2 WITH geminal vs the He engine, RAW Sturmians (matched geminal)")
print("=" * 80)
print(f"{'ns':>3}{'gamma':>8}{'this engine':>16}{'ctf12 engine':>16}{'diff':>12}{'ok':>6}")
g3_ok = True
for ns, gam in [(3, 0.4), (3, 0.7), (3, 1.1), (4, 0.4), (4, 0.7)]:
    S1, H1, dts = E.build(ns, 2, 2.0, 1.7, gam, 0, Ng=NG, with_geminal=True,
                          nx=NX, lowdin=False)
    e1, _ = E.solve(S1, H1, nd=len(dts))
    eref = ctf.he_energy(ns, 1.7, gam, n_gem=1, Ng=NG, nx=NX)[0]
    ok = abs(e1 - eref) < 5e-6
    g3_ok &= ok
    print(f"{ns:>3}{gam:>8.2f}{e1:>16.9f}{eref:>16.9f}{e1 - eref:>12.2e}"
          f"{'OK' if ok else 'FAIL':>6}")
    out["gates"].append(dict(gate=3, ns=ns, gamma=gam, mine=e1, ref=eref,
                            diff=e1 - eref, ok=bool(ok)))
print(f"GATE 3: {'PASS' if g3_ok else 'FAIL'}")

out["summary"] = dict(gate1=bool(g1_ok), gate2=bool(g2_ok), gate2b=bool(g2b_ok),
                      gate3=bool(g3_ok))
os.makedirs("debug/data", exist_ok=True)
with open("debug/data/r12ci_ne_gates.json", "w") as f:
    json.dump(out, f, indent=2)
print()
print(f"LADDER: g1={g1_ok}  g2={g2_ok}  g2b={g2b_ok}  g3={g3_ok}")
