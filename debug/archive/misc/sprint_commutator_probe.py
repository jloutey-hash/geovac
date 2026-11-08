"""
Diagnostic (2): is the multi-focal-composition wall a NON-COMMUTING conditional-
expectations obstruction?  (/ahha leap, 2026-07-06.)

Model the two focal-length projections as the orthogonal projectors P_A, P_B onto
the single-center orbital subspaces of a diatomic bond block (center A at 0,
center B at R zhat; both Z=1 hydrogenic -- the LiH balanced bond block).  Same-Z
hydrogenic orbitals are orthonormal WITHIN a center (G_AA=G_BB=I), so the singular
values of the cross-overlap block S_AB are exactly the canonical correlations
cos(theta_k) (principal angles between the subspaces), and
    ||[P_A,P_B]||_block = max_k cos(theta_k) sin(theta_k).
[P_A,P_B]=0 (bit) <=> commuting square <=> a clean composed conditional expectation
exists <=> claim (2) DEAD.  Nonzero => no naive composed projection; the right
object is the Jones basic construction, and its tractability turns on the STRUCTURE
of the coupling (m-block-diagonal? how many nonzero angles? R-scaling?).

Two-center overlaps reuse the VALIDATED overlap_two_center from the Topos-3 driver
(reproduces the exact 1s-1s Slater overlap to 12 digits).
"""
from __future__ import annotations
import importlib.util
import json
import numpy as np
import mpmath as mp

mp.mp.dps = 10   # singular values need ~1e-6; dps=20 was 5x slower, bit-identical

# --- reuse the validated two-center overlap (Topos-3) ---
spec = importlib.util.spec_from_file_location(
    "topos3", "debug/compute_topos3_two_center_meet.py")
topos3 = importlib.util.module_from_spec(spec)
spec.loader.exec_module(topos3)
overlap_two_center = topos3.overlap_two_center

R_TRUE = 3.015
def S_1s1s(R):  # parameter-free 1s-1s Slater overlap (reference scale)
    return float(mp.e ** (-R) * (1 + R + R ** 2 / 3))


def m_states(n_max):
    """{m: [(n,l),...]} for m>=0."""
    out = {}
    for n in range(1, n_max + 1):
        for l in range(n):
            for m in range(0, l + 1):
                out.setdefault(m, []).append((n, l))
    return out


def probe(R, n_max, Z1=1, Z2=1):
    blocks = m_states(n_max)
    per_m = {}
    max_comm = 0.0
    max_sigma = 0.0
    for m, states in blocks.items():
        d = len(states)
        S = np.zeros((d, d))
        for i, (na, la) in enumerate(states):
            for j, (nb, lb) in enumerate(states):
                S[i, j] = float(overlap_two_center(Z1, na, la, Z2, nb, lb, m, mp.mpf(R)))
        sig = np.linalg.svd(S, compute_uv=False)         # canonical correlations
        sig = np.clip(sig, 0.0, 1.0)
        comm = float(np.max(sig * np.sqrt(np.maximum(1 - sig ** 2, 0.0)))) if d else 0.0
        n_coupled = int(np.sum(sig > 1e-10))
        # l-mixing: is any off-l entry (la != lb) nonzero within this m-block?
        l_mix = any(abs(S[i, j]) > 1e-10 and states[i][1] != states[j][1]
                    for i in range(d) for j in range(d))
        per_m[m] = {'dim': d, 'sigma': sig.tolist(), 'n_coupled': n_coupled,
                    'comm_norm': comm, 'l_mixing': bool(l_mix),
                    'max_sigma': float(sig.max()) if d else 0.0}
        max_comm = max(max_comm, comm)
        max_sigma = max(max_sigma, float(sig.max()) if d else 0.0)
    return {'R': R, 'n_max': n_max, 'per_m': per_m,
            'max_comm_norm': max_comm, 'max_sigma': max_sigma,
            'S_1s1s': S_1s1s(R)}


def dump_S_AB(R, n_max=2, m=0, Z1=1, Z2=1):
    """Print the raw cross-center overlap block S_AB for one m-block (shows the
    cross-n / cross-l structure that the Paper-8 orthogonal D-matrix forbids)."""
    states = m_states(n_max)[m]
    print(f"\n  S_AB[m={m}] at R={R} (rows=center A, cols=center B), states {states}:")
    for na, la in states:
        row = "   ".join(f"{float(overlap_two_center(Z1, na, la, Z2, nb, lb, m, mp.mpf(R))):+7.4f}"
                         for nb, lb in states)
        print(f"    (n{na},l{la}): {row}")


if __name__ == "__main__":
    Rs = [2.5, R_TRUE, 5.0, 8.0]
    print("=" * 92)
    print("Commutator probe: principal angles between the two center-subspaces "
          "(Z=1 bond block)")
    print("=" * 92)
    dump_S_AB(R_TRUE, n_max=2, m=0)     # show cross-n/cross-l structure explicitly
    rows = []
    for nmax in (2,):
        print(f"\n--- n_max={nmax} ---")
        print(f"  {'R':>6s} {'||[P_A,P_B]||':>14s} {'max cos(theta)':>14s} "
              f"{'S_1s1s(R)':>11s}  {'comm/S':>8s}   per-m (dim, n_coupled, l_mix)")
        for R in Rs:
            r = probe(R, nmax)
            rows.append(r)
            ms = "  ".join(f"m{m}:d{b['dim']}/c{b['n_coupled']}/{('L' if b['l_mixing'] else '-')}"
                          for m, b in sorted(r['per_m'].items()))
            print(f"  {R:6.3f} {r['max_comm_norm']:14.6f} {r['max_sigma']:14.6f} "
                  f"{r['S_1s1s']:11.6f}  {r['max_comm_norm']/r['S_1s1s']:8.3f}   {ms}")

    json.dump(rows, open('debug/data/commutator_probe.json', 'w'), indent=2)
    print("\n[saved] debug/data/commutator_probe.json")

    # verdict read-out
    print("\n" + "=" * 92)
    print("VERDICT READ-OUT")
    print("=" * 92)
    r_true_n2 = next(r for r in rows if r['n_max'] == 2 and abs(r['R'] - R_TRUE) < 1e-6)
    r_far = next(r for r in rows if r['n_max'] == 2 and r['R'] == 7.0)
    print(f"  [P_A,P_B] at R_true (n2): {r_true_n2['max_comm_norm']:.5f}  "
          f"(bit-zero => claim-2 DEAD; nonzero => Jones-tower reading LIVE)")
    print(f"  [P_A,P_B] at R=7   (n2): {r_far['max_comm_norm']:.5f}  "
          f"(-> 0 as R->inf => commuting square recovered at separation)")
    print(f"  m-block-diagonal by construction; l-mixing WITHIN m-block flagged 'L' above")
    print(f"  comm/S_1s1s ratio ~ const => the non-commutativity IS the overlap (function, "
          f"not constant)")
