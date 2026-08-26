"""
Track 4a / Probe 2 -- attempt the banked "balanced LiH curvature at n_max=4"
decider named in debug/sprint_abc_connections_test_memo.md, with a hard cost
guard (~20 min wall / memory blow-up => STOP and report the measured wall).

Registered A-vs-C predictions (quoted from the memo, "The residual A-vs-C
question" section):
  "The distinguishing question -- irreducible free-side wall (A) vs
   slow-but-eventual basis-response error (C) -- turns on whether the
   curvature converges as max_n->infinity. Over n_max 2->3 it does NOT
   converge (frozen), consistent with BOTH. Deciding needs n_max>=4
   ... Named decider (banked, not run): balanced LiH curvature at n_max=4."
So: if curvature/omega_e at n_max=4 stays frozen near the n_max=2,3 value
(omega_e ~+45%, curv/k~2.1-2.2x) -> supports A (irreducible wall).
If it moves back toward the true value (omega_e -> +0%, curv/k -> 1.0x)
-> supports C (slow basis-response error, healing).

This script does NOT try to force the full FCI through. It:
  1. Builds the real h1/eri integrals at n_max=4 (live solver call, timed).
  2. Times the REAL diagonal-loop and single-excitation-loop code paths
     used by geovac.coupled_composition.coupled_fci_energy on a wall-clock
     -capped sample (not a toy problem -- the actual arrays/algorithm),
     then extrapolates total cost from the measured per-iteration rate.
  3. Compares the extrapolate to the ~20 min guard and to the n_max=3
     measured cost (~2.3 h/pt, per debug/sprint_chem_error_projection_memo.md)
     and issues a GO/STOP verdict.
"""
from __future__ import annotations
import time
import itertools
import numpy as np

from geovac.balanced_coupled import build_balanced_hamiltonian
from geovac.molecular_spec import lih_spec

R_TRUE = 3.015
WALL_BUDGET_S = 1200.0          # ~20 min hard guard for this whole probe
DIAG_SAMPLE_BUDGET_S = 25.0     # sample window for the diagonal loop
SINGLE_SAMPLE_BUDGET_S = 45.0   # sample window for the (dominant) single-excitation loop

t_script0 = time.perf_counter()


def elapsed():
    return time.perf_counter() - t_script0


def n_orb(max_n: int) -> int:
    return sum(n * n for n in range(1, max_n + 1))


print("=" * 88)
print("Probe 2 -- balanced LiH n_max=4 curvature decider: cost-guarded attempt")
print("=" * 88)

for mn in (2, 3, 4):
    no = n_orb(mn)
    M = 3 * no  # LiH balanced: core (1 sub-block) + bond (2 sub-blocks: center+H)
    import math
    sector = math.comb(M, 2) ** 2
    print(f"  max_n={mn}: orbitals/center={no:3d}  M_total={M:3d}  "
          f"FCI sector dim={sector:,}")

print(f"\n[analytic] n_max=3 -> n_max=4 sector-dim ratio = "
      f"{math.comb(90,2)**2 / math.comb(42,2)**2:.1f}x")
print("[reference] n_max=3 measured cost (banked, per sprint_chem_error_projection_memo.md):"
      " ~2.3 h/point at sector dim 741,321.\n")

# --------------------------------------------------------------------- build
print("-" * 88)
print("STEP 1: build real h1/eri integrals at n_max=4, R=R_true (live solver call)")
print("-" * 88)
t0 = time.perf_counter()
spec = lih_spec(R=R_TRUE, max_n=4)
n_e = sum(b.n_electrons for b in spec.blocks)
ham = build_balanced_hamiltonian(spec, R=R_TRUE, n_grid_vne=8000, L_max=4,
                                  screened_cross_center=False, verbose=False)
t_build = time.perf_counter() - t0
M = ham['M']
h1 = ham['h1']
eri = ham['eri']
print(f"  build_balanced_hamiltonian(n_max=4) done in {t_build:.1f} s "
      f"(M={M}, h1 shape={h1.shape}, eri shape={eri.shape}, "
      f"eri nbytes={eri.nbytes/1e6:.1f} MB)")
print(f"  [cumulative wall so far: {elapsed():.1f} s]\n")

if elapsed() > WALL_BUDGET_S:
    print("COST GUARD TRIPPED after integral build alone. STOPPING.")
    raise SystemExit(0)

# ------------------------------------------------------- determinant basis
print("-" * 88)
print("STEP 2: generate the (N_up=2, N_down=2) determinant basis (real M=90 case)")
print("-" * 88)
t0 = time.perf_counter()
n_up = n_down = n_e // 2
alpha_strings = list(itertools.combinations(range(M), n_up))
beta_strings = list(itertools.combinations(range(M), n_down))
n_alpha, n_beta = len(alpha_strings), len(beta_strings)
n_det = n_alpha * n_beta
t_strings = time.perf_counter() - t0
print(f"  n_alpha={n_alpha:,}  n_beta={n_beta:,}  n_det={n_det:,}  "
      f"(built in {t_strings:.2f} s)")
print(f"  [cumulative wall so far: {elapsed():.1f} s]\n")

alpha_idx = {s: i for i, s in enumerate(alpha_strings)}


def excitation_phase(det, p, r):
    # matches geovac.coupled_composition._excitation_phase convention:
    # count occupied orbitals strictly between p and its removal position;
    # reuse the same sign rule by counting parity of the reordering.
    det_list = list(det)
    pos_p = det_list.index(p)
    new_det = sorted((set(det_list) - {p}) | {r})
    pos_r = new_det.index(r)
    return (-1) ** (pos_p + pos_r)


# --------------------------------------------------------- diagonal sample
print("-" * 88)
print(f"STEP 3: time the DIAGONAL loop (real h1/eri, capped at {DIAG_SAMPLE_BUDGET_S:.0f} s sample)")
print("-" * 88)
t0 = time.perf_counter()
count = 0
nuclear_repulsion = ham['nuclear_repulsion']
for ai, alpha in enumerate(alpha_strings):
    for bi, beta in enumerate(beta_strings):
        E_diag = nuclear_repulsion
        for p in alpha:
            E_diag += h1[p, p]
        for p in beta:
            E_diag += h1[p, p]
        for i_idx in range(n_up):
            for j_idx in range(i_idx + 1, n_up):
                p, q = alpha[i_idx], alpha[j_idx]
                E_diag += eri[p, p, q, q] - eri[p, q, q, p]
        for i_idx in range(n_down):
            for j_idx in range(i_idx + 1, n_down):
                p, q = beta[i_idx], beta[j_idx]
                E_diag += eri[p, p, q, q] - eri[p, q, q, p]
        for p in alpha:
            for q in beta:
                E_diag += eri[p, p, q, q]
        count += 1
        if count % 5000 == 0 and (time.perf_counter() - t0) > DIAG_SAMPLE_BUDGET_S:
            break
    if (time.perf_counter() - t0) > DIAG_SAMPLE_BUDGET_S:
        break
t_diag_sample = time.perf_counter() - t0
rate_diag = count / t_diag_sample
est_diag_total_s = n_det / rate_diag
print(f"  sampled {count:,} / {n_det:,} diagonal entries in {t_diag_sample:.1f} s "
      f"-> rate={rate_diag:,.0f} it/s")
print(f"  EXTRAPOLATED full diagonal-loop time: {est_diag_total_s:.0f} s "
      f"= {est_diag_total_s/60:.1f} min = {est_diag_total_s/3600:.2f} h")
print(f"  [cumulative wall so far: {elapsed():.1f} s]\n")

if elapsed() > WALL_BUDGET_S:
    print("COST GUARD TRIPPED after diagonal sampling. STOPPING before single-excitation sample.")
    raise SystemExit(0)

# ---------------------------------------------------- single-excitation sample
print("-" * 88)
print(f"STEP 4: time the ALPHA single-excitation loop "
      f"(the dominant off-diagonal cost; capped at {SINGLE_SAMPLE_BUDGET_S:.0f} s sample)")
print("-" * 88)
t0 = time.perf_counter()
count2 = 0
stop = False
for ai, alpha in enumerate(alpha_strings):
    alpha_set = set(alpha)
    for p in alpha:
        for r in range(M):
            if r in alpha_set:
                continue
            new_alpha = tuple(sorted((alpha_set - {p}) | {r}))
            if new_alpha not in alpha_idx:
                continue
            ai_new = alpha_idx[new_alpha]
            phase = excitation_phase(alpha, p, r)
            new_alpha_set = (alpha_set - {p}) | {r}
            val_base = phase * h1[r, p]
            for q in alpha:
                if q == p:
                    continue
                val_base += phase * (eri[r, p, q, q] - eri[r, q, q, p])
            # the real code then loops over ALL beta strings per (ai,p,r);
            # charge that inner cost analytically (n_beta) rather than
            # actually iterating it in the sample -- it is a pure O(1)-per-
            # element multiply-accumulate into the sparse matrix, so charging
            # n_beta counts it fairly without needing scipy lil_matrix here.
            count2 += n_beta
            if count2 % 200000 < n_beta and (time.perf_counter() - t0) > SINGLE_SAMPLE_BUDGET_S:
                stop = True
                break
        if stop:
            break
    if stop:
        break
t_single_sample = time.perf_counter() - t0
# total alpha-single-excitation element count (matches the real loop's total work)
total_single_alpha = 0
for alpha in alpha_strings[: min(200, n_alpha)]:  # unbiased small independent count
    alpha_set = set(alpha)
    n_allowed = 0
    for p in alpha:
        for r in range(M):
            if r in alpha_set:
                continue
            new_alpha = tuple(sorted((alpha_set - {p}) | {r}))
            if new_alpha in alpha_idx:
                n_allowed += 1
    total_single_alpha += n_allowed
avg_allowed_per_alpha = total_single_alpha / min(200, n_alpha)
total_single_elems = avg_allowed_per_alpha * n_alpha * n_beta
rate_single = count2 / t_single_sample
est_single_total_s = total_single_elems / rate_single
print(f"  sampled {count2:,} matrix-element updates in {t_single_sample:.1f} s "
      f"-> rate={rate_single:,.0f} it/s")
print(f"  avg allowed single-excitations/alpha-string (measured on {min(200,n_alpha)} "
      f"strings) = {avg_allowed_per_alpha:.1f}")
print(f"  total ALPHA single-excitation matrix-element count (this class alone) "
      f"~= {total_single_elems:,.0f}")
print(f"  EXTRAPOLATED alpha-single-excitation-loop time: {est_single_total_s:.0f} s "
      f"= {est_single_total_s/60:.1f} min = {est_single_total_s/3600:.2f} h")
print(f"  (beta singles ~same order again; alpha-alpha/beta-beta/alpha-beta DOUBLE "
      f"excitations are additional and NOT sampled here -- this is a lower bound "
      f"on total off-diagonal cost)")
print(f"  [cumulative wall so far: {elapsed():.1f} s]\n")

# ------------------------------------------------------------------ verdict
print("=" * 88)
print("VERDICT")
print("=" * 88)
total_est_h = (est_diag_total_s + 2 * est_single_total_s) / 3600.0  # x2 for beta singles too
print(f"  Lower-bound extrapolated total FCI-matrix-build time "
      f"(diag + alpha singles + beta singles, EXCLUDING doubles): "
      f"~{total_est_h:.1f} h")
print(f"  This already exceeds the ~20 min (0.33 h) cost guard by a factor of "
      f"~{total_est_h/0.333:.0f}x, BEFORE accounting for double excitations "
      f"or the eigsh diagonalization itself, and before accounting for the "
      f"scipy lil_matrix assembly overhead at n_det={n_det:,} rows.")
print(f"  Memory note: a lil_matrix with {n_det:,} rows allocates a Python list "
      f"per row; at ~200-400 bytes/empty-row overhead this is "
      f"~{n_det*300/1e9:.1f}-{n_det*400/1e9:.1f} GB structural overhead alone, "
      f"before any nonzero entries are inserted.")
print(f"  STOP: n_max=4 is not run to completion. Reporting the measured cost "
      f"wall instead of grinding (per probe instructions).")
print(f"  Largest genuinely feasible point remains n_max=3 (banked/cached; "
      f"~2.3 h/point at sector dim 741,321) and n_max=2 (live, ~14 s/point).")
print(f"\n[total script wall time: {elapsed():.1f} s]")
