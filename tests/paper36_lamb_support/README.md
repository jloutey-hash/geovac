# Paper 36 Lamb-chain support (durable test backing)

Backing module for `tests/test_paper36_lamb_chain.py`, which pins the
Paper 36 (`papers/group5_qed_gauge/paper_36_bound_state_qed.tex`)
LS-1 -> LS-6a Lamb-shift chain (predicted 1052.19 MHz) and the
GeoVac-native Sturmian Bethe logarithms.

## Why this directory exists

The original computations live in **archived sprint drivers** under
`debug/archive/`:

| Sprint | Archived driver | What it computed |
|:-------|:----------------|:-----------------|
| LS-3 | `debug/archive/qed_arc/ls3_bethe_log_regularized.py` | acceleration-form Bethe logs (2S = 2.786 at N=40) |
| LS-4 | `debug/archive/precision_arc/ls4_bethe_log_drake.py` | Drake--Swainson regularized l>0 logs (2P = -0.0307, 3D = -0.005236 at N=40) |
| LS-6a | `debug/archive/precision_arc/ls6a_eides_convention.py` | Eides 3.2 convention assembly (Lamb = 1052.19 MHz) |

`debug/` is the transient clean-room directory and is pruned by design
(CLAUDE.md section 9), so tests must not import from it (the same
silent-loss hazard fixed for the WH7 cluster in v4.49.1 by
`tests/wh7_support/`).  `bethe_log_sturmian.py` is a minimal, faithful
port of the shared Sturmian machinery plus the two Bethe-log forms;
the LS-6a closed-form assembly is short enough that it is recomputed
inline in the test itself.

## Provenance and fidelity

The port preserves the archived algorithms verbatim (same basis
chi_{k,l}(r) = r^{l+1} e^{-lam r} L_k^{(2l+1)}(2 lam r), same mpmath
exact polynomial-times-exponential integrals at dps=50, same
loewdin-style canonical orthogonalization, same Drake structural
denominator D_drake = 2(2l+1) Z^4 / n^3).  The only omissions are the
LS-4 K-sweep diagnostics (the l>0 spectral sum is K-independent by
I_v = 0, so a single assembly suffices) and plotting/JSON I/O.

Validated 2026-07-02 against fresh runs of the archived drivers:
identical values at every pinned (n, l, N).
