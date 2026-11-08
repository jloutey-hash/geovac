"""I/O ladder Rung 3 — Load-vs-Generate block-encoding cost model (MODEL beat).

Diagnostic driver. Evaluates the three resource ledgers (T-count, ancilla, depth)
for the LOAD (QROM/QROAM) vs GENERATE (arithmetic + seed-load) coefficient sub-block
inside a qubitized PREPARE, and reports the Load-Generate (LG) crossover.

It invents NO physics numbers: L, G, S, b, lambda are FREE PARAMETERS supplied by
later rungs. This file only exercises the closed-form cost expressions from the memo
(debug/sprint_io_ladder_costmodel_memo.md) so the inequalities are concrete.

Cost expressions (per PREPARE query; SELECT ~ Theta(N) and reflection are common to
both schemes and drop out of the comparison):

  LOAD  (ancilla-min QROM,  Babbush et al. PRX 8, 041015 (2018)): T ~ 4L,        A ~ log2 L
  LOAD  (QROAM optimum,     Low-Kliuchnikov-Schaeffer Q 8, 1375): T ~ sqrt(L*b), A ~ sqrt(L/b)
  GEN   (arith + seed QROM):                                      T ~ G + 4S,    A ~ w_arith + log2 S
  GEN   (arith + seed QROAM):                                     T ~ G + sqrt(S*b)

Total simulation cost ~ O(lambda/eps) * (per-query T). lambda cancels in the per-query
comparison ONLY if lambda_gen == lambda_load; otherwise use the full product test (LG-3).
"""
from __future__ import annotations
import math
from dataclasses import dataclass


@dataclass
class Ledger:
    t_count: float
    ancilla: float
    label: str


def load_qrom(L: float, b: float) -> Ledger:
    """Ancilla-minimal QROM: T = 4L - 4 (word-length independent), A ~ log2 L."""
    return Ledger(t_count=4.0 * L - 4.0, ancilla=math.log2(max(L, 2.0)), label="LOAD/QROM")


def load_qroam(L: float, b: float) -> Ledger:
    """Space-time QROAM optimum: T ~ 2*sqrt(L*b), A ~ sqrt(L/b) dirty qubits."""
    return Ledger(t_count=2.0 * math.sqrt(L * b), ancilla=math.sqrt(L / b), label="LOAD/QROAM")


def gen_seedqrom(G: float, S: float, b: float, w_arith: float) -> Ledger:
    """Generate: fixed arithmetic G + ancilla-minimal seed load 4S."""
    return Ledger(t_count=G + 4.0 * S, ancilla=w_arith + math.log2(max(S, 2.0)), label="GEN/seedQROM")


def gen_seedqroam(G: float, S: float, b: float, w_arith: float) -> Ledger:
    """Generate: fixed arithmetic G + QROAM seed load ~2*sqrt(S*b)."""
    return Ledger(t_count=G + 2.0 * math.sqrt(S * b), ancilla=w_arith + math.sqrt(S / b), label="GEN/seedQROAM")


def crossover_report(L: float, G: float, S: float, b: float = 20.0,
                     w_arith: float = 40.0,
                     lam_load: float = 1.0, lam_gen: float = 1.0) -> None:
    print(f"\n=== (L={L:g}, G={G:g}, S={S:g}, b={b:g}, w_arith={w_arith:g}, "
          f"lambda_load={lam_load:g}, lambda_gen={lam_gen:g}) ===")

    lo_q = load_qrom(L, b)
    lo_a = load_qroam(L, b)
    ge_q = gen_seedqrom(G, S, b, w_arith)
    ge_a = gen_seedqroam(G, S, b, w_arith)

    for led in (lo_q, lo_a, ge_q, ge_a):
        print(f"  {led.label:14s}  T={led.t_count:12.1f}   ancilla={led.ancilla:8.1f}")

    # LG-1 (linear-QROM regime): G + S < L
    lg1 = (G + S) < L
    # LG-2 (QROAM regime): G + sqrt(S*b) < sqrt(L*b)
    lg2 = (G + math.sqrt(S * b)) < math.sqrt(L * b)
    # LG-3 (full test with per-query cost and 1-norm): lambda_gen*C_gen < lambda_load*C_load
    #   using the QROAM ledgers as the space-time-optimal representative
    lhs3 = lam_gen * ge_a.t_count
    rhs3 = lam_load * lo_a.t_count
    lg3 = lhs3 < rhs3

    print(f"  LG-1  G+S < L            : {G + S:g} < {L:g}   -> {'GEN wins' if lg1 else 'LOAD wins'}")
    print(f"  LG-2  G+sqrt(Sb)<sqrt(Lb): {G + math.sqrt(S*b):.1f} < {math.sqrt(L*b):.1f}"
          f"   -> {'GEN wins' if lg2 else 'LOAD wins'}")
    print(f"  LG-3  lam*C (QROAM)      : {lhs3:.1f} < {rhs3:.1f}"
          f"   -> {'GEN wins' if lg3 else 'LOAD wins'}   (load-bearing: uses BOTH 1-norms)")


if __name__ == "__main__":
    print(__doc__)
    print("Illustrative regimes ONLY (parameters are placeholders, not GeoVac numbers):")

    # A: large tensor, cheap local evaluator, few seeds, equal 1-norm -> generate wins cleanly
    crossover_report(L=1e6, G=500.0, S=200.0, b=20.0, lam_load=1.0, lam_gen=1.0)

    # B: same, but generate forces a 30x larger 1-norm (plane-wave-style penalty)
    #    -> LG-1/LG-2 still say "generate", but LG-3 can flip. Shows why I/O-only misleads.
    crossover_report(L=1e6, G=500.0, S=200.0, b=20.0, lam_load=1.0, lam_gen=30.0)

    # C: small tensor / expensive evaluator -> load wins (no asymptotic regime yet)
    crossover_report(L=800.0, G=5000.0, S=400.0, b=20.0, lam_load=1.0, lam_gen=1.0)

    print("\nNote: swap in real (L,G,S,b,lambda) from later rungs; this driver invents none.")
