"""Door 3 — full-coverage test of the Paper 35 pi-criterion.

The classification itself is an enumeration over Paper 34's 28 projections + the
master Mellin engine M1/M2/M3 + the Paper 32 VIII case-exhaustion list; that is
documented in debug/door3_pi_criterion_coverage_memo.md (no computation needed for
the table — every transcendental signature is already pinned in the papers).

The ONE computation worth recording is the adversarial-case-2 spot-check: the
apparatus-identity / state-side projection (Paper 34 III.28) produces transcendentals
PSLQ-disjoint from the master Mellin engine M1 u M2 u M3 (Sprint TD Track 5). The
pi-criterion is about pi specifically, so the question is: do those von Neumann
correlation entropies contain pi? Sprint TD Track 5 says NO (PSLQ-null vs pi at
100 dps, 12,312-form basis). This spot-check confirms the documented values are not
simple pi-multiples, consistent with that verdict -> apparatus identity is CONSISTENT
with the criterion (transcendental, but pi-free).
"""

import math

# Documented framework correlation entropies (Sprint TD Track 5, Paper 34 III.28).
ENTROPIES = {
    "He (n_max=2)": 0.040811051,
    "Li+ (n_max=2)": 0.011211717,
    "He (n_max=4)": 0.041879430,
}


def pi_probe(S: float) -> dict:
    """Return S divided/multiplied by small pi-powers; none should be a simple rational."""
    return {
        "S/pi": S / math.pi,
        "S*pi": S * math.pi,
        "S/pi^2": S / math.pi**2,
        "S*pi^2": S * math.pi**2,
        "S/(4pi)": S / (4 * math.pi),
        "exp(S)": math.exp(S),
    }


def near_rational(x: float, max_den: int = 64, tol: float = 1e-6) -> bool:
    """Crude rational detector: is x within tol of p/q for q <= max_den?"""
    for q in range(1, max_den + 1):
        p = round(x * q)
        if p != 0 and abs(x - p / q) < tol:
            return True
    return False


def main() -> None:
    print("Door 3 adversarial case 2: apparatus-identity von Neumann entropy pi-probe")
    print("=" * 72)
    any_pi_hit = False
    for name, S in ENTROPIES.items():
        print(f"\n{name}: S = {S}")
        for lbl, v in pi_probe(S).items():
            hit = near_rational(v)
            flag = "  <-- RATIONAL?" if hit and lbl != "exp(S)" else ""
            if hit and lbl != "exp(S)":
                any_pi_hit = True
            print(f"   {lbl:10s} = {v:.10f}{flag}")
    print("\n" + "=" * 72)
    if any_pi_hit:
        print("WARNING: a pi-probe landed on a rational -> investigate (would threaten verdict).")
    else:
        print("VERDICT: no pi-probe lands on a simple rational across all three systems.")
        print("Consistent with Sprint TD Track 5 PSLQ-null-vs-pi: von Neumann entropy is")
        print("pi-FREE. Apparatus-identity projection is CONSISTENT with the criterion")
        print("(transcendental log-class, but contains no pi -> no integration, no pi).")


if __name__ == "__main__":
    main()
