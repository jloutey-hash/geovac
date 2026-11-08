"""S^(4) k=4 o-space reduction identities — bit-exact verification.

Verifies identities (i1)-(i7) in Fraction arithmetic at N=8 (o-box [5,19])
and N=12 (o-box [5,27]), then assembles the full S^(4) o-space relation.

Pattern mirrors s3_decomp_setup.py V3.
"""
from __future__ import annotations

import json
import sys
import time
from fractions import Fraction
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

t0 = time.time()

# ---------------------------------------------------------------- helpers

def phi(n: int) -> Fraction:
    """phi(n) = 2/(n+3/2)^2 - (1/2)/(n+3/2)^4"""
    x = Fraction(2 * n + 3, 2)
    return 2 / x ** 2 - Fraction(1, 2) / x ** 4


def psi(o: int) -> Fraction:
    """psi(o) = 8(o^2-1)/o^4;  phi(n) = psi(2n+3)"""
    return Fraction(8 * (o * o - 1), o ** 4)


def build_nspace(N: int) -> dict:
    """Compute all n-space objects over rng=[1,N]."""
    rng = range(1, N + 1)
    PH = {n: phi(n) for n in rng}
    P  = sum(PH[n] for n in rng)
    Q  = sum(PH[n] ** 2 for n in rng)
    R3 = sum(PH[n] ** 3 for n in rng)
    R4 = sum(PH[n] ** 4 for n in rng)

    S_min = sum(min(a, b) * PH[a] * PH[b]
                for a in rng for b in rng)
    M2    = sum(min(a, b) * PH[a] * PH[b] ** 2
                for a in rng for b in rng)
    T31   = sum(min(a, b) * PH[a] * PH[b] ** 3
                for a in rng for b in rng)
    T22   = sum(min(a, b) * PH[a] ** 2 * PH[b] ** 2
                for a in rng for b in rng)

    M3 = sum(min(a, b) * min(b, c) * PH[a] * PH[b] * PH[c]
             for a in rng for b in rng for c in rng)
    M3e = sum(min(a, b) * min(b, c) * PH[a] * PH[b] * PH[c] ** 2
              for a in rng for b in rng for c in rng)
    M3m = sum(min(a, b) * min(b, c) * PH[a] * PH[b] ** 2 * PH[c]
              for a in rng for b in rng for c in rng)

    M4 = sum(min(a, b) * min(b, c) * min(c, d)
             * PH[a] * PH[b] * PH[c] * PH[d]
             for a in rng for b in rng for c in rng for d in rng)

    # S^(4) direct (15-term identity, verified bit-exact in scoping diag)
    S4 = (8 * M4
          - 8 * M3 * P
          - 8 * M3e
          - 4 * M3m
          + 8 * M2 * P
          + 6 * S_min * P ** 2
          + 4 * S_min * Q
          - 4 * S_min ** 2
          + 4 * T31
          + 2 * T22
          - 2 * R3 * P
          - 3 * P ** 2 * Q
          - P ** 4
          - Q ** 2
          - R4)

    return dict(P=P, Q=Q, R3=R3, R4=R4,
                S_min=S_min, M2=M2, T31=T31, T22=T22,
                M3=M3, M3e=M3e, M3m=M3m, M4=M4, S4=S4)


def build_ospace(N: int) -> dict:
    """Compute all o-space objects over odds=[5, 2N+3 step 2]."""
    O    = 2 * N + 3
    odds = range(5, O + 1, 2)
    PS   = {o: psi(o) for o in odds}

    S_min_o = sum(min(a, b) * PS[a] * PS[b]
                  for a in odds for b in odds)
    G       = sum(min(a, b) * PS[a] * PS[b] ** 2
                  for a in odds for b in odds)
    H31     = sum(min(a, b) * PS[a] * PS[b] ** 3
                  for a in odds for b in odds)
    H22     = sum(min(a, b) * PS[a] ** 2 * PS[b] ** 2
                  for a in odds for b in odds)
    C       = sum(min(a, b) * min(b, c) * PS[a] * PS[b] * PS[c]
                  for a in odds for b in odds for c in odds)
    Ce      = sum(min(a, b) * min(b, c) * PS[a] * PS[b] * PS[c] ** 2
                  for a in odds for b in odds for c in odds)
    Cm      = sum(min(a, b) * min(b, c) * PS[a] * PS[b] ** 2 * PS[c]
                  for a in odds for b in odds for c in odds)
    C4      = sum(min(a, b) * min(b, c) * min(c, d)
                  * PS[a] * PS[b] * PS[c] * PS[d]
                  for a in odds for b in odds
                  for c in odds for d in odds)

    return dict(S_min_o=S_min_o, G=G, H31=H31, H22=H22,
                C=C, Ce=Ce, Cm=Cm, C4=C4)


# ---------------------------------------------------------------- main loop

OUT = {}
all_pass = True

for N in (8, 12):
    print(f"\n=== N={N} (o-box [5,{2*N+3}]) ===")
    ns = build_nspace(N)
    os_= build_ospace(N)

    P      = ns["P"];      Q      = ns["Q"]
    R3     = ns["R3"];     R4     = ns["R4"]
    S_min  = ns["S_min"];  M2     = ns["M2"]
    T31    = ns["T31"];    T22    = ns["T22"]
    M3     = ns["M3"];     M3e    = ns["M3e"]
    M3m    = ns["M3m"];    M4     = ns["M4"]
    S4     = ns["S4"]

    S_min_o = os_["S_min_o"]
    G       = os_["G"];    H31    = os_["H31"]
    H22     = os_["H22"];  C      = os_["C"]
    Ce      = os_["Ce"];   Cm     = os_["Cm"]
    C4      = os_["C4"]

    # (i1) S_min_o = 2 S_min + 3 P^2
    lhs_i1 = S_min_o
    rhs_i1 = 2 * S_min + 3 * P ** 2
    v_i1 = lhs_i1 == rhs_i1
    print(f"  (i1) S_min_o = 2 S_min + 3 P^2 :  {'PASS' if v_i1 else 'FAIL'}")
    OUT[f"N{N}_i1"] = v_i1

    # (i2) M4 = (1/8)[C4 - 6 C P - 3 S_min_o^2 + 27 S_min_o P^2 - 27 P^4]
    rhs_i2 = Fraction(1, 8) * (C4
                                - 6 * C * P
                                - 3 * S_min_o ** 2
                                + 27 * S_min_o * P ** 2
                                - 27 * P ** 4)
    v_i2 = M4 == rhs_i2
    print(f"  (i2) M4 = (1/8)[...] :  {'PASS' if v_i2 else 'FAIL'}")
    OUT[f"N{N}_i2"] = v_i2

    # (i3) M3e = (1/4)[Ce - 3 S_min_o Q - 3 P G + 9 P^2 Q]
    rhs_i3 = Fraction(1, 4) * (Ce
                                - 3 * S_min_o * Q
                                - 3 * P * G
                                + 9 * P ** 2 * Q)
    v_i3 = M3e == rhs_i3
    print(f"  (i3) M3e = (1/4)[...] :  {'PASS' if v_i3 else 'FAIL'}")
    OUT[f"N{N}_i3"] = v_i3

    # (i4) M3m = (1/4)[Cm - 6 G P + 9 Q P^2]
    rhs_i4 = Fraction(1, 4) * (Cm - 6 * G * P + 9 * Q * P ** 2)
    v_i4 = M3m == rhs_i4
    print(f"  (i4) M3m = (1/4)[...] :  {'PASS' if v_i4 else 'FAIL'}")
    OUT[f"N{N}_i4"] = v_i4

    # (i5) T31 = (1/2)[H31 - 3 P R3]
    rhs_i5 = Fraction(1, 2) * (H31 - 3 * P * R3)
    v_i5 = T31 == rhs_i5
    print(f"  (i5) T31 = (1/2)[H31 - 3PR3] :  {'PASS' if v_i5 else 'FAIL'}")
    OUT[f"N{N}_i5"] = v_i5

    # (i6) T22 = (1/2)[H22 - 3 Q^2]
    rhs_i6 = Fraction(1, 2) * (H22 - 3 * Q ** 2)
    v_i6 = T22 == rhs_i6
    print(f"  (i6) T22 = (1/2)[H22 - 3Q^2] :  {'PASS' if v_i6 else 'FAIL'}")
    OUT[f"N{N}_i6"] = v_i6

    # (i7) k=3 controls (regression against s3_decomp_setup V3)
    v_i7a = (M3 == C / 4 - 3 * S_min * P - Fraction(9, 4) * P ** 3)
    v_i7b = (M2 == G / 2 - Fraction(3, 2) * P * Q)
    print(f"  (i7a) M3 = C/4 - 3 S_min P - (9/4)P^3 :  {'PASS' if v_i7a else 'FAIL'}")
    print(f"  (i7b) M2 = G/2 - (3/2)PQ :  {'PASS' if v_i7b else 'FAIL'}")
    OUT[f"N{N}_i7a"] = v_i7a
    OUT[f"N{N}_i7b"] = v_i7b

    # ---- Assemble S^(4) o-space relation --------------------------------
    # Substitute (i2)-(i6) and (i7a/b) into the 15-term identity.
    # Use (i1) to eliminate S_min_o in favour of S_min, P.
    #
    # 15-term identity:
    #   S^(4) = 8 M4 - 8 M3 P - 8 M3e - 4 M3m + 8 M2 P
    #           + 6 S_min P^2 + 4 S_min Q - 4 S_min^2
    #           + 4 T31 + 2 T22 - 2 R3 P - 3 P^2 Q - P^4 - Q^2 - R4
    #
    # Substitutions (with S_min_o = 2 S_min + 3 P^2):
    #
    #   8 M4  = (8/8)[C4 - 6CP - 3(2S+3P^2)^2 + 27(2S+3P^2)P^2 - 27P^4]
    #         = C4 - 6CP - 3(4S^2+12SP^2+9P^4) + 27(2SP^2+3P^4) - 27P^4
    #         = C4 - 6CP - 12S^2 - 36SP^2 - 27P^4 + 54SP^2 + 81P^4 - 27P^4
    #         = C4 - 6CP - 12S^2 + 18SP^2 + 27P^4
    #   (S ≡ S_min)
    #
    #   -8 M3 P = -8P * [C/4 - 3SP - (9/4)P^3] = -2CP + 24SP^2 + 18P^4
    #
    #   -8 M3e = -8*(1/4)*[Ce - 3(2S+3P^2)Q - 3PG + 9P^2 Q]
    #          = -2Ce + 6(2S+3P^2)Q + 6PG - 18P^2 Q
    #          = -2Ce + 12SQ + 18P^2 Q + 6PG - 18P^2 Q
    #          = -2Ce + 12SQ + 6PG
    #
    #   -4 M3m = -4*(1/4)*[Cm - 6GP + 9QP^2] = -Cm + 6GP - 9QP^2
    #
    #   8 M2 P = 8P * [G/2 - (3/2)PQ] = 4GP - 12P^2 Q
    #
    # Collect o-space parts:
    #   C4: +1
    #   C:  -6P - 2P = -8P  (from 8M4 and -8M3P)
    #   Ce: -2
    #   Cm: -1
    #   G:  +6P + 6P + 4P = +16P  (from -8M3e, -4M3m, +8M2P)
    #   (H31 from T31: 4*(1/2)*H31 = 2 H31)
    #   (H22 from T22: 2*(1/2)*H22 = H22)
    #
    # Polynomial in S_min(≡S), P, Q, R3, R4:
    #   from 8M4:       -12S^2 + 18SP^2 + 27P^4
    #   from -8M3P:     +24SP^2 + 18P^4
    #   from -8M3e:     +12SQ
    #   from -4M3m:     -9QP^2
    #   from 8M2P:      -12P^2 Q
    #   from 15-term:   +6SP^2 + 4SQ - 4S^2 - 2R3P - 3P^2Q - P^4 - Q^2 - R4
    #   from T31:       4*(1/2)*(-3PR3) = -6PR3
    #   from T22:       2*(1/2)*(-3Q^2) = -3Q^2
    #
    # Sum S^2 terms:    -12S^2 - 4S^2 = -16 S^2
    # Sum SP^2 terms:   18SP^2 + 24SP^2 + 6SP^2 = 48 SP^2
    # Sum SQ terms:     12SQ + 4SQ = 16 SQ
    # Sum P^4 terms:    27P^4 + 18P^4 - P^4 = 44 P^4
    # Sum P^2Q terms:   -9QP^2 - 12P^2Q - 3P^2Q = -24 P^2 Q
    # Sum Q^2 terms:    -Q^2 - 3Q^2 = -4 Q^2
    # Sum R3P terms:    -2R3P - 6PR3 = -8 R3 P
    # R4:               -R4
    #
    # So the assembled o-space S^(4):
    S4_ospace = (C4
                 - 8 * C * P
                 - 2 * Ce
                 - Cm
                 + 16 * G * P
                 + 2 * H31
                 + H22
                 - 16 * S_min ** 2
                 + 48 * S_min * P ** 2
                 + 16 * S_min * Q
                 + 44 * P ** 4
                 - 24 * P ** 2 * Q
                 - 4 * Q ** 2
                 - 8 * R3 * P
                 - R4)

    v_assembled = (S4_ospace == S4)
    print(f"  assembled S^(4) o-space == direct : {'PASS' if v_assembled else 'FAIL'}")
    OUT[f"N{N}_assembled"] = v_assembled

    local_pass = all(OUT[f"N{N}_{k}"] for k in
                     ["i1", "i2", "i3", "i4", "i5", "i6", "i7a", "i7b", "assembled"])
    OUT[f"N{N}_all"] = local_pass
    all_pass = all_pass and local_pass

OUT["all_pass"] = all_pass

# ---- Print assembled relation ------------------------------------------
print("""
=== Assembled S^(4) o-space relation ===

S^(4) = C4
      - 8*C*P
      - 2*Ce
      - Cm
      + 16*G*P
      + 2*H31
      + H22
      - 16*S_min^2
      + 48*S_min*P^2
      + 16*S_min*Q
      + 44*P^4
      - 24*P^2*Q
      - 4*Q^2
      - 8*R3*P
      - R4
""")

# ---- Coefficient table (exact rationals) --------------------------------
coeff_table = {
    "C4":        "1",
    "C_times_P": "-8",
    "Ce":        "-2",
    "Cm":        "-1",
    "G_times_P": "16",
    "H31":       "2",
    "H22":       "1",
    "S_min^2":   "-16",
    "S_min*P^2": "48",
    "S_min*Q":   "16",
    "P^4":       "44",
    "P^2*Q":     "-24",
    "Q^2":       "-4",
    "R3*P":      "-8",
    "R4":        "-1",
}
OUT["coeff_table"] = coeff_table

elapsed = time.time() - t0
print(f"ALL: {'PASS' if all_pass else 'FAIL'}  ({elapsed:.1f}s)")
OUT["elapsed_s"] = round(elapsed, 2)

outp = Path(__file__).parent / "data" / "s4_ospace_identities.json"
outp.write_text(json.dumps(OUT, indent=2))
print("Saved:", outp)
