"""S^(3) round-2 decomposition — Step 1: exact setup verification (Fractions).

S^(3) (Paper 28 eq:three_loop_factorized, a=4, p=1, d_q^T/mu_q = 1) is the
chain-topology three-loop sum

    S3 = sum_{n1,n2,n3} [sum_{q1} W(n1,n2,q1)] [sum_{q2} W(n2,n3,q2)]
         phi(n1) phi(n2) phi(n3),
    phi(n) = g_n/|lambda_n|^4 = 2/(n+3/2)^2 - (1/2)/(n+3/2)^4.

Verified bit-exactly in Fraction arithmetic at n_max <= 12:

  V1. Channel-count collapse (round-1c closed form, PRODUCTION W):
      sum_q so4_channel_count(n1,n2,q) = I(n1,n2) = 2 min(n1,n2) - 1 - [n1=n2]
      for n1,n2 >= 1 (0 if min = 0), so
      S3(N) [direct CG] == sum I(n1,n2) I(n2,n3) phi phi phi  [collapsed].

  V2. Nine-term expansion (same box [1,N]^3):
      S3 = 4 M3 - 4 S_min P - 4 M2 + P^3 + 2 P Q + R3,
      M3 = sum min(n1,n2) min(n2,n3) phi1 phi2 phi3,
      M2 = sum min(n1,n2) phi1 phi2^2,   S_min = sum min phi phi,
      P = sum phi, Q = sum phi^2, R3 = sum phi^3.

  V3. o-space reductions (o = 2n+3, psi(o) = 8(o^2-1)/o^4 = phi(n)):
      M3 = C/4 - 3 S_min P - (9/4) P^3,   C = sum_{o>=5} min12 min23 psi^3
      M2 = G/2 - (3/2) P Q,               G = sum_{o>=5} min psi1 psi2^2

  V4. Factorized closed-form anchor (basis of the high-precision run):
      S3(N) = sum_m phi(m) L_N(m)^2,  L_N(m) = 2 A_N(m) - P_N - phi(m),
      A_N(m) = sum_{n<=N} min(n,m) phi(n)  -- box-consistent identity.

Output: prints PASS/FAIL per check; writes debug/data/s3_decomp_setup.json
"""
from __future__ import annotations

import json
import sys
import time
from fractions import Fraction
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from geovac.qed_vertex import so4_channel_count  # noqa: E402

N = 12
OUT = {}
t0 = time.time()


def phi(n: int) -> Fraction:
    x = Fraction(2 * n + 3, 2)
    return 2 / x ** 2 - Fraction(1, 2) / x ** 4


def I_closed(n1: int, n2: int) -> int:
    if min(n1, n2) == 0:
        return 0
    return 2 * min(n1, n2) - 1 - (1 if n1 == n2 else 0)


# ---------------------------------------------------------------- V1
print("V1: channel-count collapse against production so4_channel_count ...")
I_direct = {}
bad = []
for n1 in range(0, N + 1):
    for n2 in range(0, N + 1):
        s = sum(so4_channel_count(n1, n2, q) for q in range(1, n1 + n2 + 1))
        I_direct[(n1, n2)] = s
        if s != I_closed(n1, n2):
            bad.append((n1, n2, s, I_closed(n1, n2)))
v1a = not bad
print("  I(n1,n2) closed form, n<=%d: %s" % (N, "EXACT" if v1a else bad[:5]))

PH = {n: phi(n) for n in range(0, N + 1)}
S3_direct = sum(I_direct[(a, b)] * I_direct[(b, c)] * PH[a] * PH[b] * PH[c]
                for a in range(N + 1) for b in range(N + 1)
                for c in range(N + 1))
S3_collapsed = sum(I_closed(a, b) * I_closed(b, c) * PH[a] * PH[b] * PH[c]
                   for a in range(1, N + 1) for b in range(1, N + 1)
                   for c in range(1, N + 1))
v1b = (S3_direct == S3_collapsed)
print("  S3 direct-CG == collapsed (N=%d): %s" % (N, v1b))
print("  S3(N=%d) = %s ~ %.10f" % (N, S3_direct, float(S3_direct)))
OUT["V1_closed_form_exact"] = v1a
OUT["V1_collapse_bitexact"] = v1b
OUT["S3_partial_N12"] = str(S3_direct)

# ---------------------------------------------------------------- V2
print("V2: nine-term expansion (box [1,N]^3) ...")
rng = range(1, N + 1)
P = sum(PH[n] for n in rng)
Q = sum(PH[n] ** 2 for n in rng)
R3 = sum(PH[n] ** 3 for n in rng)
S_min = sum(min(a, b) * PH[a] * PH[b] for a in rng for b in rng)
M2 = sum(min(a, b) * PH[a] * PH[b] ** 2 for a in rng for b in rng)
M3 = sum(min(a, b) * min(b, c) * PH[a] * PH[b] * PH[c]
         for a in rng for b in rng for c in rng)
rhs = 4 * M3 - 4 * S_min * P - 4 * M2 + P ** 3 + 2 * P * Q + R3
v2 = (rhs == S3_collapsed)
print("  S3 == 4 M3 - 4 S_min P - 4 M2 + P^3 + 2PQ + R3: %s" % v2)
OUT["V2_nine_term_bitexact"] = v2

# ---------------------------------------------------------------- V3
print("V3: o-space reductions ...")
O = 2 * N + 3
odds = range(5, O + 1, 2)


def psi(o: int) -> Fraction:
    return Fraction(8 * (o * o - 1), o ** 4)


C = sum(min(a, b) * min(b, c) * psi(a) * psi(b) * psi(c)
        for a in odds for b in odds for c in odds)
G = sum(min(a, b) * psi(a) * psi(b) ** 2 for a in odds for b in odds)
v3a = (M3 == C / 4 - 3 * S_min * P - Fraction(9, 4) * P ** 3)
v3b = (M2 == G / 2 - Fraction(3, 2) * P * Q)
print("  M3 == C/4 - 3 S_min P - (9/4)P^3: %s" % v3a)
print("  M2 == G/2 - (3/2) P Q: %s" % v3b)
v3c = (S3_collapsed
       == C - 2 * G - 16 * S_min * P - 8 * P ** 3 + 8 * P * Q + R3)
print("  S3 == C - 2G - 16 S_min P - 8 P^3 + 8 P Q + R3: %s" % v3c)
OUT["V3_M3_reduction"] = v3a
OUT["V3_M2_reduction"] = v3b
OUT["V3_S3_assembled"] = v3c

# ---------------------------------------------------------------- V4
print("V4: factorized closed-form anchor identity ...")
S3_fact = Fraction(0)
for m in rng:
    A = sum(min(n, m) * PH[n] for n in rng)
    L = 2 * A - P - PH[m]
    S3_fact += PH[m] * L * L
v4 = (S3_fact == S3_collapsed)
print("  S3 == sum_m phi(m) [2 A(m) - P - phi(m)]^2: %s" % v4)
OUT["V4_factorized_anchor"] = v4

ok = all(v for k, v in OUT.items() if isinstance(v, bool))
print("\nALL CHECKS: %s  (%.1fs)" % ("PASS" if ok else "FAIL",
                                     time.time() - t0))
OUT["all_pass"] = ok
outp = Path(__file__).parent / "data" / "s3_decomp_setup.json"
outp.write_text(json.dumps(OUT, indent=2))
print("Saved:", outp)
