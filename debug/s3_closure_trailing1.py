"""S^(3) closure sprint — follow-on 1: trailing-argument-1 evaluator fix (v2).

The stage-A/B evaluator for t3(b1, b2, 1) fed the log-growing odd-harmonic
factor Htilde(j) ~ (ln j)/2 into mpmath nsum(method='levin'), which
silently mis-converges on log-modulated summands (the same disease that
broke the route-B global anchor; see s3_closure_anchor.py — true value
31.5716(12) vs the two wrong stage-1 figures 30.6154 / 30.2197).

Fix: Abel summation over the odd harmonic weight.  With
a_j = (2j+1)^{-b2} 2^{-b1} zeta(b1, j+3/2),

    t3(b1,b2,1) = sum_{j>=1} a_j Htilde(j)
                = sum_{k>=1} (2k-1)^{-1} A(k),    A(k) = sum_{j>=k} a_j

(Tonelli, all terms positive).  Both the inner sum A(1) and the outer
summand A(k)/(2k-1) ~ k^{-(b1+b2-1)} are pure algebraic decay — the
summand class the same nsum machinery already handles at 220-340 dps.

v2 algorithm (the v1 backward-seed design churned: each ceiling
extension recomputed the whole array from a fresh Levin seed):
  - Hurwitz forward recurrence  zeta(b1, (j+1)+3/2) =
    zeta(b1, j+3/2) - (j+3/2)^{-b1}  => one mpf power per term, no
    special-function calls in the bulk loop.
  - Forward cumsum  A(k) = A(1) - sum_{j<k} a_j  with A(1) from ONE
    Levin-safe nsum.  Absolute error stays at the ulp level (positive
    decreasing subtractions), so terms are accurate ABSOLUTELY even
    where relative precision is lost — exactly what the outer Levin
    extrapolation needs.
  - Lazy block extension (recurrence continues), no recomputation.

Gates:
  G1 (rigorous bracket, replaces the v1 old-evaluator comparison):
     direct partial sum of a_j*Htilde(j) to J = 1e5 (strict lower
     bound, positive terms) + explicit integral tail bound.  The new
     value must sit inside; the old cached stage-A/B values are scored
     against the same bracket — quantifying their true errors.
  G2 (stuffle witnesses): R2, R6 at 220 dps, gate 1e-180
     (stage-1 failures: 6.4e-12, 9.8e-20).
  G3 (self-consistency): t3(2,1,1) recomputed at work dps 280 vs 250.
  G4 (the four pending PSLQs): t3(2,1,1) w4 dim-5; t3(4,1,1), t3(3,2,1)
     w6 dim-13; t3(3,4,1) w8 dim-30.  An imprecise value cannot
     PSLQ-identify at 180+ digit residual, so ACCEPTs double as
     precision certificates.

Output: debug/data/s3_closure_trailing1.json (checkpointed per phase).
Run with python -u.
"""
from __future__ import annotations

import json
import sys
import time
from pathlib import Path

import mpmath

HERE = Path(__file__).parent
sys.path.insert(0, str(HERE))
import s3_decomp_pslq as P  # noqa: E402  (monomials, identify, T2/T3, CACHE)

OUT_PATH = HERE / "data" / "s3_closure_trailing1.json"
OUT = {}


def checkpoint():
    OUT_PATH.write_text(json.dumps(OUT, indent=2))


B3ONE = [(2, 1), (2, 3), (3, 2), (4, 1), (4, 3), (3, 4)]  # (b1, b2), b3=1


def lam(s):
    return (1 - mpmath.power(2, -s)) * mpmath.zeta(s)


# ------------------------------------------------------------- fast evaluator
class Trailing1:
    """t3(b1, b2, 1) via Abel summation, recurrence-built terms.

    PRECISION CONTRACT (the v2 bug, caught by the G1 brackets): mpmath's
    nsum re-calls the term function at ELEVATED internal precision when
    the Levin transform escalates; serving terms frozen at the original
    precision lets the transform's cancellations amplify the frozen
    rounding error (observed: 1e-4 .. 1e-11 final errors at 40 dps,
    worst on the slowest series).  Fix: term arrays are cached PER
    WORKING PRECISION and rebuilt at the current precision on demand.
    """

    NTERMS0 = 4096
    KMAX_HARD = 4_000_000

    def __init__(self, b1, b2):
        self.b1, self.b2 = b1, b2
        self._cache = {}            # prec -> state dict
        self.max_k_requested = 0
        self.rebuilds = 0

    def _a_direct(self, j):
        j = int(j)
        return (mpmath.mpf(2 * j + 1) ** -self.b2 * mpmath.power(2, -self.b1)
                * mpmath.zeta(self.b1, mpmath.mpf(j) + mpmath.mpf(3) / 2))

    def _new_state(self):
        # everything at the CURRENT working precision
        st = {
            "A1": mpmath.nsum(self._a_direct, [1, mpmath.inf],
                              method='levin'),
            "z": mpmath.zeta(self.b1, mpmath.mpf(5) / 2),   # j = 1
            "j_done": 1,
            "cum": mpmath.mpf(0),       # sum_{j <= j_done-1} a_j
            "terms": [],                # term[k-1] = A(k)/(2k-1)
        }
        self.rebuilds += 1
        return st

    def _fill(self, st, k_target):
        if k_target > self.KMAX_HARD:
            raise RuntimeError("outer nsum requested k=%d > hard cap"
                               % k_target)
        terms = st["terms"]
        while len(terms) < k_target:
            k = len(terms) + 1               # invariant: j_done == k,
            A_k = st["A1"] - st["cum"]       # cum == sum_{j<k} a_j,
            terms.append(A_k / (2 * k - 1))  # z == zeta(b1, k+3/2)
            j = st["j_done"]
            a_j = (mpmath.mpf(2 * j + 1) ** -self.b2
                   * mpmath.power(2, -self.b1) * st["z"])
            st["cum"] += a_j
            st["z"] -= (mpmath.mpf(j) + mpmath.mpf(3) / 2) ** -self.b1
            st["j_done"] += 1

    def term(self, k):
        k = int(k)
        prec = mpmath.mp.prec
        st = self._cache.get(prec)
        if st is None:
            st = self._new_state()
            self._fill(st, self.NTERMS0)
            self._cache[prec] = st
        if k > len(st["terms"]):
            self._fill(st, 2 * k)
        if k > self.max_k_requested:
            self.max_k_requested = k
        return st["terms"][k - 1]

    def value(self):
        return mpmath.nsum(self.term, [1, mpmath.inf], method='levin')


# ------------------------------------------------------------------ G1
print("=== G1: rigorous direct-sum brackets at 40 dps ===", flush=True)
mpmath.mp.dps = 40
G1 = {}
NEW40 = {}
J_DIRECT = 100_000
for b1, b2 in B3ONE:
    t0 = time.time()
    # direct partial sum via recurrences (strict lower bound)
    z = mpmath.zeta(b1, mpmath.mpf(5) / 2)
    H = mpmath.mpf(1)            # Htilde(1) = 1
    pref = mpmath.power(2, -b1)
    partial = mpmath.mpf(0)
    for j in range(1, J_DIRECT + 1):
        a_j = mpmath.mpf(2 * j + 1) ** -b2 * pref * z
        partial += a_j * H
        z -= (mpmath.mpf(j) + mpmath.mpf(3) / 2) ** -b1
        H += mpmath.mpf(1) / (2 * j + 1)   # Htilde(j+1)
    # rigorous tail bound: a_j <= 2^{-b2-b1} Cz (j+1/2)^{1-b1}/(b1-1) j^{-b2},
    # Htilde(j) <= (ln j)/2 + 1.5;  bound by integral of x^{-m}(ln x/2 + c),
    # m = b1+b2-1 >= 2, integrand decreasing:
    #   int_J^inf x^{-m}(ln x/2 + c) dx
    #     = J^{1-m} [ (ln J/2 + c)/(m-1) + 1/(2(m-1)^2) ]
    m_exp = b1 + b2 - 1
    Cz = 1 + mpmath.mpf(b1 - 1) / (J_DIRECT + mpmath.mpf(1) / 2)
    pref_tail = (mpmath.power(2, -b1) * mpmath.power(2, -b2)
                 * Cz / (b1 - 1))
    c = mpmath.mpf(3) / 2
    J = mpmath.mpf(J_DIRECT)
    tail = (pref_tail * J ** (1 - m_exp)
            * ((mpmath.log(J) / 2 + c) / (m_exp - 1)
               + mpmath.mpf(1) / (2 * (m_exp - 1) ** 2)))
    ev = Trailing1(b1, b2)
    v_new = ev.value()
    NEW40[(b1, b2)] = v_new
    in_bracket = partial <= v_new <= partial + tail
    # score the old cached value against the same bracket
    old_key = "t3(%d,%d,1)@220" % (b1, b2)
    old = mpmath.mpf(P.CACHE[old_key]) if old_key in P.CACHE else None
    old_err = (mpmath.nstr(old - v_new, 6) if old is not None else "n/a")
    old_in = (partial <= old <= partial + tail) if old is not None else None
    G1["t3(%d,%d,1)" % (b1, b2)] = {
        "lower": mpmath.nstr(partial, 25), "tail_bound": mpmath.nstr(tail, 4),
        "new_value": mpmath.nstr(v_new, 25), "new_in_bracket": bool(in_bracket),
        "old_minus_new": old_err, "old_in_bracket": old_in,
        "outer_terms_used": ev.max_k_requested}
    print("  t3(%d,%d,1): new=%s in-bracket=%s | old-new=%s old-in=%s "
          "(outer k<=%d, %.0fs)"
          % (b1, b2, mpmath.nstr(v_new, 18), in_bracket, old_err, old_in,
             ev.max_k_requested, time.time() - t0), flush=True)
G1_PASS = all(v["new_in_bracket"] for v in G1.values())
OUT["G1"] = G1
OUT["G1_pass"] = G1_PASS
print("  G1:", "PASS" if G1_PASS else "FAIL", flush=True)
checkpoint()

# ----------------------------------------------------------- high-precision
DPS_A, DPS_B = 220, 340
NEW = {}
for tgt_dps, keys in ((DPS_A, B3ONE), (DPS_B, [(3, 4)])):
    mpmath.mp.dps = tgt_dps + 30
    for b1, b2 in keys:
        key = "t3(%d,%d,1)@%d" % (b1, b2, tgt_dps)
        t0 = time.time()
        ev = Trailing1(b1, b2)
        v = ev.value()
        NEW[key] = mpmath.nstr(v, tgt_dps - 5)
        print("  [%s recomputed: outer k<=%d, %.0fs]"
              % (key, ev.max_k_requested, time.time() - t0), flush=True)
        OUT.setdefault("outer_terms", {})[key] = ev.max_k_requested
        checkpoint()

# G3: self-consistency at a different working precision
mpmath.mp.dps = DPS_A + 60
v_alt = Trailing1(2, 1).value()
mpmath.mp.dps = DPS_A + 60
d3 = abs(v_alt - mpmath.mpf(NEW["t3(2,1,1)@%d" % DPS_A]))
OUT["G3_two_dps_diff"] = mpmath.nstr(d3, 4)
G3_PASS = bool(d3 < mpmath.mpf(10) ** -(DPS_A - 10))
OUT["G3_pass"] = G3_PASS
print("=== G3: two-working-dps self-consistency: %s -> %s ==="
      % (mpmath.nstr(d3, 4), "PASS" if G3_PASS else "FAIL"), flush=True)
checkpoint()

# ------------------------------------------------------------- cache rewrite
OUT["old_cache_values"] = {
    k: P.CACHE[k] for k in list(P.CACHE)
    if k.startswith("t3(") and ",1)@" in k}
for k, v in NEW.items():
    P.CACHE[k] = v
P.CACHE_PATH.write_text(json.dumps(P.CACHE, indent=1))
OUT["new_cache_values"] = NEW
print("Cache overwritten with %d corrected values." % len(NEW), flush=True)
checkpoint()

# ------------------------------------------------------------------ G2
print("=== G2: stuffle witness gates R2, R6 (dps %d) ===" % DPS_A, flush=True)
mpmath.mp.dps = DPS_A
t3v = {(b1, b2, 1): mpmath.mpf(NEW["t3(%d,%d,1)@%d" % (b1, b2, DPS_A)])
       for b1, b2 in B3ONE}
for k in [(2, 1, 3), (4, 1, 3)]:
    t3v[k] = P.T3(*k, dps=DPS_A)
t2v = {k: P.T2(*k, dps=DPS_A)
       for k in [(2, 1), (2, 4), (4, 1), (4, 4), (5, 1), (7, 1)]}
R2 = lam(3) * t2v[(2, 1)] - (t3v[(3, 2, 1)] + t3v[(2, 3, 1)]
                             + t3v[(2, 1, 3)] + t2v[(5, 1)] + t2v[(2, 4)])
R6 = lam(3) * t2v[(4, 1)] - (t3v[(3, 4, 1)] + t3v[(4, 3, 1)]
                             + t3v[(4, 1, 3)] + t2v[(7, 1)] + t2v[(4, 4)])
gate = mpmath.mpf(10) ** -(DPS_A - 40)
OUT["G2_R2_residual"] = mpmath.nstr(abs(R2), 4)
OUT["G2_R6_residual"] = mpmath.nstr(abs(R6), 4)
G2_PASS = bool(abs(R2) < gate and abs(R6) < gate)
OUT["G2_pass"] = G2_PASS
print("  R2 residual: %s (was 6.4e-12)" % mpmath.nstr(abs(R2), 4), flush=True)
print("  R6 residual: %s (was 9.8e-20)" % mpmath.nstr(abs(R6), 4), flush=True)
print("  G2:", "PASS" if G2_PASS else "FAIL", flush=True)
checkpoint()

# ------------------------------------------------------------------ G4: PSLQ
print("=== G4: the four pending PSLQs ===", flush=True)
mpmath.mp.dps = DPS_A
ln2, pi = mpmath.log(2), mpmath.pi
z3, z5 = mpmath.zeta(3), mpmath.zeta(5)
Li = {k: mpmath.polylog(k, mpmath.mpf(1) / 2) for k in (4, 5, 6)}
t251 = P.T2(5, 1, DPS_A)

gens4 = [("ln2", 1, ln2), ("pi2", 2, pi ** 2), ("z3", 3, z3),
         ("Li4", 4, Li[4])]
b4 = P.monomials(gens4, 4)
OUT["t3(2,1,1)"] = P.identify("t3(2,1,1)", t3v[(2, 1, 1)], b4, DPS_A)

gens6 = gens4 + [("z5", 5, z5), ("Li5", 5, Li[5]), ("Li6", 6, Li[6]),
                 ("t2(5,1)", 6, t251)]
b6 = P.monomials(gens6, 6)
OUT["t3(4,1,1)"] = P.identify("t3(4,1,1)", t3v[(4, 1, 1)], b6, DPS_A)
OUT["t3(3,2,1)"] = P.identify("t3(3,2,1)", t3v[(3, 2, 1)], b6, DPS_A)
checkpoint()

print("--- w8 at %d dps ---" % DPS_B, flush=True)
mpmath.mp.dps = DPS_B
ln2, pi = mpmath.log(2), mpmath.pi
z3, z5, z7 = mpmath.zeta(3), mpmath.zeta(5), mpmath.zeta(7)
Li = {k: mpmath.polylog(k, mpmath.mpf(1) / 2) for k in range(4, 9)}
gens8 = [("ln2", 1, ln2), ("pi2", 2, pi ** 2), ("z3", 3, z3),
         ("Li4", 4, Li[4]), ("z5", 5, z5), ("Li5", 5, Li[5]),
         ("Li6", 6, Li[6]), ("t2(5,1)", 6, P.T2(5, 1, DPS_B)),
         ("z7", 7, z7), ("Li7", 7, Li[7]),
         ("Li8", 8, Li[8]), ("t2(7,1)", 8, P.T2(7, 1, DPS_B)),
         ("t2(5,3)", 8, P.T2(5, 3, DPS_B))]
b8 = P.monomials(gens8, 8)
OUT["t3(3,4,1)"] = P.identify(
    "t3(3,4,1)", mpmath.mpf(NEW["t3(3,4,1)@%d" % DPS_B]), b8, DPS_B)

G4_PASS = all(OUT["t3(%d,%d,1)" % (b1, b2)] is not None
              for b1, b2 in [(2, 1), (4, 1), (3, 2), (3, 4)])
OUT["G4_all_identified"] = bool(G4_PASS)
print("  G4 (all four identified):", "PASS" if G4_PASS else "FAIL", flush=True)

OUT["all_gates_pass"] = bool(G1_PASS and G2_PASS and G3_PASS and G4_PASS)
checkpoint()
print("\nALL GATES:", "PASS" if OUT["all_gates_pass"] else "FAIL", flush=True)
print("Saved:", OUT_PATH, flush=True)
