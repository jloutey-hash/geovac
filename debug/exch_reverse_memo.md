# Exchange order exactly 9 for the reverse ordering a > b (the l > 0 orbital carries the larger exponent)

Date 2026-11-19. Derivation memo; no paper, test, doc, CLAUDE.md, CHANGELOG or git state touched; nothing SEALED
opened. Continues `debug/pachucki_extension_memo.md`, `debug/pachucki_m_nonzero_memo.md`,
`debug/exact9_independence_memo.md`, `debug/exact9_m_memo.md`.

Convention throughout (as in the derivation): a = exponent of the l > 0 orbital on A, b = exponent of 1s on B,
t = b/a. The four mixed integrals: 2psig = (2p0_A 1s_B|2p0_A 1s_B), 3dsig = (3d0_A 1s_B|3d0_A 1s_B),
2ppi = (2p+1_A 1s_B|2p-1_A 1s_B), 3ddel = (3d+2_A 1s_B|3d-2_A 1s_B). Covered before: t > 1. Here: t in (0, 1).

## Verdict: GO, all four integrals

| integral | span membership, a > b | order exactly 9, every a > b > 0 | hp agreement (digits) at (2,1) and (3/2,4/5) |
|---|---|---|---|
| 2p_sigma_A 1s_B | yes, same nine functions | yes | 48-49 (T24) |
| 3d_0_A 1s_B | yes | yes | 42.4, 40.4 (T18; tau-truncation limited, see table) |
| 2p_{+-1,A} 1s_B | yes | yes | 48-48.5 (T24) |
| 3d_{+-2,A} 1s_B | yes | yes | 48-48.7 (T24) |

With the earlier t > 1 results, the statement now holds for every a != b > 0. No exceptional ratio in (0, 1); no
order drop. The existing Q(t) closed forms are valid on (0, 1) **unchanged**, provided the nine functions are read
with their real branch: Ei of a positive argument = principal value, and LC = ln((a+b)/|b-a|).

## Q1. Is the Q(t) closed form valid for t in (0, 1)?

Yes. Three independent arguments, each sufficient.

**(i) Analytic continuation (a proof; uses only the t > 1 result).** For fixed R > 0, K(R; 1, t) is a parameter
integral of e^{-(r1A + t r1B + r2A + t r2B)} x polynomial / r12, absolutely and locally uniformly convergent for
Re t > 0, hence holomorphic in t on the right half-plane, in particular across t = 1. The closed form
F(t) = sum_i c_i(t, R) g_i(t, R) equals K on (1, oo). The c_i are rational in t with poles only at t = 0, +-1. Each
g_i continues holomorphically from (1, oo) into the upper half-plane near the real axis (Ei(z) = gamma + Log z +
Ein z with Ein entire; LC = Log(1+t) - Log(t-1)). By the identity theorem F = K on that region, so for real
s in (0, 1), K(s) = lim_{eps -> 0+} F(s + i eps). Boundary values, key by key (t -> s + i0, s < 1):
* Ei(-(1+t)R), Ei(-2(1+t)R) (in g2, g6, g7): argument stays negative; no change.
* psi_b = e^{-2tR}[Ei(2(t-1)R) - 2Ei((t-1)R) - LC]: the two Ei arguments go from positive to negative through the
  upper half-plane, each picks up +i pi relative to the real Ei: +i pi - 2 i pi = -i pi. LC: Log(t-1) -> ln|1-s| + i pi,
  so -LC picks up +i pi. Net 0.
* psi_a = e^{-2R}[Ei(2(1-t)R) - 2Ei((1-t)R) - LC]: Ei arguments go from negative to positive through the lower
  half-plane; with Ei(w) = gamma + Log(-w) + Ein(w) for w < 0, Log(-w) picks up +i pi each: +i pi - 2 i pi; -LC: +i pi.
  Net 0.
* f0 = int_R^oo FM1(s)/s ds, FM1 = c(g7 - g6 - g9 + g8): continuation commutes with the integral (uniform
  exponential decay on compacts), and g6..g9 have real boundary values by the above.
So every one of the nine functions has as boundary value its **real-branch** form, and K = sum c_i(s, R) g_i^real on
(0, 1). The "LC tie" (the -LC inside psi, with the coefficient of LC tied to that of Ei(2(a-b)R)) is exactly what
makes the i pi's cancel; this is why the bracket is branch-free.

**(ii) Decay, key by key, and uniqueness at t < 1 (`exch_reverse_decay.py`).** Every function key of every
intermediate jet of the frozen Q(t) solves (audit lists in `tests/data/paper58_pachucki_ext_2026_11_18.json` and
`tests/data/paper58_exchange_m_2026_11_18.json`) has rate linear in t; the same 12 keys + f0 in all four cases:

| key | rate at R -> oo | decays on (0,1) | Ei argument sign t<1 / t>1 |
|---|---|---|---|
| e^{-2tR}, e^{-2tR} LC | -2t | yes | - |
| e^{-2R}, e^{-2R} LC | -2 | yes | - |
| e^{-(1+t)R} | -(1+t) | yes | - |
| e^{-2tR} Ei(2(t-1)R) | -2 | yes | - / + |
| e^{-2tR} Ei((t-1)R) | -(1+t) | yes | - / + |
| e^{-2R} Ei(2(1-t)R) | -2t | yes | + / - |
| e^{-2R} Ei((1-t)R) | -(1+t) | yes | + / - |
| Ei(-(1+t)R) | -(1+t) | yes | - / - |
| e^{2tR} Ei(-2(1+t)R) | -2 | yes | - / - |
| e^{2R} Ei(-2(1+t)R) | -2t | yes | - / - |
| f0 | -2 min(1, t) | yes | - |

All decay on the whole of (0, 1) (and of (1, oo)); Ei arguments vanish only at t = 1 (or t = -1). The decay verdict
does not change across t = 1, although four Ei arguments change sign (a positive argument is a principal value and
grows like e^{mu R}/(mu R); the key's rate lam + mu is still negative). Denominators of all intermediates: t, t - 1,
t + 1 only. Hence the Q(t) residual identities (identically zero, existing memos) specialise at every t0 in (0, 1);
the derivative rule Ei(mu R)' = e^{mu R}/R holds for the principal value too, and LC is a constant; the homogeneous
solutions (A + B ln r for (r u')' = src; A + B r for u'' = src) do not decay; so the specialised jets are the decaying
solutions **provided the sources are right at t0 < 1**. The sources use Pachucki's ln|x| and real Ei; (i) and (iii)
confirm them end to end.

**(iii) Fresh solves over Q at a > b (`exch_reverse_closed.py`), not from the stored form.** The reviewed helpers
(`tests/_pachucki_extension.solve_all`, `tests/_exchange_m_nonzero.solve_case`) were run over Q at
(a, b) = (2, 1), (3/2, 4/5), (3, 1), (7, 3), i.e. t = 1/2, 8/15, 1/3, 3/7, with every `solve_L0` / `solve_D2` wrapped
so that its output is re-checked as an exact identity (L0 u = src, resp. u'' = src) and its free unknowns counted:

| case | solves per pair | residual failures | free unknowns (max) | non-decaying keys | K in nine-span | fresh = stored Q(t) at t = b/a (rescaled), exact |
|---|---|---|---|---|---|---|
| 2psig | 143 | 0 | 0 | none | yes, all nine present | yes, 4/4 pairs |
| 3dsig | 575 | 0 | 0 | none | yes, all nine | yes, 4/4 |
| 2ppi | 159 | 0 | 0 | none | yes, all nine | yes, 4/4 |
| 3ddel | 735 | 0 | 0 | none | yes, all nine | yes, 4/4 |

Rescaling used, for the coefficient of R^p g_i: a^{p+1} c_{i,p}(1, t) for i = 2..9 and a^{p+3} c_{0,p}(1, t)
for f0 (K(R; a, b) = a K(aR; 1, t), f0(R; a, b) = a^{-2} f0(aR; 1, t)). The exact equality is coordinate by
coordinate over Q. (t = 1/3 is the coincidence -(a+b) = 2(b-a) of Ei arguments named in the independence memo; it
changes nothing.) Solve times per pair: 2psig ~6 s, 2ppi ~8 s, 3dsig ~40 s, 3ddel ~60 s.

**Numerical agreement (closed form, real branch, dps 60, f0 by one 1-D quadrature) vs `geovac.qfd_core.exchange_hp`
(Neumann sum over tau, dps 50; never sees the closed form)** (`exch_reverse_compare.py`):

| case | (a, b) | R | closed form | exchange_hp | abs diff | digits | hp last term (tau_max) |
|---|---|---|---|---|---|---|---|
| 2psig | (2,1) | 7/5 | 0.21367060454688371794197541596666500058880902378683 | 0.213670604546883717941975415966665000588809023787 | 1.7e-49 | 48.1 | (24) |
| 2psig | (2,1) | 3 | 0.040608113546844355350888472346160800646877846053398 | 0.0406081135468443553508884723461608006468778460534 | 2.3e-51 | 49.2 | (24) |
| 2psig | (3/2,4/5) | 23/10 | 0.14010439248053493905322062816435580112334122300457 | 0.140104392480534939053220628164355801123341223005 | 4.3e-49 | 47.5 | (24) |
| 3dsig | (2,1) | 7/5 | 0.12342591693187075402340762686153857506668050026616 | 0.1234259169318707540234076268615385750666804454 | 5.5e-44 | 42.4 | 1.1e-40 (18) |
| 3dsig | (3/2,4/5) | 23/10 | 0.1071294410961232734366125770198190175104772281078 | 0.1071294410961232734366125770198190175104731157 | 4.1e-42 | 40.4 | 6.3e-39 (18) |
| 2ppi | (2,1) | 7/5 | -0.053414407564864070731721969857544101833810078199181 | -0.0534144075648640707317219698575441018338100781992 | 1.9e-50 | 48.5 | -5.1e-65 (24) |
| 2ppi | (2,1) | 3 | -0.005221878549229800645919110448752139698383675023487 | -0.00522187854922980064591911044875213969838367502349 | 3.0e-51 | 48.2 | -1.4e-51 (24) |
| 2ppi | (3/2,4/5) | 23/10 | -0.025561712337534865064785475286099417236581602558428 | -0.0255617123375348650647854752860994172365816025584 | 2.8e-50 | 48.0 | -1.4e-62 (24) |
| 3ddel | (2,1) | 7/5 | 0.021538789173312192283418778011746166652730919221796 | 0.0215387891733121922834187780117461666527309192218 | 4.3e-51 | 48.7 | 2.4e-62 (24) |
| 3ddel | (2,1) | 3 | 0.0026715017576626612511866741357537887519491336736568 | 0.00267150175766266125118667413575378875194913367366 | 3.2e-51 | 47.9 | 2.5e-49 (24) |
| 3ddel | (3/2,4/5) | 23/10 | 0.010531578613204001824969129654174240375705081700188 | 0.0105315786132040018249691296541742403757050817002 | 1.2e-50 | 48.0 | 5.8e-60 (24) |

Extra rows at tau_max = 18 (3ddel, both pairs) agree to 42.6 and 40.6 digits, i.e. to the size of the dropped tau
tail, and converge to the T24 values: the residual is hp truncation, not the closed form. The 3dsig T24 runs hit
the 60-min cap and were killed; the 3dsig record is the T18 rows, which clear 40 digits and whose differences match
the size of the dropped tail (last term x ~1e-3).
The closed forms at (2,1) are not the (1,2) values (2psig: 0.213671 at R = 7/5 vs 0.220068 at (1,2), R = 1.4), so
these are genuinely new points, not a relabelling.

## Q2. The function list for a > b

Unchanged: f0, Ei(-(a+b)R), e^{-(a+b)R}, e^{-2aR}, e^{-2bR}, e^{2aR}Ei(-2(a+b)R), e^{2bR}Ei(-2(a+b)R),
psi_a = e^{-2aR}[Ei(2(a-b)R) - 2Ei((a-b)R) - LC], psi_b = e^{-2bR}[Ei(2(b-a)R) - 2Ei((b-a)R) - LC], with
**LC = ln((a+b)/|b-a|) = ln((a+b)/(a-b))** and principal-value Ei for the (now positive) arguments of psi_a. f0 is
defined by the same label-60 F^(-1) with the same LC (the helper already writes it as ln|(w-u)/(w+u)| = -LC).
No psi swap is needed in the definitions; what swaps is the asymptotic class: for a > b, psi_a ~ e^{-2bR}/(2(a-b)R)
and psi_b ~ -LC e^{-2bR}, mirror images of a < b. The Ei-argument coincidence -(a+b) = 2(b-a) moves from b = 3a to
a = 3b (t = 1/3; covered by the (3,1) fresh solve). The derivative table (`tests/_exact_order.hand_table`) is
branch-free and identical.

## Q3. Independence over C(R) for a > b

The lemma's proof is symmetric in a <-> b and goes through verbatim: Step 1 log coefficients at R = 0 are
A_2 = 1, A_6 = e^{2aR}, A_7 = e^{2bR}, A_8 = -e^{-2aR}, A_9 = -e^{-2bR}, A_1 = -2c[Cinh(2bR) - Cinh(2aR)], A_3..5 = 0
(ln|x| has the same log coefficient for either sign of x; the real branch is what enters); Step 3 needs the five
rates {0, +-2a, +-2b} distinct and c = 1/(2(a^2 - b^2)) != 0 (both hold for a != b); Step 4 needs {a+b, 2a, 2b}
distinct (holds for a != b). Numerical check of Step 1 with the real branch (`exch_reverse_lemma.py`, dps 50, fifth
difference on 1e-3..6e-3) at (2,1), (3/2,4/5), (3,1), (sqrt 5, 1): every g_i - A_i ln R is at the analytic floor
(<= 8e-12, the same size as for the pure exponentials) while the removed ln parts are 0.087 (g2, g6..g9) and 1.2e-7 (f0,
whose A_1 = O(R^2)). Confirmed.

## Q4. Determinant leg on (0, 1)

`exch_reverse_det.py`: input = the frozen Q(t) closed form (a = 1, b = t). Algorithm A: rows K..K^(8) by the helper
calculus, cleared by t^i (t-1)^j (t+1)^k, det at rational R by fraction-free Bareiss at integer t + exact Newton
interpolation (checked on extra nodes), R in {1, 2, 7/5, 1/3, 5}. Algorithm B: rows by the hand connection matrix
(checked equal to A's rows identically over Q(t)), entries at R cleared by their lcm, det by sympy DomainMatrix over
ZZ[t], R in {3/2, 11/7, 9}. Residual = det numerator with t, t - 1, t + 1 divided out.

| case | det degree | mult t, t-1, t+1 | residual degree | residual real roots in (0,1) per R (A: 1,2,7/5,1/3,5; B: 3/2,11/7,9) | gcd over 8 R | gcd roots in (0,1) | max pairwise gcd degree |
|---|---|---|---|---|---|---|---|
| 2psig | 151 | 59, 33, 22 | 37 | 1,3,2,1,1; 3,3,1 | constant | 0 | 0 |
| 3dsig | 216 | 91, 42, 28 | 55 | 3,1,1,1,2; 1,1,0 | constant | 0 | 0 |
| 2ppi | 146 | 59, 34, 22 | 31 | 1,1,1,0,1; 1,1,1 | constant | 0 | 0 |
| 3ddel | 206 | 91, 44, 28 | 43 | 1,1,1,1,1; 1,1,2 | constant | 0 | 0 |

Degrees and multiplicities agree with `exact9_independence_memo.md` / `exact9_m_memo.md` (same determinant), and
both algorithms agree at every R. The gcd of the residuals is a constant: even any two R values share no root. The
residuals do have real roots in (0, 1) at each single R (as in (1, oo)); they move with R, so they are zero curves of
det(t, R), not order drops. Since no entry has a pole in (0, 1) (denominators t, t +- 1), det M(t0, R) is not
identically zero for any t0 in (0, 1). With (L) this gives minimal order exactly 9 for every t in (0, 1).

## Q5. A symmetry a > b <-> a < b?

None that is clean, and none needed. The family is not symmetric under a <-> b (values above). Swapping the centres
(reflection z -> R - z) maps (phi_A(a) 1s_B(b)|...) to (1s_A(b) phi_B(a)|...): the l > 0 orbital still carries a, so
the ratio t is unchanged; swapping electron labels changes nothing for an exchange integral. The only link between
the two orderings is the analytic continuation in t of Q1(i), through the upper (or lower) half-plane around t = 1:
the same rational coefficients, the real branch of the same nine functions.

## Where this could be wrong / not established

1. Q1(ii) relies on Pachucki's source formula (ln|x|, real Ei) being valid at a > b; that is not proved here
   separately, but Q1(i) does not need it, and Q1(iii)'s fresh solves + the 40-49 digit agreement test it end to end.
2. Q1(i) is a short hand argument (holomorphy of the parameter integral in t, identity theorem, boundary values of
   Log); unreviewed. The i pi bookkeeping is the whole content: one sign slip would leave +-2 i pi e^{-2aR} terms,
   which the real numerical agreement would expose (it does not).
3. Exactly 9 also needs the operator to exist (order <= 9): the span is closed under d/dR (branch-free table), so
   it holds at t < 1 as at t > 1.
4. t = 1 is excluded (order 4, the equal-exponent collapse, already covered by the limit).
5. The 3dsig numerical rows use tau_max 18 (40.4 and 42.4 digits, truncation-limited); the T24 runs hit the 60-min
   cap. A third R for 3dsig was not obtained (the other three cases have three (pair, R) points at 47.5-49 digits).

## Files written

Scripts: `debug/exch_reverse_core.py` (real-branch evaluator, stored-form specialisation with scaling, wrapped
fresh solver), `debug/exch_reverse_closed.py` (Q1 (iii)), `debug/exch_reverse_hp.py` (exchange_hp references),
`debug/exch_reverse_compare.py`, `debug/exch_reverse_decay.py` (Q1 (ii)), `debug/exch_reverse_lemma.py` (Q3),
`debug/exch_reverse_det.py` (Q4). Data: `debug/data/exch_reverse_closed_<case>_<a>_<b>.json`,
`debug/data/exch_reverse_hp_<case>_<a>_<b>_<R>_T<tau>.json`, `debug/data/exch_reverse_compare.json`,
`debug/data/exch_reverse_decay.json`, `debug/data/exch_reverse_lemma.json`, `debug/data/exch_reverse_det_<case>.json`.
Logs: `debug/exch_reverse_scratch/`.
