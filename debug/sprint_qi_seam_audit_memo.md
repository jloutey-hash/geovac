# Sprint memo — adversarial audit of the Q(i) / conductor-4 seam between Paper 56 (Dirac/Kramers) and Paper 59 (T2)

Date: 2026-08-21 | Branch: work/sparsity-boundary | **No papers modified** (audit only)
Discipline invoked: `memory/feedback_audit_numerical_claims.md` (this is an "X matches Y" claim);
`memory/feedback_verify_current_state.md` (owning sections read directly, not from memos).

Sources read: P56 `prop:hodge_cm_point` / `rem:cm_explains_hb` / `rem:qi_triangulation` /
`rem:paper59_cm` (lines ~1670-1760); P59 `sec:modular` (lines 667-895); P28 `thm:chi4` +
`rem:chi4_motivic` (lines 450-560); `debug/sprint_T2_zhou_bd_adaptation_memo.md`;
`debug/_kms_mt_torus.py` (re-run: 12/12 PASS); `tests/test_paper59_coarea_reduction.py`;
`debug/sprint_hodge_sl2_audit_qi_triangulation.py` (re-run).

---

## 1. VERDICT — **CONVERGENT**, not shared mechanism

Two structurally independent routes land on the same *arithmetic invariant* — a Q-rational complex
structure J with J^2 = -I, i.e. a Q(i)-structure on a 2-dimensional Q-vector space, with the
associated conductor-4 / mu_4 data — but **on two different underlying Q^2's, with no comparison map
exhibited between them**, and living in **two different Tannakian categories**:

| | route (a) Dirac/Kramers | route (b) T2 / Legendre |
|:--|:--|:--|
| carrier of the Q(i)-structure | `V_fund = Q^2` = a j=1/2 Dirac doublet | `H_1(E_lambda, Q)` of the Legendre fibre |
| origin of the "4" | half-integer shift (3/2) x parity-2 split => **quarter-integer** Hurwitz args | 2-torsion (lambda) x its **square root** => level Gamma(2) -> Gamma(4) |
| period category | mixed **Tate** over Q(i), `MT(Z[i,1/2])` (P28 `rem:chi4_motivic`) | mixed **elliptic** + **irregular** (Gevrey-1, Borel radius 2) over X(2) |
| conductor-4 content | **measured** (beta(2), beta(4) in `D_even - D_odd`, PSLQ 80-dig) | **not measured in T2**; measured only in the *family* (see 3.3) |

The fields coincide; the categories do not; the objects are not connected. That is the definition of
CONVERGENT. Two of the four candidate mechanisms survive as *real structure* (4.1, 4.2), but each
strengthens **one side only** — neither supplies a common cause.

**However**: the audit found a genuinely promotable structural upgrade on route (b) (3.3) that is
*stronger* than what P56 `rem:paper59_cm` currently claims, and **one MATERIAL defect** in P59's
load-bearing "requires G" argument (5.1). Both are stated below; neither applied.

---

## 2. Route (a): where chi_-4 is actually forced (P28)

`thm:chi4`: `D_even(s) - D_odd(s) = 2^{s-1}(beta(s) - beta(s-2))`, beta(s)=L(s,chi_-4). Traced to its
root, the forcing is a **two-factor 2x2**, and both factors are checkable:

- Dirac on S^3: `m = n + 3/2` (**half-integer** — the spin double cover). Splitting on n-parity gives
  Hurwitz arguments `3/4` and `5/4` -> **quarter**-integer shifts -> `zeta(s,3/4) - zeta(s,1/4) = -4^s beta(s)`.
  Verified here: `zeta(4,3/4) - zeta(4,1/4) = -253.16980524572296604...` = `-4^4 beta(4)` to 30 digits.
- The scalar null is *structural, not just empirical*: with the integer shift `m = n+1`, parity
  splitting can only produce **half**-integer Hurwitz arguments {1/2, 1}, whose values are
  `zeta(s,1/2) = (2^s - 1)zeta(s)` and `zeta(s)` — level <= 2, pure Tate over Q. Checked numerically
  (odd-m sum = `(1 - 2^-s)zeta(s)`).

So: **conductor 4 = (half-integer index from spin) x (a Z_2 parity split).** This is a genuine
forcing, and the P56 `rem:qi_triangulation` "spin-gated" reading is *correct for route (a)*.

*Audit note on the backing script.* `debug/sprint_hodge_sl2_audit_qi_triangulation.py` leg C compares
beta(2) against `zeta(2) = pi^2/6` — a **different observable**, not the scalar analogue of the same
even/odd parity computation. The conclusion survives (the half-vs-quarter shift argument above is the
real proof), but the script's leg C is weaker than the claim it is cited for.

---

## 3. Route (b): the critical tension, resolved

### 3.1 T2 is spinless — so the spin gate cannot apply, and it does not have to

T2 is a scalar Coulomb three-centre ERI. There is no half-integer index anywhere in it. The
`rem:qi_triangulation` reading ("both sides gated by half-integer spin") therefore **does not extend
to route (b)** — and P56 `rem:paper59_cm` does not claim it does; it claims only that the two routes
"reach the same CM field Q(i) from independent directions". That claim is literally true. The
question is whether it is *informative*.

### 3.2 The tau=i fibre hit is close to vacuous — P56 carries the corroboration without P59's own deflation

- `lambda = 1 - rho` sweeps (0,1), so tau sweeps the whole imaginary axis, which contains **every**
  imaginary-quadratic CM point `tau = i*sqrt(n)`. Passing "through every CM fibre" is a property of any
  family sweeping the full real modulus locus — it carries no arithmetic information.
- P59 itself states the deflation: *"landing on the Legendre family is expected for a four-branch-point
  curve, so the content here is the rational modulus map and the in-domain complex-multiplication
  fibres, not an exotic level."* **P56 `rem:paper59_cm` does not carry this sentence.**
- Worse for the "central fibre" framing: `rho = 1/2` is **not** physically distinguished. The physically
  distinguished point is `rho = 1` (coincident Fock scales) — which is the **cusp** `lambda=0`, the
  genus-0 degeneration, and the fixed point of the *physical* symmetry `Phi(rho) = Phi(1/rho)/rho^2`
  (read from `tests/test_paper59_coarea_reduction.py`). `rho = 1/2` is an arbitrary 2:1 scale ratio;
  it is the midpoint of the lambda-range, nothing more.
- Also: `tau = i` is **not** an elliptic point of Gamma(2) (Gamma(2) is torsion-free). The order-4
  stabilizer `S = [[0,-1],[1,0]]` lies in SL_2(Z) but not in Gamma(2); on X(2) it descends only to the
  order-2 Fricke involution `lambda <-> 1-lambda`. The full mu_4 is visible **only on H_1 of the CM
  fibre**, where S acts as multiplication by i on the lattice Z[i].

**Deflation verdict:** as written, `rem:paper59_cm`'s corroboration rests on the weakest available
leg. A *stronger* and fibre-independent leg exists, and the audit found it:

### 3.3 The real conductor-4 content of route (b): K = (pi/2) x (theta series of Z[i])   [PROMOTABLE]

P59 already records `K(m) = (pi/2) theta_3(0,q)^2` (validated 40 digits). The audit adds the arithmetic
identification the paper does not make:

> **theta_3(0,q)^2 is the theta series of the Gaussian integers.**
> `theta_3^2 = sum_{m,n in Z} q^{m^2+n^2} = sum_n r_2(n) q^n` with, by Jacobi's two-squares theorem,
> `r_2(n) = 4 sum_{d|n} chi_-4(d)`. So theta_3^2 **is** the weight-1, level-4 Eisenstein series with
> character chi_-4 — the CM/Hecke theta series of Q(i) — and its Dirichlet series is
> `sum r_2(n) n^-s = 4 zeta(s) beta(s)`.

Verified here to machine-exact: `theta_3(0,q)^2 - (1 + 4 sum_n (sum_{d|n} chi_-4(d)) q^n) = 0.0`
(220 terms, dps 30); and `K(lambda) - (pi/2) theta_3^2 = 0.0` at a generic tau.

Consequences (all favour route (b) being conductor-4 by construction, **not** by fibre accident):

1. The conductor-4 datum sits in the **fibre period at every fibre**, as a property of the Legendre
   family, not at the single point tau=i.
2. It explains *why* beta-values are the natural L-values of this family (`L(theta_3^2, s) = 4 zeta(s) beta(s)`)
   — i.e. it supplies the missing justification for admitting `G = beta(2)` into the PSLQ ring, in
   place of the defective one P59 currently gives (5.1).
3. **Second, independent level-4 witness in the physics.** The production twist's decay rate is
   `A = sqrt(c1) + sqrt(c2) = sqrt(u)(1 + sqrt(rho))` (read from `_Jsum`/`_four_bs` in the backing
   test), and the Stokes action is `sqrt(z*) = b + i sqrt(c1)`. So **sqrt(rho) = sqrt(1-lambda)** is a
   physical variable. But `sqrt(1-lambda) = theta_4^2/theta_3^2` (verified, 30 digits) is a
   **Gamma(4)-level** modular function, whereas lambda is Gamma(2). The Bessel twist therefore
   *promotes the level from 2 to 4* — and 4 = |disc Q(i)| = cond(chi_-4) = the level-4 field Q(mu_4)
   of Deligne 2010 cited in P56.

This is the audit's one net-positive finding, and it is exactly candidate mechanism **2**.

---

## 4. The four candidate mechanisms, adjudicated

### 4.1 Candidate 1 — half-integer modular weight <-> half-integer spin: **PARTLY REAL, but not a shared cause**

*Real content.* On the modular side the analogy is not a pun: half-integer weight forces level 4.
theta_3 lives on Gamma_0(4) with the theta multiplier `(c/d) eps_d^-1`, `eps_d = 1` or `i` according to
`d mod 4` — a mu_4-valued, conductor-4 cocycle. And theta_3^2 is exactly the chi_-4 weight-1 form (3.3).
So "half-integer index => 4 appears" holds on **both** sides:
route (a) half-integer *spin* -> half-integer Hurwitz shift -> quarter shifts -> chi_-4;
route (b) half-integer *weight* -> metaplectic/theta multiplier -> level 4 -> chi_-4.

*Against.* The two double covers are covers of **different groups** — `Mp_2(R) -> SL_2(R)`
(metaplectic, symplectic/Weil) versus `SU(2) -> SO(3)` (spin, compact). No homomorphism between them is
exhibited, and none is natural here. What is genuinely shared is the *arithmetic pattern* 4 = 2 x 2: a
Z_2-graded (double-cover) datum combined with a second Z_2 (parity split / square root of the modulus)
yields mu_4 as the relevant root-of-unity group, hence Q(i) as the field of definition.

**Verdict: shared pattern, not shared mechanism.** More than a pun (both 2's are forced, and forced in
the same "4 = 2x2" way), less than a cause (the 2's are structurally independent).

### 4.2 Candidate 2 — level tower Gamma(2) > Gamma(4), lambda <-> 2-torsion, Q(mu_4): **STRONGEST CANDIDATE; real on route (b), one-sided**

Substantiated in 3.3(3): the physical Bessel twist makes `sqrt(rho) = sqrt(1-lambda) = theta_4^2/theta_3^2`
a physical variable, which is Gamma(4)-level. Combined with `theta_3^2 = E_1(chi_-4)` this is a genuine,
checkable statement that T2's arithmetic home is level 4 = Q(mu_4) = Q(i), *independently of any CM fibre*.

**But it is a route-(b) statement only.** Route (a)'s level-4 arrives via Deligne–Glanois descent on
mixed-Tate periods (`MT(Z[i,1/2])`, P28 `rem:chi4_motivic`); route (b)'s arrives via a level structure
on an elliptic family. Both name the same field; neither is derived from the other, and route (b)'s
object is explicitly **not** mixed-Tate (P59: genuinely elliptic, and irregular/resurgent).

**Verdict: real mechanism on route (b); does not bridge.**

### 4.3 Candidate 3 — mu_4 torsion of the KMS/BW circle (this session's `_kms_mt_torus.py`): **RESTATES ROUTE (a); cannot be the common cause**

Re-ran: 12/12 exact. The findings are correct and non-trivial *as a cross-thread identification*: the
WH7 modular generator `K = diag(2 m_j)`, transported to the Kramers frame, is literally P56's Hodge
circle `cos t I + sin t J`; J is the beta/4 quarter-period point; the four quarter-points are
{I, J, -I, -J} = mu_4; and `e^{i pi K} = -I` on the spinor doublet versus `+I` on the scalar
(`K_scalar = diag(2,0,-2)`) — the spin double cover made visible **inside the thermal circle at beta/2**.

*Deflation the PI already named, and it bites:* any two maximal tori of SU(2) are conjugate, so the
existence of *some* frame in which `e^{itK}` is the Hodge circle is automatic once spec K = {+1,-1} and
J^2 = -I. The content is the identification of the corpus's *specific* operators, not the conjugacy.

*The decisive point against candidate 3 as a common cause:* **T2 has no KMS flow, no time, no thermal
circle, no beta.** Nothing in route (b) is a temporal compactification. Candidate 3 illuminates route
(a) (it re-derives the spin gate as a mu_4-vs-mu_2 torsion statement on the flow) and connects WH7 to
`prop:hodge_cm_point`, which is a real internal seam — but it is silent on route (b).

**Verdict: genuine intra-corpus seam (WH7 <-> P56); zero explanatory power for T2.**

### 4.4 Candidate 4 — real forms: compact SO(2)/CM/Q(i) vs split SO(1,1)/pure-Tate/Q: **PUN as stated; falsified by the corpus's own scalar sector**

*Seductive part.* The Q-group `SO(2) = {a^2+b^2=1}` is the norm-1 torus of Q(i); its character module is
Z with Galois acting by -1, whose Artin L-function is `L(s, chi_-4) = beta(s)`. So "compact Q-torus =>
conductor 4" looks forced.

*Why it fails.* The splitting field of the torus is **not** the field of the periods. The corpus's own
scalar sector is the counterexample: `K_scalar = diag(2,0,-2)` also generates a compact circle
(`e^{2 pi i K} = I`), i.e. also a Q-form of SO(2) — yet the scalar S^3 periods are pure Tate over Q
(conductor 1). Compactness gives you a torsion structure; it does not give you Q(i) periods.

*What actually discriminates* is whether the **mu_4 torsion point acts non-trivially** — i.e. whether
the representation is odd (half-integer weight): `e^{i pi K} = -I` (spinor, faithful through mu_4) vs
`+I` (scalar, factors through mu_2). That is candidate 1 again, not candidate 4.

Additionally, "Wick rotation = base change to Q(i)" is at odds with the corpus's own settled result:
WH7's Lorentzian closure (2026-06-19) is that the truncated BW boost is **compact** (integer spectrum,
`e^{2 pi i K} = I` bit-exact), so at finite cutoff there is no split torus to be the Lorentzian
partner — the split/non-compact side is strictly a continuum limit. The proposed dichotomy has no
finite-cutoff realisation.

**Verdict: pun. Do not promote.**

---

## 5. MATERIAL findings (not applied — papers untouched)

### 5.1 P59 `sec:modular`: the justification for admitting G into the corrected PSLQ ring is stated in the **wrong convention**, and in the paper's own convention both cited integrals are rational

P59 writes (justifying the corrected ring `{pi, K(1/2), 1/K(1/2), G}`):

> "...the weight-two Eisenstein L-value `G = beta(2) = L(2, chi_-4)` (Catalan's constant; native to
> precisely this class, `int_0^1 K(k) dk = 2G` and `int_0^1 E(k) dk = G + 1/2`)..."

Those two identities are in the **modulus** convention (`K(k)`, `m = k^2`). P59's own convention is the
**parameter** convention — its `K(1/2) = Gamma(1/4)^2/(4 sqrt(pi))` is `K(m=1/2)` (confirmed:
`mp.ellipk(0.5)` matches that Gamma-value to 40 digits). In the paper's own convention:

| | parameter measure dm (= the paper's K, and its physical d(rho)) | modulus measure dk (as quoted) |
|:--|:--|:--|
| `int_0^1 K` | **2** (exact, rational) | `2G` |
| `int_0^1 E` | **4/3** (exact, rational) | `G + 1/2` |

(All four verified to 40 digits; the first is also the paper's own "int_0^1 K = 2" measure validation
recorded in the v4.100.0 CHANGELOG / CLAUDE.md §2 bullet.) And the physical measure of the co-area
reduction is `d(rho) = -d(lambda)` — the **parameter** measure (read from `_Phi` in
`tests/test_paper59_coarea_reduction.py`).

=> **In the paper's own convention and its own physical measure, the two identities cited to make G
"native to precisely this class" produce no Catalan at all.** The G-in-the-ring decision may well be
right — 3.3 gives a *better*, convention-free reason (`theta_3^2 = E_1(chi_-4)`,
`L(theta_3^2, s) = 4 zeta(s) beta(s)`; plus the Gamma(4)-level sqrt(rho)) — but the reason currently
printed does not survive the audit. This is load-bearing: the whole "weight-three and requires G"
conclusion is downstream of it.

### 5.2 The "requires G" result is **conditional and under-powered**, not proven

The task framing calls it "the proven ... weight-3 and REQUIRES G result". The paper is more careful,
and correctly so:

- The conclusion is explicitly conditional: *"the finite closed form, **should it exist**, is
  weight-three and requires G."*
- It is about the **D->0 shadow**, not T2: *"and it is then the closed form of the D->0 shadow, the
  residual question being whether the physical D=1 value itself collapses onto that ring element
  (undecided at ~19 digits, awaiting ~32)."*
- Its second leg is **the same regime the paper elsewhere calls under-powered.** The conclusion rests
  on (i) a decisive weight<=2 negative in the corrected ring, plus (ii) a *period-only weight<=3*
  negative. But two paragraphs earlier: *"The wider searches, through weight three ... are only
  consistent with the negative at present: those bases are over-determined at ~19 digits — the
  spurious low-height hits differ from one precision to the next, and are matched by the decoy."*
  Leg (ii) is that search.

=> Any downstream use of "T2 requires Catalan" as evidence for the seam is **not supported**. There is
at present **no measured conductor-4 content in T2 itself** — only in the Legendre family it reduces to.

### 5.3 P56 `rem:paper59_cm` states the corroboration without P59's own deflation

`rem:paper59_cm` presents the tau=i hit as corroboration of `rem:cm_explains_hb`. P59 `sec:modular`
carries the deflation ("expected for a four-branch-point curve"; "not an exotic level"); P56 does not.
Per 3.2 the tau=i leg is close to vacuous, and per 3.3 a stronger fibre-independent leg is available.
Recommended (PI-gated, **not applied**): replace the tau=i leg with the theta_3^2 = Z[i]-theta-series
leg, and import P59's deflation sentence. Also: "the central fibre rho=1/2" invites a physical reading
it does not have — the physically central fibre is rho=1 (the cusp).

---

## 6. What would upgrade CONVERGENT -> SHARED MECHANISM, and the test

The seam is currently: *two isomorphic Q(i)-structures on two unrelated Q^2's.* Any two such are
abstractly isomorphic, so isomorphism alone is not evidence. An upgrade requires **a map**, not a
coincidence of invariants. Two concrete, falsifiable targets:

**T-1 (comparison map).** Exhibit a Q-linear, J-equivariant comparison between `(V_fund, J, Q)` — the
Dirac j=1/2 doublet with the Kramers structure — and `H_1(E_i, Q)` with its Z[i] action, that is
*natural in the GeoVac data* (i.e. constructed from the Fock projection / the momentum-space reduction,
not chosen). Falsifier: no such natural map exists, or the only one available is the tautological
"both are Q(i)-lines" identification. **Prediction if the seam is real:** the map should carry the
polarisation `QJ = I` to the elliptic curve's Riemann form; that is a checkable normalisation, and it
is exactly what a tautological identification would fail to fix.

**T-2 (period-level test — the decidable one).** The seam predicts T2's finite closed form contains
beta(2). This is currently undecided at ~19 digits. The paper's own gate is ~32 digits for the D=1
question and ~40 for the multi-fibre decision. Running the guarded, decoy-controlled PSLQ against the
corrected ring at >=32 digits **decides it**: a genuine G in T2 (route (b) conductor-4 realised in the
physical observable, not just in the family) is the strongest evidence the seam is more than
convergence; a decisive negative at 32+ digits *falsifies the period-level seam entirely* and leaves
`rem:paper59_cm` resting on the tau=i leg alone, which 3.2 shows is vacuous. Blocked by the same
Borel–Lambert second-cusp resummation already named as the frontier — i.e. **the seam question and the
T2 closed-form question are the same question**, which is itself a useful finding: the seam is not
independently testable at the period level until T2 is.

**Promotable now (no new work needed), P59 `sec:modular`, [MEASURED]:**

> `K(lambda) = (pi/2) theta_3(0,q)^2`, and `theta_3(0,q)^2` is the theta series of Z[i] — the weight-1
> level-4 Eisenstein series with character chi_-4 (Jacobi: `r_2(n) = 4 sum_{d|n} chi_-4(d)`), with
> `L(theta_3^2, s) = 4 zeta(s) beta(s)`. The conductor-4 arithmetic of this family is therefore carried
> by the fibre period at *every* fibre, not by the tau=i CM point; and the physical Bessel twist, whose
> decay rate `A = sqrt(c1)+sqrt(c2)` makes `sqrt(rho) = sqrt(1-lambda) = theta_4^2/theta_3^2` a physical
> variable, promotes the level Gamma(2) -> Gamma(4) = Q(mu_4) = Q(i).

Backing test (all four legs verified in this audit, exact or >=30 digits; would need writing as
`tests/test_paper59_theta_chi4.py`): (i) `theta_3^2 - (1 + 4 sum r_2 q^n) = 0`;
(ii) `K(lambda) = (pi/2) theta_3^2`; (iii) `sqrt(1-lambda) = theta_4^2/theta_3^2`;
(iv) the production evaluator's `A = sqrt(u)(1+sqrt(rho))` (structural read of `_Jsum`).

---

## 7. One-line summary

**CONVERGENT.** Route (a)'s conductor 4 is spin-forced (half-integer Dirac shift x parity split ->
quarter Hurwitz shifts -> chi_-4) and mixed-Tate; route (b)'s is level-forced (theta_3^2 = the Z[i]
theta series; the Bessel twist's sqrt(rho) promotes Gamma(2) -> Gamma(4)) and mixed-elliptic/irregular.
Same field, same 4 = 2x2 pattern, isomorphic Q(i)-structures — different carriers, different Tannakian
categories, no map. The tau=i fibre hit that P56 currently cites is close to vacuous; the KMS/mu_4
result is a real WH7 <-> P56 seam but restates route (a) only; the real-forms framing is falsified by
the corpus's own compact-but-pure-Tate scalar sector. Two MATERIAL items flagged for the PI: P59's
G-in-the-ring justification is stated in the wrong (modulus) convention — in the paper's own parameter
convention and physical d(rho) measure both cited integrals are rational (2 and 4/3) — and the
"requires G" conclusion is conditional, about the D->0 shadow, and rests partly on a PSLQ leg the paper
itself calls over-determined at ~19 digits.

---

## ADDENDUM (2026-08-21, v4.103.0 — post-dating the memo; see debug/sprint_aha_cross_corpus_memo.md)

Two corrections to §6 from the T2a involution track (driver `debug/aha_t2a_involution.py`, 44/44):

1. **T-1's polarisation discriminator is VOID.** h(Q(i)) = 1, so the principally polarised rank-2
   Q(i)-CM Hodge structure is unique up to isomorphism (unique up to the norm-1 torus = the Hodge
   circle itself when the polarisation is preserved) — EVERY identification, tautological ones
   included, carries QJ = I to the Riemann form. T-1 is restated (P56 rem:paper59_cm, v4.103.0):
   only NATURALITY/functoriality discriminates, e.g. presenting V_fund as H_1 of a
   GeoVac-constructed abelian variety.
2. **The τ=i deflation of §3.2 is re-derived from the group side and sharpened.** The only physical
   involution on the modulus is ρ↦1/ρ (the s↔t exchange; λ↦λ/(λ−1), coset T·Γ(2)); ρ=1/2 is half
   of the orbit {1/2, 2}, and the involution's own Z[i] fixed point (λ=2, ρ=−1, j=1728) is off the
   physical contour. Positive salvage: ⟨Γ(2), M⟩ = Γ_0(2) with M²=−I ≡ T mod 2, and the co-area
   fold IS the involution statement, so **T2's modular home is X_0(2)** (Paper 59 sec:modular,
   `tests/test_paper59_gamma0_2.py`).

Base-rate companion (T2b census, `debug/data/aha_period_census.md`): blind searches 0/30 on
disc −4; forced non-pure-Tate constructions 3/3. CONVERGENT stands, quantified: Q(i) is the first
CM field the substrate's rational arithmetic reaches (order-2 → order-4 lift: spin, metaplectic,
or modular torsion — conductor-4 periods; CM Hodge structure additionally needs the ±-paired
Q-rational carrier, T4b / Paper 24 Layer 7).

## ADDENDUM 2 (2026-08-21, v4.104.0): T-2 has FIRED — NEGATIVE

The pre-registered T-2 test (measured beta(2) in the physical integrated observable) has been run.
A new exact factorization (Paper 59 eq:kw: T2 = (8/pi) int dk int dw cos(kw) R(k,w)^2 — the (s,t)
double integral is the square of a 1-D integral) certified T2 to 66 cross-validated digits, and the
guarded, decoy-calibrated PSLQ at 64 working digits is DECISIVE-NEGATIVE for the corrected wt<=3
ring {pi, K(1/2)^{+-1}, G} at height <= 10, for disc-4 wt<=3 and disc-8 wt<=2 at heights up to
1e12, and for eight targeted Catalan probes. NO measured conductor-4 content exists in the
integrated observable at clean heights: the Q(i) arithmetic is a property of the FAMILY (theta_3^2)
and of the D->0 shadow, not of the physical value. Riders: height budgets scale 10^(D/n) in ring
dimension n (this memo's "decidable at >= 32 digits" was never true for the dim-20 ring); the
negative excludes clean closed forms, not large-height ones. With T-1 tautology-blocked (Addendum 1)
and T-2 negative, the seam stands CONVERGENT on all tested legs. See
debug/sprint_beta2_exchange_davidson_memo.md + debug/beta2_track_a_findings.md.
