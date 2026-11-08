"""DoD criteria name the GATE THAT OWNS a number, never the number.

PI-approved 2026-09-11, after /qa paper_60 stopped at protocol step 1: the
pre-registered goalposts named RETIRED values as canonical, which inverts the
test -- W1 ranked HIGHEST asserted "Canonical = 0.84" for an exponent retired on
2026-09-07 in favour of 0.82, so reviewing against it would have marked the
CORRECTED paper MATERIAL and passed a paper that still carried the old value.

Same class the 2026-09-08 paper_61 run named (a .done.md ratifying a retired
claim).  Second instance => a pattern, and no gate covers it: C16 scopes to
papers and their code, not to criteria records.

The fix is structural rather than a repair.  A DoD freezes CRITERIA -- relations,
structures, prohibitions, which are stable.  MEASUREMENTS belong to the gates
that already own them (C21's numeric registry, C17's headline families), both of
which carry provenance, retirement tracking and a check.  A DoD literal has
neither, and CLAUDE.md Sec.13.11 rule 8 ("one canonical record per fact") already
forbids the duplication: that exponent lived in the paper, the registry, a C17
family AND the DoD -- and the copy that rotted was the only one with no gate.

Delegation is also STRONGER: "abstract == body == registry key" checks three loci
agree where the literal checked two.

Historical/narrative sections keep their dated numbers -- those are chronicle,
not goalposts, and are correct as records of what was true when written.
"""
import io

CRIT = "docs/qa/criteria.md"
DOD = "docs/qa/paper_60.done.md"

# --------------------------------------------------------------- convention
c = io.open(CRIT, encoding="utf-8").read()
ANCHOR = "## Branch-specific criteria (C14+)"
CONVENTION = """## DoD criteria name the OWNING GATE, never the value (added 2026-09-11, PI direction)

**In a `.done.md`, criteria sections state relations and delegate every
load-bearing number to the gate that owns it. They never write the literal.**

- a C21 registry key (`p60_onenorm_exponent`), or
- a C17 headline family id (`paper60-molecular-lambda-exponent`), or
- the paper's own equation label (`eq:sublinear`) when the criterion is
  "these two loci agree".

So instead of *"Canonical = the labelled body equation, 0.84; the abstract must
match"*, write *"the abstract's K-exponent, `eq:sublinear`, and registry key
`p60_onenorm_exponent` must all agree; any disagreement = MATERIAL."*

**Why, measured.** `/qa paper_60` (2026-09-11) stopped at protocol step 1: the
frozen goalposts asserted a RETIRED value as canonical. W1 -- the watch-note
ranked HIGHEST -- required the abstract to match `K^0.84`, an exponent retired on
2026-09-07 as a truncated-radial-domain artifact and replaced by `K^0.82`. The
paper had been corrected; the criteria had not. **Reviewing against them would
have graded the correct paper MATERIAL and passed a wrong one** -- the gate
returning the wrong answer at full confidence. Its C8 block additionally froze
two claims retired or withdrawn since: "residual = basis incompleteness" (it is
the scale lock, `eq:scale_lock`) and a comparison to Avery's 102-configuration
figure (withdrawn -- different posing, `eq:no_selection`).

This is the second instance of the class the 2026-09-08 `paper_61` run named --
a `.done.md` ratifying a retired claim, a re-infection vector because a criteria
file is the goalpost the next certifying run measures against. **No deterministic
gate covers it:** C16 scans papers and their contained code/test modules, not
`docs/qa/*.done.md`.

**Two objections, answered.**

*Doesn't delegation reopen goalpost-moving, since the registry can change?* No.
A registry change is **disciplined and gated** -- C21 blocks retired values and
names the key replacing them, Sec.15 rule 3 forbids registering anything
unmeasured, Sec.15 rule 2 forces a prose re-read at every locus the gate names.
A DoD literal is **undisciplined and frozen**. Delegation moves the number to the
only place equipped to keep it honest. The freeze that matters -- criteria do not
move *during* a run -- is untouched; pin the registry commit at run start if you
want belt-and-braces.

*Doesn't this weaken pre-registration?* It strengthens it. The criterion becomes
a relation over three loci instead of a literal over two, and a relation cannot
go stale.

**Scope.** This governs *criteria* sections only -- branch-defining criteria,
watch-notes, the C8 headline enumeration. **Historical and narrative sections
keep their dated numbers**: a run's FAIL account and its measurement tables are
chronicle, correct as records of what was true when written, and rewriting them
would destroy the audit trail. The distinction is whether a number is asserted as
*canonical for the next run* (delegate) or *reported as what was measured then*
(leave).

"""
assert ANCHOR in c and "OWNING GATE" not in c
c = c.replace(ANCHOR, CONVENTION + ANCHOR, 1)
io.open(CRIT, "w", encoding="utf-8").write(c)
print("criteria.md: delegation convention added")

# ------------------------------------------------------------ paper_60 DoD
d = io.open(DOD, encoding="utf-8").read()

OLD_BRANCH = """1. **The atomic sublinearity is a CONFIGURATION-count (K) statement, not a qubit-count (Q)
   one.** ‖M‖1∼K^0.84 is sublinear in CI configuration count for a single atom; the paper
   itself flags it as not-yet-mapped to Q. No prose may imply a qubit-count sublinearity or
   a many-electron win. **exact ≠ accurate (inherited Paper 58 W1):** metric-free / pi-free /
   pure-number structure buys *encoding cost*, not accuracy.
2. **Molecules are polynomial, and the paper must say so.** The N-electron interacting
   molecular block-encoding 1-norm is ∼n_orb^2.2 (standard second-quantization ballpark, no
   advantage over DF/THC); the paper concedes "not a sublinear matrix." Any prose implying a
   molecular many-electron 1-norm advantage = MATERIAL."""

NEW_BRANCH = """1. **The atomic sublinearity is a CONFIGURATION-count (K) statement, not a qubit-count (Q)
   one.** The exponent of `eq:sublinear` is sublinear in CI *configuration* count for a
   single atom — value owned by C21 key `p60_onenorm_exponent`, and it is a WINDOW fit, not
   a regime. No prose may imply a qubit-count sublinearity or a many-electron win.
   **exact ≠ accurate (inherited Paper 58 W1):** metric-free / pi-free / pure-number
   structure buys *encoding cost*, not accuracy.
2. **Molecules are polynomial, and the paper must say so.** The N-electron interacting
   molecular block-encoding 1-norm is polynomial (standard second-quantization ballpark, no
   advantage over DF/THC) — value owned by C17 family `paper60-molecular-lambda-exponent`;
   the paper concedes "not a sublinear matrix." Any prose implying a molecular many-electron
   1-norm advantage = MATERIAL."""
assert OLD_BRANCH in d
d = d.replace(OLD_BRANCH, NEW_BRANCH, 1)

OLD_W1 = """- **W1 — abstract K-exponent drift [HIGHEST, headline-number].** Abstract (line 46)
  ‖M‖1∼K^0.78; body `eq:sublinear` (labelled, line 247) ‖M‖1∼K^0.84; CLAUDE §2 + the
  lit-comparison memo both use 0.84. **Canonical = the labelled body equation, 0.84.** The
  abstract figure must match. Mismatch = MATERIAL (C8/C17)."""
NEW_W1 = """- **W1 — K-exponent agreement [HIGHEST, headline-number].** The abstract's K-exponent,
  the labelled body equation `eq:sublinear`, and C21 key `p60_onenorm_exponent` must **all
  three agree**. Any disagreement = MATERIAL (C8/C17/C21). *No value is written here by
  design* — see "DoD criteria name the OWNING GATE" in `criteria.md`. This note previously
  froze `0.84` as canonical and was still asserting it after that value was retired
  (2026-09-07), which is what stopped the 2026-09-11 run at protocol step 1."""
assert OLD_W1 in d
d = d.replace(OLD_W1, NEW_W1, 1)

OLD_W6 = """- **W6 — He is DELIBERATELY low-accuracy.** −2.897 (spdf, K=164) vs exact −2.90372 (~7 mHa)
  is worse than STO-3G-class; the Goscinskian basis is deliberately poor for the He GS. The
  validation proves the *machinery is correct*, NOT that it is accurate. "Accurate helium"
  reading = MATERIAL."""
NEW_W6 = """- **W6 — He is DELIBERATELY low-accuracy, and the RESIDUAL IS THE SCALE LOCK.** The
  helium ladder is worse than STO-3G-class; the validation proves the *machinery is correct*,
  NOT that it is accurate. "Accurate helium" reading = MATERIAL. **Superseded 2026-09-08:**
  this note used to attribute the residual to *basis incompleteness*. It is not — a
  variational CI over the identical span reaches far closer (C21 key
  `p60_span_deficit_spdf`), so the residual is the scale lock, `eq:scale_lock`. Any prose
  still calling it basis incompleteness = MATERIAL."""
assert OLD_W6 in d
d = d.replace(OLD_W6, NEW_W6, 1)

OLD_C8_HEAD = """## C8 headlines (enumerated, with tiers — the frozen goalposts; canonical = BODY values)"""
NEW_C8_HEAD = """## C8 headlines (enumerated, with tiers — the frozen goalposts)

> **Values are NOT written here.** Each headline names the claim, its tier, and the
> **gate that owns its number** (C21 registry key / C17 family / equation label). See
> "DoD criteria name the OWNING GATE" in `criteria.md`. The pre-2026-09-11 version of this
> block froze three claims that were later retired or withdrawn — `K^0.84` as canonical,
> "residual = basis incompleteness", and the Avery 102-configuration comparison — and would
> have failed a corrected paper."""
assert OLD_C8_HEAD in d
d = d.replace(OLD_C8_HEAD, NEW_C8_HEAD, 1)

OLD_C8_3 = """3. **Helium validation [MEASURED].** Single Goscinskian config −2.847 Ha = textbook
   variational (−2.84766; the pure number 5/8·sqrt2^−1); bare (no V′) matrix −4.0 (two He+ 1s,
   exact non-interacting); multiconfig −2.847(1s^2)→−2.873(s)→−2.894(+p)→−2.897(spdf,K=164) →
   exact −2.90372; ~7 mHa residual = basis incompleteness (deliberately poor basis; Avery &
   Avery reach −2.90250 with 102 configs), every point above exact (no overshoot)."""
NEW_C8_3 = """3. **Helium validation [MEASURED].** Single Goscinskian config reproduces the textbook
   single-exponent variational value from the pure number `5/8·sqrt2^-1`; the bare (no V′)
   matrix returns two non-interacting He+ 1s exactly; the multiconfiguration ladder descends
   monotonically toward the exact non-relativistic energy with **every point above it**.
   **The residual is the SCALE LOCK, not basis incompleteness** (`eq:scale_lock`; span
   deficit owned by C21 `p60_span_deficit_spdf`). **The Avery 102-configuration comparison is
   WITHDRAWN** — `eq:no_selection` proves that figure unreachable in the locked posing at any
   K or selection, and a relayed consultation attributes it to a scale-optimized
   (ordinary variational CI) calculation. Any prose reinstating either = MATERIAL."""
assert OLD_C8_3 in d
d = d.replace(OLD_C8_3, NEW_C8_3, 1)

OLD_C8_4 = """4. **Atomic sublinear 1-norm [MEASURED].** ‖M‖1∼**K^0.84** (`eq:sublinear`, full s+p+d+f) in
   configuration count K — sublinear, opposite of the naive inflation. (CANONICAL 0.84;
   abstract 0.78 = the W1 drift.)"""
NEW_C8_4 = """4. **Atomic sublinear 1-norm [MEASURED].** `eq:sublinear` — ‖M‖₁ grows more slowly than
   the configuration count K on a **converged radial domain**, opposite of the naive
   inflation. Exponent owned by C21 `p60_onenorm_exponent` (s-only companion
   `p60_onenorm_exponent_sonly`); split owned by `p60_T0_asymptotic_exponent` /
   `p60_Tprime_full_exponent`. **It is a WINDOW fit, not an asymptotic regime** — the local
   slope falls monotonically and no window value is stable; prose asserting a regime, or any
   value measured on a truncated domain, = MATERIAL."""
assert OLD_C8_4 in d
d = d.replace(OLD_C8_4, NEW_C8_4, 1)

NEW_EQUATIONS = """
11. **Metric-free ⟺ the scale is locked to the eigenvalue [INTERNAL THEOREM].**
    `eq:W_diagonal` (the one-body Coulomb metric is exactly diagonal) and `eq:scale_lock`
    (metric-free ⟺ E = −λ²/2 ⟺ λ = p_κ). The posing IS the variational problem of its own
    span at the one scale where the L² metric cancels. Backing
    `tests/test_paper60_scale_lock.py`. **Consequences that must travel with it:** the
    variational bound is automatic, not fortunate; freeing the scale reaches chemical
    accuracy but hands back BOTH advertised cost risks and the encoding advantage
    (`p60_freescale_set_sonly`, a MATCHED SET on one ladder — quoting one member against a
    value from another ladder = MATERIAL); and the floor is a **ground-state pathology**
    (`p60_posing_cost_ground` vs `p60_posing_cost_exc`, `p60_gnd_ratio_k452` vs
    `p60_exc_ratio_k452`), not a property of the method.
12. **No selection rescues the locked posing [INTERNAL THEOREM].** `eq:no_selection` — T′
    being pure numbers makes any sub-family's secular matrix exactly a principal submatrix,
    so Cauchy interlacing gives E(A) ≥ E(M). Measured corollary owned by C21
    `p60_best102_locked`. Backing `tests/test_paper60_no_selection.py`. The property that
    makes the encoding attractive is what supplies the bound.
13. **The general-V₀ form, and what the molecular metric IS [INTERNAL THEOREM].**
    `eq:general_v0` — V·C = V₀·B·C, every L² overlap cancelling for ANY local V₀, orthonormal
    configurations or not; atomic specialisation reproduces `eq:secular`. Backing
    `tests/test_paper60_general_v0.py`. **Provenance split that must stay intact:** the
    Shibuya–Wulfman integrals are Avery's; the identification of that matrix *as* the
    V₀-weighted overlap is ours, and the variational bound for the fixed-scale metric-free
    problem is not in his canon. Crediting either to Avery = MATERIAL.
14. **The accuracy mechanism is AVERY'S, the price in qubits is OURS [ESTABLISHED, from
    Avery].** In-out radial correlation; split-shell 1s1s′ carries two independent exponents;
    a Goscinskian 1s² pins both electrons to one exponent for *any* weighting potential.
    Presenting this mechanism as a GeoVac discovery = MATERIAL.
"""
OLD_SEED_HEAD = "## Seeding plan (worktree only; never touches the real corpus)"
assert OLD_SEED_HEAD in d
d = d.replace(OLD_SEED_HEAD, NEW_EQUATIONS + "\n" + OLD_SEED_HEAD, 1)

io.open(DOD, "w", encoding="utf-8").write(d)
print("paper_60.done.md: branch criterion + W1 + W6 + C8 head + C8.3/C8.4 delegated;"
      " 4 new equations pre-registered as C8.11-14")
