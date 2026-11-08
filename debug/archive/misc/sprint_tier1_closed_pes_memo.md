# Sprint: Tier-1 zero-quadrature extensions — the closed-form H₂ PES + certified table growth + u2 — 2026-08-21

**Origin:** PI approved the Tier-1 plan ("anything else zero-quadrature?"). Main-session centerpiece +
two agents (D1 helium, D2 diatomic table) + the u2 compute landing. Canonical memo; details:
`debug/data/h2_closed_pes.json`, `debug/qfd_he_findings.md`, `debug/qfd_table_findings.md`.

## The centerpiece — the H₂ potential curve as ONE formula (main session)

In the minimal 1s/1s basis the singlet-Σg⁺ CI space is {|g²⟩,|u²⟩}, so the FCI ground state is the
quadratic-formula root of a 2×2 with closed-form entries ⇒ **E(R) is a single symbolic expression,
atoms exactly {exp, E₁, log, γ_E}** (`geovac.qfd_assemble.h2_closed_form_E`; verified vs assembled
FCI to 1e-42 at R=1.4 and 2.0). Its exact derivatives certify the equilibrium constants by
closed-form Newton (dps 45 vs 60 agree 9.5e-46):
- R_eq = 1.667999966972748724927046410307264513349804147483 bohr
- D_e  = 0.11865036209829531185349776433388467169190968328627 Ha (dissociation to exactly −1 at 3e-49
  ⇒ D_e itself closed-form)
- k    = 0.25470393430796998115272793572435567267825178990899 Ha/bohr²
Honest frame: minimal-basis model constants (exact R_eq = 1.401; ω_e conversion = Layer-2). Captured:
Paper 58 sec:qfd eq:h2_pes; tests test_paper58_qfd.py::test_closed_form_pes_* (slow-marked); the
PES functions promoted into geovac/qfd_assemble.py.

## D1 — certified helium graph-native CI (GO)

n_max=1 EXACT: E = −11/4 Ha. n_max=2..5 at 35 digits each with a RIGOROUS residual bound
(‖Hv−λv‖ ≤ 1e-75 + rounding perturbation carried to the exact matrix), second eigensolver to 75
digits, float pipeline 1e-16-consistent; all 145 hard-coded rationals in casimir_ci re-derived, 0
mismatches. Honest-scope caveat + Pekeris/Drake % (5.29→0.250) test-enforced per entry.
`benchmarks/certified_reference/entries_helium.py`. Perf: float64-preconditioned residual refinement
(O(n²)/iter) replaced the O(n³) LU high-precision eigensolver.

## D2 — certified minimal-diatomic table (GO)

| system | R | Nₑ | τ series | digits | E_total (Ha) |
|:--|--:|--:|:--|--:|:--|
| H₂⁺ | 2.0 | 1 | terminates | 40 | −0.5537714953184827365067633613319649091848 |
| H₂⁺ | 1.4 | 1 | terminates | 40 | −0.4713457017330697952317131946245614997863 |
| HeH⁺ | 1.46 | 2 | ∞, τ_max16 | 40 | −2.895950902302325174310463610145323586925 |
| He₂²⁺ | 1.3 | 2 | terminates | 40 | −3.586503080128554308997786255254462141274 |
| BeH⁺ | 2.5 | 4 | ∞, τ 22/18/16 | 31 | −14.69882770595311489573042272632 |

Homonuclear termination verified as symbolic zeros through τ=8; heteronuclear tails bounded with
per-system amplification (HeH⁺ 41 digits not binding; BeH⁺ 1.31e-30 Ha = the binding constraint).
Validations 10–13 orders inside gates. Findings: (a) debug/qfd_lih.py's exchange-vs-quadrature line
was mislabelled (compared full closed form vs τ≤4 reference — the residual was omitted τ terms, not
quadrature error; label FIXED; LiH's 30-digit claim unaffected — the τ-bound nets it); (b) the l>0
hybrid Z_B<Z_A gap is UNREACHABLE from s-only bases (equal-rate (2,2) works via the direct route).

## u2 + PSLQ (main session)

u2 (dps 150, K=560, disjoint config) assembled; **u1 and u2 agree to 83 significant digits**
(6.9e-85) ⇒ cross-validated T2 = 0.395355765901713964325229296804847564260563977867082108935234265469515508685067297425 (83 digits; stored debug/data/t2_83digit.txt). Guarded PSLQ re-run at 83 digits
(driver's height formula confirmed corrected to 10^(D/n); stale docstring fixed): the pre-registered
T-2 ring DECISIVE-NEG extends h≤10 → **h≤10²**; wt≤2 rings → h≤10⁶–10⁷; disc-4 wt≤1 / probes →
h≤10¹²; +ln2 ring still UNDERPOWERED. The h≤10⁴ goal on the dim-20 ring needs ~120 cross-validated
digits (a u3 at higher node counts — optional follow-on).

## Production fixes (fix-on-sight)

- `geovac/casimir_ci.py::_wigner3j`: phase exponents `(-1)**(negative int)` returned float ±1.0,
  silently poisoning exact-Fraction reuse (D1's catch) — both phase sites now use `% 2`. 215
  adjacent tests green.
- `debug/qfd_lih.py` diagnostic label corrected (above).

## Artifact state

**64 entries / 7 categories**, `--check` 0 mismatches, 12 tests green. New since v4.105.0: 5 helium
(incl. the exact −11/4), 5 diatomics, 3 H₂-PES constants (digit-claim discipline enforced itself:
the schema test rejected 49-printed/45-claimed until truncated — the artifact self-polices).

## Follow-ons
- u3 for ~120 digits → h≤10⁴ on the pre-registered ring (optional; machinery ready).
- l>0 one-electron closed forms → s+p diatomics (the approved next build).
- Isoenergetic circuit-level reach; PSD-shift amplification pricing (standing queue).
