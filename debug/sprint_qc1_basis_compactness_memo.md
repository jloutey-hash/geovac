# QC-1 — the basis-compactness claim, measured. Verdict: **LOSS**

**Date:** 2026-08-12 · branch `work/sparsity-boundary`
**Driver:** `debug/qc1_basis_compactness.py` (both legs; log is `*.log`-gitignored,
so this memo is the record)

## The claim, and why it mattered

Asserted in `memory/native_two_center_eri_engine.md` since 2026-08-09, never
measured:

> Slater functions have the correct nuclear cusp and correct exponential tail;
> Gaussians have neither. Fewer functions for equal accuracy ⇒ fewer spatial
> orbitals ⇒ **fewer qubits**, since qubit count scales with orbital count and
> integral evaluation is offline preprocessing.

This is the *entire* quantum-resource case for the native two-center ERI engine.
Sparsity was already known not to improve (N3b: cross-center kills Gaunt
l-selection, tensor inflates 13.8–15×). Compactness was the only remaining
mechanism.

## Leg 1 — the naive comparison, which looked like a decisive win

H₂ at R = 1.4 bohr, 2e FCI, s-only, identical pipeline
(`integral_set_md` → Löwdin → `fci_ground` → JW), both families variationally
optimised at every size.

| M | qubits | Pauli | Slater E | Gaussian E | gap |
|---|---|---|---|---|---|
| 2 | 4 | 26 | −1.147769 | −0.991219 | 0.157 |
| 4 | 8 | 360 | −1.152755 | −1.122064 | 0.031 |
| 6 | 12 | — | −1.155028 | −1.148588 | 0.0064 |

Pre-registered WIN condition (Slater at M beats Gaussian at M+2) satisfied at
**both** testable points. Read naively: Slater at 8 qubits / 360 Pauli beats
Gaussian at 12 qubits / 1818 Pauli.

## Leg 2 — the control that invalidates Leg 1

Leg 1's "Gaussian" family is **uncontracted single primitives**. That is not what
a Gaussian basis set is, and it is the wrong comparison for qubit cost:

- qubit count = 2M with M the number of **contracted** basis functions;
- contraction depth does not enter M, so it is **free** in qubit terms;
- a Slater function **is** a contracted Gaussian — that is literally how the
  script evaluates it, as a k-primitive fit of e^{−ζr}.

So the real question is: at fixed M, how much does contraction depth buy?
Held at M = 2 — one contracted function per H, hence **4 qubits and 26 Pauli
terms for every row**:

| k primitives | ⟨fit\|STO⟩ | E (Ha) | gap to k=10 |
|---|---|---|---|
| 1 | 0.97840439 | −0.991219 | +1.57e-01 |
| 2 | 0.99841970 | −1.114097 | +3.37e-02 |
| 3 | 0.99983474 | −1.138757 | +9.01e-03 |
| 4 | 0.99997812 | −1.145181 | +2.59e-03 |
| 6 | 0.99999938 | −1.147480 | +2.90e-04 |
| 10 | 1.00000000 | −1.147769 | 0 |

**Leg 1's entire M = 2 "Slater advantage" was 1.57e-01 Ha. Leg 2 reproduces
exactly that number by varying only the contraction depth, at fixed M, fixed
qubits, fixed Pauli count.** The k=1 row *is* Leg 1's Gaussian family; the k=10
row *is* Leg 1's Slater family. The gap was never a basis-family effect.

## Verdict: LOSS on the pre-registered gate

The gate read: *LOSS if Gaussian matches or beats Slater at matched M.* A
contracted Gaussian set matches Slater at matched M — it **is** the Slater
function — so the compactness argument yields **no qubit advantage and no Pauli
advantage**.

What Slater genuinely buys is accuracy **per primitive**, which is exactly the
axis that does not enter qubit cost. And primitive count is where Slater is
*worse* off classically: it needs a special integral engine, which is the thing
the last two days built.

**Consequence.** With sparsity already ruled out (N3b) and compactness now ruled
out, the native two-center ERI engine has **no established quantum-resource
advantage**. The memory's framing must move from "asserted, unquantified" to
**tested-negative**.

## Honest bounds on this result

- H₂, 2 electrons, s-only, minimal sizes. It does not establish scaling, and it
  does not exclude a residual advantage at high accuracy or for heavy atoms,
  where standard Gaussian contractions (optimised for atoms) may fit molecular
  environments less well.
- It *does* establish that the stated mechanism — cusp/tail ⇒ fewer contracted
  functions ⇒ fewer qubits — does not operate, because contraction absorbs it.
- The k=3 row is the STO-3G construction, and it already recovers 94% of the gap.
  Standard minimal Gaussian bases are built exactly this way, so real basis sets
  already capture most of what the claim attributed to Slater.

## Process note

Leg 1 passed its own pre-registered gate at both testable points and would have
been reported as a clean win. It was wrong because the *comparison* was unfair,
not because the measurement was wrong — the same failure mode as the increment-3b
quartet that returned zero and validated nothing. A gate can only test what the
comparison lets it see; pre-registering the threshold does not protect against
choosing the wrong control.
