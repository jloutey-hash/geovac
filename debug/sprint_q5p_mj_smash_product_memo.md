# Sprint Q5'-mJ-Smash-Product (T4) — bit-exact $\Delta m_J \in \{-1, 0, +1\}$ grading and $U(1)$ z-rotation eigenspace decomposition of the m_J-resolved OffDiag substrate at $n_{\max} = 2$

**Date:** 2026-06-06 (T4 of the Q5'-HardParts-Round3 sprint)
**Driver:** `debug/compute_q5p_mj_smash_product.py`
**Data:** `debug/data/sprint_q5p_mj_smash_product.json`
**Wall time:** 0.01 s
**Discipline:** bit-exact `sympy.Rational` throughout; no PSLQ; no floats.

---

## 1. TL;DR

**Verdict: POSITIVE-FOUNDATION.** The $m_J$-resolved OffDiag substrate at $n_{\max} = 2$ admits a **bit-exact $\mathbb{Z}_3$-grading by $\Delta m_J \in \{-1, 0, +1\}$** (Wigner-Eckart selection for the vector operator $A = \kappa A_{\mathrm{graph}}$) with **symmetric dimension count $30 / 40 / 30$**:

$$\mathcal{A}^{m_J\text{-OD}, (n_{\max}=2)} \;=\; \mathcal{A}_{-1} \oplus \mathcal{A}_0 \oplus \mathcal{A}_{+1}, \qquad \dim \mathcal{A}_{\pm 1} = 30, \qquad \dim \mathcal{A}_0 = 40, \qquad \mathrm{total} = 100.$$

The $U(1)$ z-rotation $U_\theta |n, \kappa, m_J\rangle = e^{i m_J \theta} |n, \kappa, m_J\rangle$ acts diagonally on the $m_J$-resolved transitions:
$$U_\theta \, T_{(n', \kappa', m_J') \to (n, \kappa, m_J)} \, U_\theta^{-1} \;=\; e^{i (\Delta m_J) \theta} \, T_{(n', \kappa', m_J') \to (n, \kappa, m_J)}.$$

This **exactly reproduces L5's three phase factors $\{1, e^{i\theta}, e^{-i\theta}\}$** (Sprint Q5'-Combined-Substrate, v3.63.0) — the $m_J$-resolved refinement DIAGONALIZES the $U(1)$ action that L5 identified as broken on the basic substrate.

**Critical scope finding:** the $j$-content at $n_{\max} = 2$ includes $j = 3/2$ states (from sector $(2, 1)$ with $\kappa = -2$ and sector $(2, 2)$ with $\kappa = +2$). Therefore the **L2 decorated-PW substrate at $j_{\max} = 1/2$ is INSUFFICIENT** for the full SU(2) action on the $m_J$-resolved OffDiag at $n_{\max} = 2$; **L2 must be extended to $j_{\max} \ge 3/2$** for the smash-product construction. This is named as a sprint-scale follow-on (~1 day).

---

## 2. Verdict against decision gate

| Gate | Selected? | Why |
|:-----|:---------:|:----|
| **POSITIVE-FOUNDATION** | **selected** | $\Delta m_J \in \{-1, 0, +1\}$ grading verified bit-exactly across all 100 non-zero off-diagonal Dirac entries; symmetric $\pm 1$ dimension count $30 = 30$; $U(1)$ z-rotation action explicitly diagonalised; L5's three-phase finding reproduced bit-exactly. Foundation for multi-year smash-product construction laid. |
| L5 already closed | not applicable | L5 ruled out smash-product at BASIC (non-$m_J$-resolved) substrate; T4 verifies the $m_J$-refined substrate IS the right load-bearing data. |
| STOP | not applicable | No structural obstruction at the $m_J$-resolved level. L2 extension is named as a sprint-scale follow-on, not a blocker. |

---

## 3. Bit-exact panel at $n_{\max} = 2$

### 3.1 Dirac state enumeration

The 16-dim Hilbert space at $n_{\max} = 2$ decomposes into 16 Dirac states $|n, \kappa, m_J\rangle$ across 5 sectors:

| Sector $(n, l)$ | dim | $j$-content | $m_J$ states |
|:----------------:|:---:|:----------:|:-----:|
| $(1, 0)$ | 2 | $j = 1/2$ ($\kappa = -1$) | $m_J \in \{-1/2, +1/2\}$ |
| $(1, 1)$ | 2 | $j = 1/2$ ($\kappa = +1$) | $m_J \in \{-1/2, +1/2\}$ |
| $(2, 0)$ | 2 | $j = 1/2$ ($\kappa = -1$) | $m_J \in \{-1/2, +1/2\}$ |
| $(2, 1)$ | 6 | $j = 3/2$ ($\kappa = -2$) + $j = 1/2$ ($\kappa = +1$) | $m_J \in \{-3/2, -1/2, +1/2, +3/2\}$ ($j=3/2$) + $\{-1/2, +1/2\}$ ($j=1/2$) |
| $(2, 2)$ | 4 | $j = 3/2$ ($\kappa = +2$) | $m_J \in \{-3/2, -1/2, +1/2, +3/2\}$ |

**Total dim H = 16 ✓**

### 3.2 Off-diagonal Dirac transitions

The off-diagonal Dirac $A = \kappa \cdot A_{\mathrm{graph}}$ has **100 non-zero entries** at $n_{\max} = 2$. Each entry corresponds to a transition $T_{(n', \kappa', m_J') \to (n, \kappa, m_J)}$ with bit-exact value $A[i, j] \in \kappa \cdot \mathbb{Q}$.

### 3.3 $\Delta m_J$ grading

Bit-exact panel: every non-zero entry has $2 \Delta m_J \in \{-2, 0, +2\}$, i.e., $\Delta m_J \in \{-1, 0, +1\}$. **Wigner-Eckart selection for vector operator passes 100/100 bit-exactly.**

| $\Delta m_J$ | # transitions | Phase under $U_\theta$ |
|:-----------:|:------:|:--------:|
| $-1$ | 30 | $e^{-i\theta}$ |
| $0$ | 40 | $1$ |
| $+1$ | 30 | $e^{+i\theta}$ |
| **Total** | **100** | (sum: $\{1, e^{i\theta}, e^{-i\theta}\}$ matches L5) |

### 3.4 Symmetric $\pm 1$ dimension count

The bit-exact dimensions $\dim \mathcal{A}_{-1} = \dim \mathcal{A}_{+1} = 30$ are forced by the $J$-reality (Sprint Q5'-Stage1-Prosystem §3.2: $J$ permutes $m_J \to -m_J$, sending transitions $T_{m_J \to m_J + 1} \mapsto T_{-m_J \to -m_J - 1}$).

The asymmetric $\dim \mathcal{A}_0 = 40$ vs $\dim \mathcal{A}_{\pm 1} = 30$ reflects the fact that $m_J$-preserving transitions are more numerous (every same-$m_J$ transition across sectors counts, vs $m_J$-changing transitions which require both $\Delta m_J = \pm 1$ AND the appropriate sector connectivity).

---

## 4. Structural readings

### 4.1 L5's three-phase finding is the $\Delta m_J$ grading

L5 (Sprint Q5'-Combined-Substrate, v3.63.0) found: the natural $SU(2)$ z-rotation $U_\theta = \mathrm{diag}(e^{i\theta m_J})$ produces three distinct phase factors $\{1, e^{i\theta}, e^{-i\theta}\}$ on the conjugate $U_\theta T U_\theta^{-1}$. **T4 identifies this as exactly the $\Delta m_J \in \{-1, 0, +1\}$ grading.** The $m_J$-resolved refinement that L5 named as "the load-bearing data for a non-trivial smash-product" is the substrate where the $U(1)$ action is diagonal.

### 4.2 The $U(1)$ embedding into $SL_2$

The $U(1) \subset SL_2$ corresponds to the maximal torus of $SU(2)$. The $\Delta m_J$ grading on $\mathcal{A}^{m_J\text{-OD}}$ is the weight grading under this $U(1)$, and the full $SL_2$ action additionally includes the ladder operators $J_+, J_-$ which raise/lower $m_J$ by $\pm 1$.

For the full $SL_2$ action on $\mathcal{A}^{m_J\text{-OD}}$:
- $J_z$ acts diagonally by $\Delta m_J$ on transitions.
- $J_+$ raises $\Delta m_J$ by 1 (sends $\mathcal{A}_k \to \mathcal{A}_{k+1}$).
- $J_-$ lowers $\Delta m_J$ by 1.

This gives the $\mathcal{A}^{m_J\text{-OD}}$ substrate the structure of an $\mathfrak{sl}_2$-module. The smash product $\mathcal{A}^{m_J\text{-OD}} \rtimes \mathcal{O}(SL_2)$ would carry both the $\mathfrak{sl}_2$-module structure AND the Hopf-coproduct from the OffDiag substrate (Sprint Q5'-OffDiag-Dirac, v3.62.0).

### 4.3 L2 extension to $j_{\max} = 3/2$ is the sprint-scale prerequisite

The decorated-PW substrate $\mathcal{H}_{\mathrm{dec}}^{(j_{\max})}$ (Sprint Q5'-Decorated-PW, v3.63.0) at $j_{\max} = 1/2$ covers only spin-$1/2$ matrix coefficients. The $m_J$-resolved OffDiag at $n_{\max} = 2$ contains spin-$3/2$ states from sectors $(2, 1)$ and $(2, 2)$, requiring matrix coefficients $\pi^{3/2, k}_{m, n}$ in the L2 substrate.

Sprint-scale extension: $\mathcal{H}_{\mathrm{dec}}^{(j_{\max} = 3/2)}$ has dimension $3 \cdot (1 + 4 + 9) = 42$. The Hopf-axiom panel for the extended substrate is sprint-scale ~1 day (analogous to the original L2 verification at $j_{\max} = 1$ in v3.63.0).

---

## 5. Verification gate compliance (§13.4)

- **Test gate ✓** — no production code touched; bit-exact `sympy.Rational` throughout.
- **Dead-end gate ✓** — no §3 match.
- **Prime directive gate ✓** — no discrete-structure modifications.
- **Consistency gate ✓** — reproduces v3.63.0 L5 finding (three phase factors) bit-exactly via direct $\Delta m_J$ computation; reproduces FockSpectralTriple state structure (16 dim, 5 sectors, $j$-content $\{1/2, 3/2\}$).
- **Equation gate ✓** — Wigner-Eckart selection $|\Delta m_J| \le 1$ verified 100/100; $\mathbb{Z}_3$-grading dimensions $30 + 40 + 30 = 100$ bit-exact; $U(1)$ z-rotation phase factors $\{1, e^{i\theta}, e^{-i\theta}\}$ symbolic match to L5.

---

## 6. Honest scope

### 6.1 Closed at theorem grade (bit-exact at finite cutoff $n_{\max} = 2$)

- $\mathbb{Z}_3$-grading by $\Delta m_J \in \{-1, 0, +1\}$.
- Symmetric $\pm 1$ dimension count $30 = 30$.
- $U(1)$ z-rotation diagonalisation.
- L5 three-phase finding reproduction.
- L2 extension prerequisite identification ($j_{\max} \ge 3/2$ required).

### 6.2 Open follow-ons

- **Sprint-scale (~1 day):** Extend L2 substrate to $j_{\max} = 3/2$. Verify Hopf-axiom panel bit-exactly on the extended substrate. Predicted dim: 42 at $j_{\max} = 3/2$, $n_k = 3$.
- **Sprint-scale (~2-3 days):** Construct $\mathcal{A}^{m_J\text{-OD}} \rtimes \mathcal{O}(SL_2)$ at $n_{\max} = 2$, $j_{\max} = 3/2$. The smash product has dimension $\dim \mathcal{A}^{m_J\text{-OD}} \cdot \dim \mathcal{O}(SL_2)^{\otimes 3}|_{j_{\max} = 3/2} = 100 \cdot 42 = 4200$ at the substrate level (massive expansion).
- **Sprint-scale (~1-2 days):** Bit-exact verification of Hopf-axiom panel for the smash product. Confirms or falsifies the non-trivial smash-product extension beyond Levi.
- **Multi-year:** Generalise to $n_{\max} = 3$ (where j-content extends to $j = 5/2$ from sector $(3, 2)$ $\kappa = -3$; needs $j_{\max} \ge 5/2$ in L2). Generalise the smash-product to the full cosmic-Galois $U^*$ extension (Sprint T5 of this sprint).

### 6.3 Hard prohibitions (§13.5) clean

- No changes to natural geometry hierarchy.
- No fitted/empirical parameters.
- No deletion of negative results.
- Paper 2 combination-rule "conjectural" label unchanged.

### 6.4 Curve-fit audit (`feedback_audit_numerical_claims`)

The $\Delta m_J \in \{-1, 0, +1\}$ grading is FORCED by Wigner-Eckart for the vector operator $A$; no fitting. The $30/40/30$ count is bit-exact from the spectral triple definition. No free parameters.

### 6.5 Tag-transcendentals (`feedback_tag_transcendentals`)

No transcendentals introduced. All matrix elements in $\kappa \cdot \mathbb{Q}$.

### 6.6 Discrete-for-skeleton compliance

All bit-exact `sympy.Rational`. No PSLQ.

### 6.7 WH1 PROVEN unaffected

This sprint constructs the substrate for a multi-year extension; does not modify any propinquity result.

---

## 7. Files

### Produced
- `debug/compute_q5p_mj_smash_product.py` — driver (~200 lines, 0.01 s wall, bit-exact `sympy`).
- `debug/data/sprint_q5p_mj_smash_product.json` — bit-exact data dump.
- `debug/sprint_q5p_mj_smash_product_memo.md` — this memo.

### Used (load-bearing inputs)
- `debug/sprint_q5p_combined_substrate_memo.md` (v3.63.0 L5 — three-phase finding).
- `debug/sprint_q5p_decorated_pw_memo.md` (v3.63.0 L2 — $\mathcal{H}_{\mathrm{dec}}$ substrate).
- `debug/sprint_q5p_offdiag_dirac_memo.md` (v3.62.0 T3b — OffDiag substrate).
- `geovac/spectral_triple.py` (FockSpectralTriple class with $m_J$-resolved labels).

---

## 8. One-line verdict

The $m_J$-resolved OffDiag substrate at $n_{\max} = 2$ admits a bit-exact $\mathbb{Z}_3$-grading by $\Delta m_J \in \{-1, 0, +1\}$ with symmetric dim count $30/40/30$ over 100 non-zero off-diagonal Dirac entries, **diagonalising the $U(1)$ z-rotation action** that L5 (v3.63.0) identified as broken on the basic substrate via three phase factors. The j-content $\{1/2, 3/2\}$ at $n_{\max} = 2$ requires extending the L2 decorated-PW substrate to $j_{\max} \ge 3/2$ for full $SU(2)$ action on the $m_J$-resolved substrate — a sprint-scale prerequisite (~1 day) to the multi-year smash-product construction $\mathcal{A}^{m_J\text{-OD}} \rtimes \mathcal{O}(SL_2)$. **Foundation laid bit-exactly for the non-trivial Hopf extension beyond Levi-decomposition.**
