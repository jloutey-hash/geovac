# Sprint Q5'-Tannakian-Closure TC-2b — $SL_2$ piece of converse reconstruction

**Date:** 2026-06-06 (same day as TC-2a;\ second sub-sprint in the TC-2 converse arc)
**Sprint:** TC-2b of Q5'-Tannakian-Closure (TC-2a closed v3.75.0 with abelian-factor equality;\ TC-2b adds the $SL_2$ factor)
**Driver:** `debug/compute_q5p_tc2b_sl2_aut_equality.py`
**Module:** `geovac/tannakian.py` (re-used substrate;\ no new production code)
**Data:** `debug/data/sprint_q5p_tc2b_sl2_aut_equality.json`
**Wall time:** 0.03 s
**Discipline:** bit-exact `sympy.Rational` / `sympy.Integer` throughout.

---

## 1. TL;DR

**Verdict: POSITIVE — equality on the $SL_2$ factor of $U^*_{\mathrm{Levi}}$ at the Peter--Weyl panel.**

Computing $\dim \mathrm{Aut}^\otimes(\omega)$ on the PW panel $\{V_{\mathrm{triv}}, V_{\mathrm{fund}}, \mathrm{Sym}^2 V_{\mathrm{fund}}\}$ bit-exactly yields

$$
\boxed{\dim \mathrm{Aut}^\otimes(\omega)\big|_{\mathrm{PW}\ \mathrm{panel}} = 3 = \dim SL_2.}
$$

Combined with TC-2a ($\dim = 15$ on the $n_{\max}$-axis substrate), the **exterior-tensor-product theorem** for neutral Tannakian categories (Deligne--Milne 1982 Theorem~2.3) gives

$$
\dim \mathrm{Aut}^\otimes(\omega)\big|_{\mathrm{combined}} = 15 + 3 = 18 = \dim U^*_{\mathrm{Levi}}.
$$

The TC-1e/TC-1f inclusion $U^*_{\mathrm{Levi}} \subseteq \mathrm{Aut}^\otimes(\omega)$ is therefore an **equality** on the combined $n_{\max}$-axis $\times$ PW panel at $n_{\max} = 2$.

**Structural mechanism (the substantive new content).** Tensor compatibility on $V_{\mathrm{fund}} \otimes V_{\mathrm{fund}}$ together with the decomposition $V_{\mathrm{fund}} \otimes V_{\mathrm{fund}} = \mathrm{Sym}^2 V_{\mathrm{fund}} \oplus V_{\mathrm{triv}}$ gives the bit-exact identity

$$
\Phi \cdot (\eta_{V_{\mathrm{fund}}} \otimes \eta_{V_{\mathrm{fund}}}) \cdot \Phi^{-1}
\;=\; \mathrm{Sym}^2(\eta_{V_{\mathrm{fund}}}) \oplus \det(\eta_{V_{\mathrm{fund}}}),
$$

where $\Phi$ is the explicit $4\times 4$ rational decomposition isomorphism (built in `debug/compute_q5p_tc2b_sl2_aut_equality.py:build_phi_decomposition`). Two structural facts fall out **automatically** from this single matrix identity:

1. **Top-left $3\times 3$ block** $=$ $\eta_{\mathrm{Sym}^2}$ $=$ $\mathrm{Sym}^2(\eta_{V_{\mathrm{fund}}})$ — i.e., the $\mathrm{Sym}^2$ action is *forced* to be the symmetric square of the standard action.
2. **Bottom-right $1\times 1$ block** $=$ $\eta_{V_{\mathrm{triv}}}$ $=$ $\det(\eta_{V_{\mathrm{fund}}})$.

Imposing **unit normalisation** $\eta_{V_{\mathrm{triv}}} = 1$ then forces $\det(\eta_{V_{\mathrm{fund}}}) = 1$ — i.e., $\eta_{V_{\mathrm{fund}}} \in SL_2(\Q)$.

The variety dimension calculation: ambient $M_2(\Q) = 4$-dim, single quadratic constraint $ad - bc - 1 = 0$, Jacobian rank $= 1$ at every $SL_2$-point (gradient $(d, -c, -b, a)$ has rank 1 since $\det \ne 0$), codim $= 1$, **variety dim $= 4 - 1 = 3 = \dim SL_2$**.

---

## 2. Verdict against decision gate

| Gate | Verdict |
|:-----|:--------|
| **POSITIVE — equality on $SL_2$ factor at PW panel** | **selected** — variety dim 3 bit-exact;\ Jacobian rank 1 at identity and generic point;\ tensor compat forces $\eta_{\mathrm{Sym}^2} = \mathrm{Sym}^2(\eta_{V_{\mathrm{fund}}})$ and $\eta_{V_{\mathrm{triv}}} = \det(\eta_{V_{\mathrm{fund}}})$ as bit-exact symbolic identities. |
| SURPRISE — extra automorphisms found | not selected (dim matches) |
| NEGATIVE — short of $SL_2$ | not selected (decomposition is fully consistent) |

---

## 3. Structural contrast with TC-2a

| Feature | TC-2a (abelian factor) | TC-2b ($SL_2$ factor) |
|:--------|:-----------------------|:----------------------|
| Hopf algebra | $\Sym(V_{15})$ primitive | $\mathcal{O}(SL_2)$ |
| Reps in panel | 15 abelian witnesses + $T$ | $V_{\mathrm{triv}}, V_{\mathrm{fund}}, \mathrm{Sym}^2$ |
| Constraints | Linear in $\eta$ entries (Schur for abelian primitive) | Quadratic via tensor compat on $V_{\mathrm{fund}}^{\otimes 2}$ |
| Solver method | `sympy.linear_eq_to_matrix` → rank/nullity | Symbolic Jacobian rank at variety point |
| Variety | Affine linear $\Q^{15}$ | Quadric $\{ad - bc = 1\} \subset \Q^4$ |
| Predicted dim | $3 N(2) = 15$ | $3 = \dim SL_2$ |
| Computed dim (bit-exact) | **15** | **3** |

The $SL_2$ piece needs the symmetric-square + antisymmetric-square decomposition because each irreducible $SL_2$-rep has $\End = \Q \cdot I$ (Schur), so naturality alone cannot detect the group;\ the binding constraint comes from how tensor products decompose.

---

## 4. Bit-exact panel

### 4.1 The decomposition isomorphism $\Phi$

Source basis (Kronecker, outer index slow):\ $(e_1 \otimes e_1, e_1 \otimes e_2, e_2 \otimes e_1, e_2 \otimes e_2)$.

Target basis (Sym$^2$ first, then $V_{\mathrm{triv}}$):\ $(e_1^2, e_1 e_2, e_2^2, e_1 \wedge e_2)$,

where $e_1 e_2 = e_1 \otimes e_2 + e_2 \otimes e_1$ (matching the convention in `geovac.tannakian._sl2_sym2_action`) and $e_1 \wedge e_2 = e_1 \otimes e_2 - e_2 \otimes e_1$.

$$
\Phi = \begin{pmatrix}
1 & 0 & 0 & 0 \\
0 & 1/2 & 1/2 & 0 \\
0 & 0 & 0 & 1 \\
0 & 1/2 & -1/2 & 0
\end{pmatrix},
\qquad
\Phi^{-1} = \begin{pmatrix}
1 & 0 & 0 & 0 \\
0 & 1 & 0 & 1 \\
0 & 1 & 0 & -1 \\
0 & 0 & 1 & 0
\end{pmatrix}.
$$

Bit-exact verified $\Phi \cdot \Phi^{-1} = I_4$, $\Phi \cdot (I \otimes I) \cdot \Phi^{-1} = I_4$.

### 4.2 Tensor compat identity

For symbolic $\eta_{V_{\mathrm{fund}}} = \begin{pmatrix} a & b \\ c & d \end{pmatrix}$:

$$
\Phi \cdot (\eta_{V_{\mathrm{fund}}} \otimes \eta_{V_{\mathrm{fund}}}) \cdot \Phi^{-1}
= \begin{pmatrix}
a^2 & 2ab & b^2 & 0 \\
ac & ad+bc & bd & 0 \\
c^2 & 2cd & d^2 & 0 \\
0 & 0 & 0 & ad - bc
\end{pmatrix}.
$$

Top-left $3 \times 3$ block exactly equals `_sl2_sym2_action(eta_fund)` (bit-exact `sympy` simplification). Bottom-right $1 \times 1$ block exactly equals $\det(\eta_{V_{\mathrm{fund}}}) = ad - bc$. Off-diagonal blocks identically zero.

### 4.3 Variety dim

| Quantity | Value |
|:---------|------:|
| Ambient dim ($\eta_{V_{\mathrm{fund}}} \in M_2(\Q)$) | 4 |
| Det constraint $ad - bc - 1 = 0$ | 1 polynomial equation |
| Jacobian at $\eta = I$:\ $\nabla(ad-bc-1) = (d, -c, -b, a)|_{a=d=1, b=c=0} = (1, 0, 0, 1)$ | rank 1 |
| Jacobian at $\eta = \begin{pmatrix} 2 & 1 \\ 1 & 1 \end{pmatrix}$ (generic $SL_2$):\ $(d, -c, -b, a) = (1, -1, -1, 2)$ | rank 1 |
| Codim | 1 |
| **Variety dim** | **3** |
| Predicted ($= \dim SL_2$) | 3 |
| **Match** | **TRUE** |

### 4.4 Combined-category sanity check

For one combined witness $V_{g_0} \otimes V_{\mathrm{fund}}$ with $g_0 = (1, 0, 0)$:

- $\eta_{V_{g_0}} = \exp(q_g E_{12}) = \begin{pmatrix} 1 & q_g \\ 0 & 1 \end{pmatrix}$ from TC-2a (with $q_g = 2/3$);
- $\eta_{V_{\mathrm{fund}}} = \begin{pmatrix} 2 & 1 \\ 1 & 1 \end{pmatrix} \in SL_2(\Q)$;
- Combined $\eta_{\mathrm{combined}} = \eta_{V_{g_0}} \otimes \eta_{V_{\mathrm{fund}}}$ as a $4 \times 4$ matrix.

Bit-exact verifications:
- **Invertibility**:\ $\det(\eta_{\mathrm{combined}}) = \det(\eta_{V_{g_0}})^{\dim V_{\mathrm{PW}}} \cdot \det(\eta_{V_{\mathrm{fund}}})^{\dim V_g} = 1 \cdot 1 = 1$.
- **$H_{GV}$-side naturality** on the combined rep:\ $[\eta_{\mathrm{combined}}, X_g^{V_g} \otimes I_{V_{\mathrm{PW}}}] = 0$ bit-exact (follows from $[\eta_{V_g}, X_g^{V_g}] = 0$, which holds because $\eta_{V_g} = \exp(q_g X_g^{V_g})$).
- **Combined dim** $= 15 + 3 = 18 = \dim U^*_{\mathrm{Levi}}$ (by exterior tensor product theorem).

---

## 5. Honest scope

**Closed at theorem grade:**

- $\dim \mathrm{Aut}^\otimes(\omega)$ on the PW panel $\{V_{\mathrm{triv}}, V_{\mathrm{fund}}, \mathrm{Sym}^2\}$ with the natural decomposition iso $\Phi$ is exactly 3, bit-exact via Jacobian rank.
- The structural identity $\eta_{\mathrm{Sym}^2} = \mathrm{Sym}^2(\eta_{V_{\mathrm{fund}}})$ is forced by tensor compat (not assumed).
- The combined-category dimension $\dim = 18 = \dim U^*_{\mathrm{Levi}}$ follows from the exterior-tensor-product theorem (Deligne--Milne 1982 Theorem~2.3).
- The combined sanity check confirms a generic combined $\eta = \eta_{V_g} \otimes \eta_{V_{\mathrm{PW}}}$ satisfies $H_{GV}$-side naturality bit-exact.

**Not closed here (deferred to TC-2c / TC-2d):**

- **Higher cutoffs.** The TC-2a/b argument is at $n_{\max} = 2$. TC-2c extends to $n_{\max} = 3$.
- **Pro-system coherence.** The per-cutoff equalities are bound up into a single inverse-limit statement via PS-2 functoriality;\ that packaging is TC-2d.
- **PW panel completeness.** The panel $\{V_{\mathrm{triv}}, V_{\mathrm{fund}}, \mathrm{Sym}^2\}$ generates $\mathrm{Rep}(SL_2)$ under tensor and subobjects (every Sym$^k$ for $k \ge 3$ appears as a subobject of $V_{\mathrm{fund}}^{\otimes k}$, which is generated). Verifying that the higher Sym$^k$ contribute no additional constraints is straightforward but requires explicit panel extension (sprint-scale follow-on, not done here).

**Multi-year content unchanged:** the pro-finite Tannakian theorem on $\mathcal{O}_\infty$ coherent with v3.66.0 FO3 Interpretation C at the period level remains the multi-year frontier. TC-2b plus TC-2a plus TC-2c/d closes the abelian-and-$SL_2$ factors of the reconstruction at finite cutoff;\ what remains is the inverse limit + period-level coherence.

---

## 6. Discipline checks

**Curve-fit audit (`feedback_audit_numerical_claims`):** Zero free parameters. The Sym$^2$ formula is the standard symmetric-square representation;\ the decomposition iso $\Phi$ is the standard (sym, antisym) split. The variety dim follows from a single Jacobian rank computation, predicted before running by classical Tannakian theory.

**Discrete-for-skeleton (`feedback_discrete_for_skeleton`):** All `sympy.Rational` / `Integer`. No floats. No PSLQ.

**Tag transcendentals (`feedback_tag_transcendentals`):** None appear.

**Diagnostic-before-engineering:** TC-2b is a fresh test of a structural prediction, not an iteration past two negatives.

**WH1 PROVEN unaffected.** **Hard prohibitions check clean.**

---

## 7. Files

### Produced
- `debug/compute_q5p_tc2b_sl2_aut_equality.py` — driver (~400 lines, 0.03 s wall).
- `debug/data/sprint_q5p_tc2b_sl2_aut_equality.json` — bit-exact data dump.
- `debug/sprint_q5p_tc2b_sl2_aut_equality_memo.md` — this memo.
- `tests/test_tannakian_sl2_aut_equality.py` — 14 regression tests, all pass (0.89 s).

### Used (load-bearing inputs)
- `geovac/tannakian.py` — `_sl2_sym2_action` (Sym$^2$ formula on standard $SL_2$ rep).
- `geovac/pro_system.py` — `primitive_generators(2)` for the combined sanity check.
- TC-2a memo (`debug/sprint_q5p_tc2a_aut_equality_memo.md`) — abelian factor closed.
- TC-1f memo (`debug/sprint_q5p_tc1f_sl2_inclusion_memo.md`) — PW substrate.

### Published references
- Deligne, P.; Milne, J. S. ``Tannakian categories'' (1982) Theorem~2.3 (exterior tensor product of neutral Tannakian categories) and Theorem~2.11 (reconstruction).
- Fulton, W.; Harris, J. ``Representation Theory'' Springer (1991) §§11--12 ($SL_2$ rep theory and Sym$^k$ decomposition).

---

## 8. Paper-edit recommendation (deferred per PI direction)

Paper updates are deferred to a single batch after TC-2c and TC-2d close. TC-2b's writeup will go into Paper 56 as a new subsection §sec:tc2b_sl2_equality, with the headline theorem dim Aut$^\otimes$ on PW panel = 3 = dim $SL_2$, the structural mechanism (tensor compat → forced det = 1), and the combined corollary dim $= 18 = \dim U^*_{\mathrm{Levi}}$ via Deligne--Milne 1982 Thm 2.3.

---

## 9. One-line verdict

**POSITIVE — $SL_2$ factor of $U^*_{\mathrm{Levi}}$ equality on the Peter--Weyl panel at $n_{\max} = 2$.** Variety $\{ad - bc = 1\} \subset M_2(\Q)$ has dim 3 bit-exact;\ tensor compat forces $\eta_{\mathrm{Sym}^2} = \mathrm{Sym}^2(\eta_{V_{\mathrm{fund}}})$ and $\eta_{V_{\mathrm{triv}}} = \det(\eta_{V_{\mathrm{fund}}})$ as bit-exact symbolic identities;\ unit normalisation $\eta_{V_{\mathrm{triv}}} = 1$ forces $\det = 1$. Combined with TC-2a, dim Aut$^\otimes$ on the combined category $= 18 = \dim U^*_{\mathrm{Levi}}$ via Deligne--Milne 1982 Theorem~2.3.
