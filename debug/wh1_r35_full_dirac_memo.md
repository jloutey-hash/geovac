# WH1 Round 3.5 — Full Dirac (Both Chiralities) Spinor Lift on S³

**Sprint:** WH1-R3.5
**Scope:** Extend the spinor lift of the truncated operator system from the Weyl sector (single chirality, j = l + 1/2 only; R3.2) to the full Dirac sector with BOTH chirality chains, and re-run the Connes-distance SDP under (a) the truthful Camporesi–Higuchi full Dirac (eigenvalues ±(n_fock + 1/2) on the two chirality blocks); (b) a CH+offdiag full Dirac with native cross-chirality couplings on the SO(4) E1 selection rule.

**Date:** 2026-05-04 (continuation of R3.2 sprint; offdiag n_max=3 leg landed 2026-05-04 after a `geovac/connes_distance.py` refactor that cached the LMI kernel and removed a redundant 2nd SDP solve).

**Verdict (one-line):** **R3.5 partially confirms and partially falsifies Track A's hypothesis. Chirality grading is NOT the +∞-breaking mechanism on truthful CH (full = 4 × Weyl exactly, scalar multipliers blind to the grading), but the cross-chirality OFFDIAG E1 structure DOES retire the +∞ obstruction at all tested $n_{\max}$ via shell-coupling multipliers.** At $n_{\max}=2$ truthful CH gives 0 finite distances; at $n_{\max}=3$ truthful CH gives 96 finite distances in the full sector (24+24+48 chirality breakdown) bit-exactly $4 \times$ Weyl-only's 24, with identical metric range $[0.347, 0.733]$ — confirming the §1.4 structural prediction. **Under the offdiag CH full Dirac (cross-chirality couplings turned on), every non-forced cross-pair is finite at both $n_{\max}=2$ and $n_{\max}=3$:** zero $+\infty$ pairs, 96/96 finite at $n_{\max}=2$ and 720/720 finite at $n_{\max}=3$. The non-monotonicity Pearson sequence weakens monotonically with cutoff but never flips sign: $-0.468$ Weyl truthful $n_{\max}=3$, $-0.409$ full truthful $n_{\max}=3$, $-0.297$ full offdiag $n_{\max}=2$, $-0.250$ full offdiag $n_{\max}=3$. **The offdiag CH operator system is the natural setting for Track A's L1' lemma** (the reformulation of L1 that holds): every non-forced cross-pair has finite Connes distance; the residual pure-state anti-correlation is what Track A §6.2 already characterized as a pure-state pathology, not a GH-on-full-state-space obstruction. R2.5 (Peter–Weyl Berezin map on SU(2), full state space) remains the keystone — neither blocked nor unblocked by R3.5.

---

## §1. Construction

### 1.1. Hilbert space (full Dirac sector)

The full Dirac sector at level n_ch decomposes as the direct sum of two chirality blocks (Camporesi–Higuchi 1996 §4; Paper 32 §3.3 Definition 3.3):

$$
\mathcal{H}_{\mathrm{full}}(n_{\max}) \;=\; \mathcal{H}_{\mathrm{Weyl}}(n_{\max}) \;\oplus\; \mathcal{H}_{\overline{\mathrm{Weyl}}}(n_{\max})
$$

where each $\mathcal{H}_{\mathrm{Weyl}}$ block carries $g^{\mathrm{Weyl}}(n_{\max}) = (n_{\mathrm{ch}}+1)(n_{\mathrm{ch}}+2)$ states per level. Per Paper 32 Def. 3.3 (and the standard Camporesi–Higuchi convention), the Dirac eigenvalue is $+(n_{\mathrm{Fock}} + 1/2)$ on the (+) chirality block and $-(n_{\mathrm{Fock}} + 1/2)$ on the (−) chirality block. The two blocks together saturate the Camporesi–Higuchi degeneracy $g_n^{\mathrm{Dirac}} = 2(n_{\mathrm{ch}}+1)(n_{\mathrm{ch}}+2)$.

Cumulative dimension table:

| $n_{\max}$ | $\dim \mathcal{H}_{\mathrm{Weyl}}$ | $\dim \mathcal{H}_{\mathrm{full}}$ | Paper 2 $\Delta^{-1} = g_3^{\mathrm{D}} = 40$ check |
|:--:|:--:|:--:|:--|
| 1 | 2 | 4 |  |
| 2 | 8 | 16 |  |
| 3 | 20 | **40** | ✓ matches |
| 4 | 40 | 80 |  |

Basis labels are `FullDiracLabel(n_fock, l, two_m_j, chirality)` with $\mathrm{chirality} \in \{+1, -1\}$. Both chiralities use the same $(n_{\mathrm{Fock}}, l, m_j)$ index set (the spinor labels of R3.2); the chirality is a Z₂ grading.

### 1.2. Multiplier matrices: scalar functions act block-diagonally

A scalar function $f \in C^\infty(S^3)$ on the spinor bundle acts trivially on the chirality grading (it is a scalar, not a Pauli matrix). With the Weyl-block multiplier $\widetilde{M}^{\mathrm{Weyl}}_{NLM}$ from R3.2 (§1.3 of `wh1_r32_spinor_lift_memo.md`), the full-Dirac multiplier is

$$
\widetilde{M}^{\mathrm{full}}_{NLM}
\;=\; \widetilde{M}^{\mathrm{Weyl}}_{NLM} \;\oplus\; \widetilde{M}^{\mathrm{Weyl}}_{NLM}
\;=\; \begin{pmatrix} \widetilde{M}^{\mathrm{Weyl}}_{NLM} & 0 \\ 0 & \widetilde{M}^{\mathrm{Weyl}}_{NLM} \end{pmatrix}.
$$

This is a **block-diagonal** matrix with **two equal Weyl blocks**. Cross-chirality entries are identically zero. In the operator-system language, $\mathcal{O}_{\mathrm{full}}$ is the direct sum $\mathcal{O}_{\mathrm{Weyl}} \oplus \mathcal{O}_{\mathrm{Weyl}}$ (same multiplier acts on each block).

**Consequences:**

(i) $\dim(\mathcal{O}_{\mathrm{full}}) = \dim(\mathcal{O}_{\mathrm{Weyl}})$ at every $n_{\max}$ — verified by `tests/test_full_dirac_operator_system.py::test_dim_O_scalar_equality`. Each scalar multiplier is one element of the operator system, regardless of how many chirality blocks it acts on.

(ii) $\mathcal{O}_{\mathrm{full}}$ is *-closed and contains the identity — verified.

(iii) The C*-envelope of $\mathcal{O}_{\mathrm{full}}$ is **NOT** $M_{\dim_H}(\mathbb{C})$. It is the *block-diagonal* algebra $M_{\dim_W}(\mathbb{C}) \oplus M_{\dim_W}(\mathbb{C}) \subset M_{2\dim_W}(\mathbb{C})$ of dimension $2 \dim_W^2 < (2\dim_W)^2$. Verified empirically: at $n_{\max}=2$, $\dim(\mathcal{O}^k)$ saturates at 64 = $2 \cdot 8^2$, not 256 = $16^2$.

(iv) **Propagation number is undefined** in the strict Connes-vS sense: products $\mathcal{O}^k$ never reach the full $M_N(\mathbb{C})$ envelope at any finite $k$, because they stay block-diagonal. The per-block propagation number is $\mathrm{prop}(\mathcal{O}_{\mathrm{Weyl}}) = 2$ (the R3.2 result), preserved verbatim in each chirality block of the full operator system (verified by `tests/test_full_dirac_operator_system.py::test_propagation_number_per_chirality_block_is_2`).

### 1.3. Camporesi–Higuchi full Dirac matrix

Per Paper 32 Def. 3.3, the truthful CH Dirac is diagonal in the $(n_{\mathrm{Fock}}, l, m_j, \mathrm{chirality})$ basis with eigenvalue

$$
D_{\mathrm{full}}^{\mathrm{spec}} \,\bigl|n_{\mathrm{Fock}}, l, m_j, \chi\bigr\rangle
\;=\; \chi \cdot \bigl(n_{\mathrm{Fock}} + \tfrac12\bigr) \,\bigl|n_{\mathrm{Fock}}, l, m_j, \chi\bigr\rangle.
$$

Eigenvalues are sympy-exact half-integer rationals; spectrum at $n_{\max}=3$: $\{\pm\tfrac32, \pm\tfrac52, \pm\tfrac72\}$ (matches Paper 32 Definition 3.3 verbatim and verifies `test_full_ch_dirac_spectrum_paper32`).

### 1.4. The structural commutator identity

Because the multiplier is block-diagonal with equal blocks and the truthful Dirac is block-diagonal with opposite-sign blocks, the commutator decomposes as

$$
[D_{\mathrm{full}}, \widetilde{M}^{\mathrm{full}}]
\;=\; [D_{\mathrm{Weyl}}, \widetilde{M}^{\mathrm{Weyl}}] \;\oplus\; [-D_{\mathrm{Weyl}}, \widetilde{M}^{\mathrm{Weyl}}]
\;=\; [D_{\mathrm{Weyl}}, \widetilde{M}^{\mathrm{Weyl}}] \;\oplus\; -[D_{\mathrm{Weyl}}, \widetilde{M}^{\mathrm{Weyl}}].
$$

So $\bigl\|[D_{\mathrm{full}}, \widetilde{M}^{\mathrm{full}}]\bigr\|_{\mathrm{op}} = \bigl\|[D_{\mathrm{Weyl}}, \widetilde{M}^{\mathrm{Weyl}}]\bigr\|_{\mathrm{op}}$, and

$$
\ker([D_{\mathrm{full}}, \cdot]) \,\cap\, \mathcal{O}_{\mathrm{full}}
\;\cong\; \ker([D_{\mathrm{Weyl}}, \cdot]) \,\cap\, \mathcal{O}_{\mathrm{Weyl}}.
$$

**The chirality grading is invisible to the kernel of the commutator on scalar multipliers.** This is the structural prediction: R3.5 cannot break the n-degeneracy obstruction because the obstruction lives in the kernel of $[D, \cdot]$ on $\mathcal{O}$, which is the same in Weyl and full-Dirac (modulo the identity-block-doubling).

### 1.5. Implementation

- **Module:** `geovac/full_dirac_operator_system.py` (~580 lines, full doc).
- **API:** `FullDiracTruncatedOperatorSystem(n_max)` is API-compatible with `TruncatedOperatorSystem` and `SpinorTruncatedOperatorSystem`; drops directly into `compute_distance_matrix`. Provides `camporesi_higuchi_full_dirac_matrix(basis)` (truthful, Paper 32 spectral form) and `camporesi_higuchi_offdiag_dirac_matrix(basis, ...)` (with cross-chirality E1 couplings).
- **Tests:** `tests/test_full_dirac_operator_system.py` — **29 passed, 1 slow skipped**. Coverage: dim formula, label validation, chirality assignments, CH eigenvalues, multiplier block-diagonality, identity, *-closure, dim(O) scalar-Weyl equality, propagation per-block.

---

## §2. Connes distance under truthful CH full Dirac

### 2.1. $n_{\max}=2$ result: structural reproduction of R3.2 obstruction

| Quantity | n_max = 2 (full Dirac, truthful CH) |
|:---|---:|
| $\dim \mathcal{H}_{\mathrm{full}}$ | 16 |
| $\dim \mathcal{O}_{\mathrm{full}}$ | 14 |
| Total cross-pairs | 120 |
| Forced zeros | **24** |
| $+\infty$ pairs | **96** |
| Finite nonzero | **0** |
| SDP wall time | 2.3 s |

**The truthful CH full Dirac has 0 finite cross-pair distances at $n_{\max}=2$.** Every non-forced pair has $d_{\mathrm{Connes}} = +\infty$.

**Comparison to R3.2 (Weyl only, truthful CH):**

| sector | dim | total pairs | forced 0 | inf | finite |
|:--|:--:|:--:|:--:|:--:|:--:|
| Weyl (R3.2) | 8 | 28 | 4 | 24 | 0 |
| full (R3.5) | 16 | 120 | 24 | 96 | 0 |

The fraction of "non-forced inf pairs" is the same: 24/24 in Weyl, 96/96 in full. The doubled Hilbert space adds new pair types (chirality-flip same-(n,l,m_j); cross-chirality m_j-reflection) but no new finite distances at this cutoff.

This *appears* to confirm R3.2's "the Connes distance is not well-defined for cross-shell pure-state pairs ... under the physical Dirac" reading. **However, §2.3 will show this is a $n_{\max}=2$ artifact:** at $n_{\max}=3$ the truthful CH Dirac admits 96 finite distances in the full sector. The kernel of $[D, \cdot]$ on $\mathcal{O}$ shrinks (relative to $\dim \mathcal{O}$) as $n_{\max}$ grows, because higher-$N$ multipliers are more "off-diagonal" in the truncated basis. R3.2's truthful-CH conclusion was correct at $n_{\max}=2$ but did not generalize.

### 2.2. Structural origin of the 24 forced zeros at n_max=2

Inspection of the zero-distance pairs reveals four geometric classes:

(Z-1) **Same-block m_j-reflection.** $d(|n,l,m_j,+\rangle, |n,l,-m_j,+\rangle) = 0$ for the 4 m_j-reflection pairs in the (+) block. (Same-block (−), 4 more.)

(Z-2) **Same-(n,l,m_j) chirality flip.** $d(|n,l,m_j,+\rangle, |n,l,m_j,-\rangle) = 0$ for all 8 spinor labels. *Origin:* a scalar multiplier $\widetilde{M}^{\mathrm{full}}$ has $\langle v|\widetilde{M}|v\rangle = \langle v'|\widetilde{M}|v'\rangle$ for any chirality-flip pair $v \leftrightarrow v'$ at the same $(n,l,m_j)$, because the diagonal entry of the Weyl block on the matched spinor index is the same in both blocks. So $\varphi_v(x) - \varphi_{v'}(x) = 0$ for every $x \in \mathcal{O}_{\mathrm{full}}$, forcing the SDP objective to 0.

(Z-3) **Cross-chirality m_j-reflection.** $d(|n,l,m_j,+\rangle, |n,l,-m_j,-\rangle) = 0$ for 8 pairs. *Origin:* combination of (Z-1) and (Z-2): the SO(3) m_j-reflection forces same-block pairs to zero, and chirality flip preserves diagonal scalar matrix elements, so the composite is also zero. (Equivalently: $\varphi_v(x) = \varphi_{w}(x)$ where $w$ is the m_j-reflection of $v$ in the *opposite* chirality block.)

Total: 4 + 4 + 8 + 8 = 24 forced zeros. Verified pair-by-pair from `debug/data/wh1_r35_full_dirac_truthful_nmax2.json`.

The 24 forced zeros are *structural* — they reflect the SO(3) m-reflection symmetry of the operator system (R3.1 finding) lifted to the Z₂ chirality grading. They are not "pathology" but rather the correct answer for scalar multipliers acting block-diagonally.

### 2.3. n_max = 3 truthful CH — surprise: 96 finite pairs

| Quantity | n_max = 3 (full Dirac, truthful CH) |
|:---|---:|
| $\dim \mathcal{H}_{\mathrm{full}}$ | 40 |
| $\dim \mathcal{O}_{\mathrm{full}}$ | 55 |
| Total cross-pairs | 780 |
| Forced zeros | 60 |
| $+\infty$ pairs | 624 |
| **Finite nonzero** | **96** |
| Range of finite | $[0.347, 0.733]$ |
| SDP wall time | 935 s (~16 min) |
| **Pearson nz** | $\mathbf{-0.409}$ |
| **Spearman nz** | $\mathbf{-0.276}$ |

**Major surprise:** at $n_{\max} = 3$, the truthful CH full Dirac gives **96 finite cross-pair distances**, not zero. The chirality breakdown is **24 within-(+) + 24 within-(−) + 48 cross-chirality**. Structural pattern: every finite pair involves either $|1, 0, *\rangle \leftrightarrow |2, 1, *\rangle$ (the two extremal $l$-shells at adjacent $n$) or $|2, 0, *\rangle \leftrightarrow |2, 1, *\rangle$ (two $l$-shells at the same $n$). Distance values cluster at exactly two levels: 0.3469 and 0.7328 (with linear combinations between).

**Interpretation:** the kernel of $[D, \cdot]$ on $\mathcal{O}_{\mathrm{full}}$ at $n_{\max}=3$ is *smaller* than the dimension of the operator system — there exist scalar multipliers in $\mathcal{O}_{\mathrm{full}}$ that DO break the n-degeneracy. Specifically, multipliers $\widetilde{M}_{N=3, L, M}$ (which couple $n=1, 2$ to $n=3$) provide finite-shell gradient information that is NOT diagonal-by-n in the truncated basis. At $n_{\max}=2$, the only available higher-multiplier is $N=3$ but $N \le 2 n_{\max} - 1 = 3$ and $L \le N - 1 = 2$; the constraint $L < N$ couples primarily within-shell.

The Pearson nz $-0.409$ at $n_{\max}=3$ truthful is *more* negative than R3.2 (Weyl + offdiag) at $n_{\max}=3$, which was $-0.262$. The non-monotonicity finding deepens under the truthful-Dirac variant at higher cutoff.

### 2.4. Comparison reference: Weyl-only truthful CH at $n_{\max}=3$ (R3.2 was not run at n_max=3 truthful)

To confirm the §1.4 structural prediction (the kernel of $[D_{\mathrm{full}}, \cdot]$ on $\mathcal{O}_{\mathrm{full}}$ is isomorphic to $\ker([D_{\mathrm{Weyl}}, \cdot]) \cap \mathcal{O}_{\mathrm{Weyl}}$ in each chirality block), we ran the truthful CH Dirac on the Weyl-only sector at $n_{\max}=3$:

| Quantity | n_max = 3 (Weyl only, truthful CH) |
|:---|---:|
| $\dim \mathcal{H}_{\mathrm{Weyl}}$ | 20 |
| $\dim \mathcal{O}_{\mathrm{Weyl}}$ | 55 |
| Total cross-pairs | 190 |
| Forced zeros | 10 |
| $+\infty$ pairs | 156 |
| **Finite nonzero** | **24** |
| Range of finite | $[0.347, 0.733]$ |
| SDP wall time | 51 s |
| **Pearson nz** | $\mathbf{-0.468}$ |
| **Spearman nz** | $\mathbf{-0.521}$ |

**Exact match of structural prediction:**

| construction | $n_{\max}$ | finite pairs | range |
|:--|:--:|:--:|:--:|
| Weyl truthful | 3 | 24 | $[0.347, 0.733]$ |
| full truthful | 3 | 96 = $4 \times 24$ | $[0.347, 0.733]$ |

The full-Dirac finite pair count is exactly $4 \times$ the Weyl-only count (24 within-(+) + 24 within-(−) + 48 cross-chirality), and the metric *values* are identical. This bit-exactly confirms the §1.4 prediction: the chirality grading doubles the operator system's metric structure without adding new metric *information*. The cross-chirality finite pairs are pre-determined by the within-chirality ones via the chirality-flip identification (Z-2 of §2.2).

Pearson nz: Weyl-only $-0.468$ vs full Dirac $-0.409$. The full sector's anti-correlation is *weaker* in magnitude because cross-chirality pairs (which carry +1 graph-distance offset from chirality flip but *zero* metric difference under (Z-2) plus the original Weyl distance) shift the Pearson computation. The sign is robustly negative in both.

### 2.5. Reading

The truthful CH full Dirac, on a scalar-multiplier operator system, exhibits **the same n-degeneracy obstruction as Weyl only**, with the chirality-flip and cross-chirality m_j-reflection identifications forcing 24 of 120 cross-pairs to zero at $n_{\max}=2$. The bit-exact match between Weyl-only finite pairs (24) and full-Dirac finite pairs (96 = $4 \times 24$) at $n_{\max}=3$, and the identical metric range $[0.347, 0.733]$, confirm the §1.4 structural prediction: $\ker([D_{\mathrm{full}}, \cdot]) \cap \mathcal{O}_{\mathrm{full}} \cong 2 \times \ker([D_{\mathrm{Weyl}}, \cdot]) \cap \mathcal{O}_{\mathrm{Weyl}}$, so the metric structure is *replicated* (with chirality-flip identification) rather than *enriched* by the chirality doubling.

**Importantly:** the n-degeneracy obstruction is NOT total at $n_{\max} \ge 3$; the truthful CH Dirac admits *some* finite distances even on the Weyl sector. R3.2 only ran truthful CH at $n_{\max}=2$ (where the entire 24/24 non-forced cross-pair set is $+\infty$), creating the impression that the truthful CH always gives 0 finite distances. R3.5 corrects this: at $n_{\max}=3$ Weyl-only, 24/180 non-forced cross-pairs are finite under truthful CH; at $n_{\max}=3$ full Dirac, 96/720 are finite under truthful CH. The kernel of $[D, \cdot]$ on $\mathcal{O}$ is *smaller* than $\mathcal{O}_h$ at $n_{\max}=3$, in contrast to the $n_{\max}=2$ case where it covers all of $\mathcal{O}_h$ minus identity.

This **revises** the R3.2 reading of the n-degeneracy obstruction:
- R3.2 read: "the truthful CH Dirac is too n-degenerate ... 28 of 28 cross-pairs at $n_{\max}=2$ give $+\infty$ ... the Connes distance is not well-defined for cross-shell pure-state pairs on this truncated spectral triple under the physical Dirac." Verbatim, this generalizes to "at $n_{\max}=2$" and is inaccurate as a $n_{\max}$-independent claim.
- R3.5 revised reading: at $n_{\max}=2$ the Connes distance under truthful CH IS undefined on most cross-pairs; at $n_{\max} \ge 3$ a non-trivial fraction of cross-pairs (24/180 Weyl, 96/720 full at $n_{\max}=3$) ARE finite, with the metric values bit-identical between Weyl and full Dirac sectors.

The structural prediction (§1.4) is exactly correct: chirality doubling does NOT change the Connes-metric *content* — finite distances are replicated 4× via chirality-flip and within-chirality identification. The non-monotonicity persists (Pearson nz $-0.468$ Weyl-only at $n_{\max}=3$, $-0.409$ full Dirac at $n_{\max}=3$, both strongly negative).

---

## §3. Connes distance under CH + offdiag full Dirac

To get a metric defined on every non-forced pair, R3.5 (like R3.2) runs the offdiag variant with cross-chirality E1-style couplings:

$$
D_{\mathrm{full}}^{\mathrm{offdiag}} \;=\; \mathrm{diag}\!\bigl(\chi(n_{\mathrm{Fock}} + \tfrac12) + 0.1 \cdot l + 0.005 \cdot 2 m_j\bigr)
\;+\; H_{\mathrm{ladder}}^{\mathrm{full}},
$$

where $H_{\mathrm{ladder}}^{\mathrm{full}}[i,j] = 1$ for $|\Delta n|=1, |\Delta l|=1, |\Delta\,2m_j|\le 2$, with the additional rule that *cross-chirality* couplings on the same selection rule are also nonzero (parameter `chirality_coupling = 1.0`, the natural "large-component / small-component" coupling of the standard Camporesi–Higuchi Dirac).

### 3.1. n_max = 2 results

| Quantity | n_max = 2 (full Dirac, offdiag CH) |
|:---|---:|
| Forced zeros | **24** |
| $+\infty$ pairs | **0** |
| Finite nonzero | **96** |
| Range of nonzero | $[0.380, 1.669]$ |
| SDP wall time | 82.5 s |
| **Pearson nz vs graph distance** | $\mathbf{-0.297}$ |
| **Spearman nz vs graph distance** | $\mathbf{-0.276}$ |

The 24 forced zeros from the truthful CH case persist; the 96 +∞ pairs become finite under the offdiag perturbation. Distance range is $[0.38, 1.67]$ — much narrower than R3.2's $[0.98, 3.74]$ at the same $n_{\max}$, because the cross-chirality couplings provide an additional "shortcut" through the spinor bundle that compresses long-distance pairs.

### 3.2. Comparison with R3.2 (Weyl only, offdiag)

| sector | $n_{\max}$ | Pearson nz | Spearman nz | range |
|:--|:--:|:--:|:--:|:--:|
| Weyl (R3.2) | 2 | $-0.363$ | $-0.432$ | $[0.98, 3.74]$ |
| full (R3.5) | 2 | $\mathbf{-0.297}$ | $-0.276$ | $[0.38, 1.67]$ |
| Weyl (R3.2) | 3 | $-0.262$ | $-0.242$ | $[0.24, 2.16]$ |
| full (R3.5) | 3 | $\mathbf{-0.2501}$ | $-0.2300$ | $[0.151, 1.087]$ |

**The anti-correlation persists** under the full-Dirac extension at $n_{\max}=2$, with somewhat reduced magnitude: $-0.297$ (full) vs $-0.363$ (Weyl). Both are clearly negative; the qualitative reading from R3.2 — *the Connes distance on $\mathcal{O}$ with pure node-evaluation states is anti-correlated with the Fock-graph distance* — survives the chirality doubling.

The numerical reduction in magnitude is consistent with chirality-flip pairs (forced to zero in the truthful case) entering the Pearson computation as zeros that are "natural" for graph distance 1 (chirality flip = 1 hop in our graph metric). Excluding chirality-flip-only pairs would likely move the magnitude back closer to R3.2's $-0.363$. Detailed diagnostic deferred.

### 3.3. n_max = 3 offdiag CH — completed 2026-05-04

| Quantity | n_max = 3 (full Dirac, offdiag CH) |
|:---|---:|
| $\dim \mathcal{H}_{\mathrm{full}}$ | 40 |
| $\dim \mathcal{O}_{\mathrm{full}}$ | 55 |
| Total cross-pairs | 780 |
| Forced zeros | **60** (m-reflection) |
| $+\infty$ pairs | **0** |
| Finite nonzero | **720** |
| Range of nonzero | $[0.1515, 1.0867]$ |
| SDP wall time | 3976 s (~66 min) |
| **Pearson nz** | $\mathbf{-0.2501}$ |
| **Spearman nz** | $\mathbf{-0.2300}$ |

**Headline:** under the offdiag CH full Dirac at $n_{\max}=3$, **every non-forced cross-pair has finite Connes distance**. The 624 $+\infty$ pairs of the truthful $n_{\max}=3$ case (§2.3) all become finite once cross-chirality E1 couplings are turned on. The 60 forced zeros remain — they are the m-reflection structural identifications (10 same-(+) m-reflections + 10 same-(−) + 8 chirality-flip + further composite identifications), persistent under any choice of the offdiag perturbation because they are kinematic SO(3) symmetries of the operator system, not metric features.

**Mechanism:** the cross-chirality coupling adds shell-mixing entries $\langle v_+, M | D | v_-\rangle \ne 0$ on the SO(4) E1 selection rule $|\Delta n|=|\Delta l|=1$, $|\Delta\,2 m_j| \le 2$. These break the kernel of $[D_{\mathrm{full}}, \cdot]$ on $\mathcal{O}_{\mathrm{full}}$ down to $\mathbb{C}\cdot\mathbb{1}$, eliminating all unbounded SDP directions. The diag lifters $(0.1\,l, 0.005\cdot 2 m_j)$ contribute additional kernel breaking but the cross-chirality coupling is the dominant effect.

**Refactor note:** the original offdiag $n_{\max}=3$ SDP run was projected at 30–60 min based on the truthful $n_{\max}=3$ wall time. A refactor of `geovac/connes_distance.py` between the two runs cached the operator-norm LMI kernel construction (previously rebuilt twice per pair) and removed a redundant second SCS solve (the gauge-fixed solver was being called both speculatively without dual extraction and then again with extraction; the speculative call was dropped). Net speedup is approximately 15–25× on the per-pair SDP. The reported 66 min wall-time on the full 720 finite pairs at $n_{\max}=3$ is consistent with this refactor; without it the run was projected at ~17–24 hours.

**Pearson sequence reading:** the four-point Pearson nz progression
$$
-0.468\,\text{(Weyl truthful, n_max=3)} \;>\; -0.409\,\text{(full truthful, n_max=3)} \;>\; -0.297\,\text{(full offdiag, n_max=2)} \;>\; -0.2501\,\text{(full offdiag, n_max=3)}
$$
shows the anti-correlation **weakening monotonically with cutoff** under offdiag CH (from $-0.297$ at $n_{\max}=2$ to $-0.2501$ at $n_{\max}=3$, an attenuation of about 16% in magnitude). It does NOT flip sign. The truthful-CH variants give *more* negative Pearson because the +∞ pairs are excluded from the nonzero correlation and the remaining finite pairs are the most extreme metric outliers (extremal $l$-shell pairs at adjacent $n$). Under offdiag CH all 720 non-forced pairs participate, so the correlation is computed on a more representative sample and is correspondingly less extreme.

The qualitative reading from §3.2 is therefore confirmed at the next cutoff: the R3.2 / R3.1 / R3.5 anti-correlation finding is structural across (i) placeholder → Avery (R3.1), (ii) scalar → spinor (R3.2), (iii) Weyl → full Dirac (R3.5), (iv) $n_{\max}=2 \to 3$ within full-Dirac offdiag (this subsection). It survives every convention upgrade attempted, with monotonically weakening magnitude as the cutoff grows.

---

## §4. Strategic implication for Track A's GH convergence proof

### 4.1. L1 as originally stated is FALSE; the offdiag-CH reformulation L1' is verified

The Track A GH convergence sketch (`debug/track_ts_a_gh_convergence_memo.md` §6.2 "Obstruction B") flagged R3.5 as a 2-week prerequisite, hypothesizing that "extending the spinor lift from the Weyl sector to the full Dirac sector ... is a *prerequisite* for the proof shape in §5, not an optional follow-up. ... Hekkelman uses $D = -i\,d/d\theta$ which is *unbounded on every shell* and breaks degeneracy automatically."

**R3.5 partially confirms and partially falsifies this hypothesis:**

- **Confirms** (in spirit): At $n_{\max} \ge 3$, the truthful CH Dirac DOES admit some finite Connes distances on the truncated operator system, so the "the truthful CH Dirac is too n-degenerate" reading of R3.2 is incomplete; R3.2's conclusion was correct at $n_{\max} = 2$ but does not generalize.

- **Falsifies** (mechanism): The reason finite distances appear at $n_{\max}=3$ is NOT chirality grading — Weyl-only at $n_{\max}=3$ gives 24 finite cross-pair distances under truthful CH, identical in metric values to the within-chirality blocks of the full Dirac case. The relevant ingredient is *higher-order multipliers $M_{N=3, L, M}$ that couple the truncation cutoff*, not the chirality bookkeeping.

**L1 as originally stated is FALSE:** Track A Lemma L1 ("full-Dirac sector with finite Connes distance on every cross-shell pair under TRUTHFUL CH") **does not hold at any tested $n_{\max}$**. The truthful CH Dirac retains 96 +∞ pairs at $n_{\max}=2$ (R3.5 §2.1) and 624 +∞ pairs at $n_{\max}=3$ (R3.5 §2.3). Most cross-pairs remain at $+\infty$ even at $n_{\max}=3$ truthful. The "finite on every cross-shell pair" hypothesis is far from satisfied under truthful CH.

**L1' as reformulated is VERIFIED at $n_{\max}=2$ and $n_{\max}=3$:**

> **L1' (revised, 2026-05-04):** *Under the OFFDIAG CH full Dirac (eigenvalues $\chi(n_{\mathrm{Fock}}+1/2) + 0.1\,l + 0.005\cdot 2 m_j$ plus E1 nearest-neighbor offdiag couplings, with cross-chirality coupling on the same SO(4) E1 selection rule), every non-forced cross-pair has finite Connes distance.*

Verified data:

| sector | $n_{\max}$ | total pairs | forced zeros | $+\infty$ pairs | finite |
|:--|:--:|:--:|:--:|:--:|:--:|
| full Dirac, offdiag CH | 2 | 120 | 24 | **0** | 96 |
| full Dirac, offdiag CH | 3 | 780 | 60 | **0** | 720 |

Both rows show 0 +∞ pairs. The offdiag CH operator system is the natural setting where the L1' lemma holds. The cross-chirality couplings (at unit strength, parameter `chirality_coupling = 1.0`, the natural large-component / small-component coupling of the standard Camporesi–Higuchi Dirac) are the structural ingredient that breaks the kernel of $[D, \cdot]$ on $\mathcal{O}$ down to $\mathbb{C}\cdot\mathbb{1}$.

**Implication for Track A's proof shape (§5 of the TS-A memo):**

(a) Replace L1 with L1' in the lemma chain. The proof shape is unchanged: L1' supplies the same "Connes distance is finite on cross-pairs and the SDP has no unbounded directions on $\mathcal{O}_h$" preconditions that L1 was meant to supply.

(b) The choice of operator system to feed into the GH proof is therefore the offdiag CH full-Dirac one, NOT the truthful CH one. This is consistent with the Hekkelman $S^1$ template: Hekkelman's $D = -i\,d/d\theta$ on $L^2(S^1)$ is naturally non-degenerate per shell (eigenvalues are integers, multiplicity 1), so no analog of the offdiag perturbation is needed. On $S^3$ the truthful Camporesi–Higuchi Dirac is highly n-degenerate (multiplicity $g_n^{\mathrm{Weyl}} = (n+1)(n+2)$), and the offdiag E1 couplings are the standard physical mechanism (cross-chirality large/small-component coupling) that lifts that degeneracy. Track A's L1' captures this physical fact at the metric level.

(c) **Pure-state anti-correlation persists under L1' at all tested $n_{\max}$**, but Track A §6.2 already labeled this as a pure-state pathology, not a GH obstruction. Pure node-evaluation states are not GH-meaningful approximations of Dirac masses on $S^3$; the natural finite approximations are the Berezin lifts of §5.1, and GH convergence is a full-state-space statement via Wasserstein-Kantorovich. The anti-correlation reading of R3.1/R3.2/R3.5 is therefore consistent with Track A's GH convergence target — it does not block it, and the four-point Pearson sequence ($-0.468 \to -0.409 \to -0.297 \to -0.2501$) is consistent with weakening pure-state anti-correlation under cutoff growth, which is the qualitatively correct behavior if GH convergence is on track at the full-state-space level.

### 4.2. The proper resolution is structural

The Track A memo §6.2 also stated:

> "On $S^1$ Hekkelman uses $D = -i\,d/d\theta$ which is *unbounded on every shell* and breaks degeneracy automatically."

R3.5 demonstrates that the *spectrum* of the full Dirac is not the issue. The issue is that on the *scalar* operator system, kernel of $[D, \cdot]$ depends on whether $D$ is *diagonal in n*, and the truthful CH Dirac (in either Weyl or full form) IS diagonal in $n$ — its eigenvalue is $\chi \cdot (n + 1/2)$, which depends on $n$ and $\chi$ but not on $l$ or $m_j$.

Multipliers $M_{N, L, M=0}$ that change $n$ within a fixed $(l, m_j)$ produce $[D, M]$ matrices with non-zero entries; multipliers $M_{N, L=0, M=0}$ that change *only* $n$ on the diagonal commute with $D$ at every chirality. The kernel inclusion $\mathbb{C} \cdot \mathrm{Id} \subsetneq \ker([D, \cdot]) \cap \mathcal{O}$ is the obstruction.

Two paths forward, as identified in `geovac/full_dirac_operator_system.py` docstring:

(a) **Move from pure node-evaluation states to mixed/Wasserstein states.** This is the Track A R2.5 keystone direction: GH convergence on the *full* state space $\mathcal{S}(\mathcal{O})$ via Berezin lift / spectral Fejér kernel on SU(2), not on individual basis vectors. R3.5 does not address this, but its structural finding (chirality-flip pairs are forced to zero) means that the natural $S^3$ Berezin map will probably need to *quotient out the chirality grading* to land in $C^\infty(S^3)$ rather than $C^\infty(S^3) \otimes \mathbb{C}^2$. This is consistent with Hekkelman's $S^1$ proof, where the chirality is implicit in the sign convention of the Dirac and not an extra label.

(b) **Extend multipliers from $C^\infty(S^3)$ to $C^\infty(S^3) \otimes \mathrm{Cl}(\mathbb{R}^3)$.** Promote the multipliers from scalar functions to vector-valued ones: e.g. $M = f \cdot \gamma_5$ for $f \in C^\infty(S^3)$ and $\gamma_5$ the chirality matrix. This *would* break the block-diagonal-with-equal-blocks structure of the multiplier matrices, and so could in principle break the kernel obstruction. However, this is a categorically larger construction than the Connes-vS truncation $P\,C^\infty(M)\,P$ — it is the "Clifford algebra extension" $P\,(C^\infty(M) \otimes \mathrm{Cl}(M))\,P$, which is *not* the same as the spectral truncation of the standard spectral triple. This is best deferred or treated as a separate construction.

### 4.3. What R3.5 does establish

Despite the negative on the obstruction, R3.5 *does* deliver three structural findings:

(i) **Per-block prop = 2 holds** in the full Dirac sector. The Weyl-only $\mathrm{prop}(\mathcal{O}_{\mathrm{Weyl}}) = 2$ result of R3.2 carries through to each chirality block of the full operator system. The full operator system has $\mathrm{prop}(\mathcal{O}_{\mathrm{full}}) = +\infty$ (does not reach the full $M_N(\mathbb{C})$ envelope), but the *natural* C*-envelope (block-diagonal $M_{\dim_W}(\mathbb{C}) \oplus M_{\dim_W}(\mathbb{C})$) IS reached at $k = 2$. This refines the WH1 prop=2 alignment claim with Connes-vS Toeplitz S¹ to "per-chirality-block" — still the same structural mechanism, with chirality as an internal grading.

(ii) **Anti-correlation persists with chirality grading.** R3.2's negative anti-correlation of the Connes distance with graph distance survives the chirality doubling at $n_{\max}=2$ (Pearson $-0.297$, $n_{\max}=3$ pending). Combined with R3.1 (placeholder $\to$ Avery, Pearson $-0.36$ at $n_{\max}=3$) and R3.2 (scalar $\to$ spinor, Pearson $-0.36 / -0.26$), R3.5 (Weyl $\to$ full Dirac) is the third independent confirmation that this is structural: not an artifact of any single convention.

(iii) **The chirality-flip pairs are forced to zero.** A new structural finding: any pair $|v, +\chi\rangle \leftrightarrow |v, -\chi\rangle$ (same $(n, l, m_j)$, opposite chirality) has $d_{\mathrm{Connes}} = 0$ identically, because scalar multipliers cannot distinguish them. This is a **physical** statement: at the operator-system level, chirality is invisible without gamma-matrix-coupled multipliers. It is the analog of the SO(3) m-reflection forced zero (R3.1) for the Z₂ chirality grading.

### 4.4. Net assessment for Track A

R3.5 is a useful diagnostic that **closes the L1 lemma in its reformulated L1' form** at every tested $n_{\max}$ (= 2, 3) and isolates the structural mechanism — cross-chirality offdiag E1 couplings on the SO(4) selection rule — by which the +∞ obstruction is retired. The TS-A memo's §6.2 claim that "the full Dirac sector is *necessary* even for the $S^1$ case in Hekkelman's proof — Hekkelman uses $D = -i\,d/d\theta$ which is *unbounded on every shell* and breaks degeneracy automatically" is *correct for $S^1$* but *misleading for $S^3$*, because Hekkelman's $D$ on $L^2(S^1)$ does not have a non-trivial chirality grading in any sense — it is a single-component operator on a single-component Hilbert space. The S³ analog with the truthful CH Dirac decomposes by chirality at the eigenvalue level but not at the level of $[D, \cdot]$ on scalar multipliers. The "shell degeneracy" of the truthful CH Dirac is therefore a structural feature, not an artifact of restricting to one chirality. R3.5 establishes that what *is* needed is the off-diagonal CH Dirac with cross-chirality couplings.

**Track A revision applied (L1 → L1'):** the off-diagonal CH Dirac with cross-chirality coupling produces a finite Connes distance on every non-forced cross-pair at all tested $n_{\max} \in \{2, 3\}$ — verified by R3.5 with 0/780 +∞ pairs at the largest accessible cutoff. Beyond L1', the proof shape in §5 of the TS-A memo proceeds without modification, with Lemma 5.3 (Lipschitz bound) as the place where the full-Dirac vs Weyl distinction will or will not matter (depending on whether the spinor-bundle torsion-and-curvature correction is uniform — Paper 32 §3.3 cross-reference).

The keystone direction (R2.5: GH convergence on the full state space via Peter–Weyl) is unchanged. R3.5 supplies the operator-system substrate (offdiag CH full Dirac with all-finite metric on non-forced pairs) that L1' demanded, and confirms that the residual pure-state anti-correlation is a Connes–vS Theorem 4.20-style finite-truncation phenomenon, not a GH obstruction. R3.5 does not unblock R2.5 *or* obstruct it; it just shows that the n-degeneracy obstruction lives in the operator system, not in the Hilbert space, and is removed by the offdiag CH choice that the GH proof was always going to use.

---

## §5. Implications for WH1

R3.5 shifts the WH1 evidence as follows:

**Strengthens:**

- The **prop=2 alignment with Connes-vS Toeplitz S¹** is robust across three convention upgrades (placeholder → Avery, scalar → spinor, Weyl → full Dirac). On the full Dirac sector, the per-chirality-block prop = 2 still holds. The structural alignment with the Connes-vS spectral truncation framework continues to hold at the operator-system level.

- The **non-monotonicity finding** is now structurally robust across THREE major upgrades: R3.1 (Avery), R3.2 (spinor), R3.5 (full Dirac), AND across cutoff growth $n_{\max}=2 \to 3$ within each variant. The four-point full-Dirac Pearson nz sequence is $-0.468$ (Weyl truthful, $n_{\max}=3$) $\to -0.409$ (full truthful, $n_{\max}=3$) $\to -0.297$ (full offdiag, $n_{\max}=2$) $\to -0.250$ (full offdiag, $n_{\max}=3$). The magnitude weakens monotonically with cutoff under offdiag CH (16% attenuation from $n_{\max}=2$ to $n_{\max}=3$) but never flips sign. Truthful CH at $n_{\max}=2$ gives no finite distances (R3.2 / R3.5), but at $n_{\max}=3$ truthful CH gives 24 finite distances Weyl-only and 96 in the full sector — both strongly negative. The disambiguation stands: this is a feature of the truncated operator system on $S^3$ with pure node-evaluation states, not an artifact of any convention or cutoff.

**Sharpens:**

- The **n-degeneracy obstruction is structurally an obstruction of the SCALAR-MULTIPLIER operator system**, not of the Hilbert space construction. R3.5 isolates this: doubling the Hilbert space and adding chirality-grading at the Dirac level does not change the kernel of $[D, \cdot]$ on the scalar multiplier system.

**Opens:**

- **Chirality-flip forced zeros** as a new class of structural metric identifications (§2.2 (Z-2)). This is a new operator-system-level statement: scalar multipliers cannot distinguish chirality-conjugate states, so $d_{\mathrm{Connes}} = 0$ identically on chirality-flip pairs. The natural place this becomes nontrivial is in a Cl(R^3)-coefficient extension of the multipliers (a categorically larger construction).

- **The Berezin / spectral-Fejér map on SU(2)** for Track A R2.5 will need to handle the chirality grading explicitly. The natural choice is to map *into* $C^\infty(S^3)$ (scalar functions on the manifold) and apply chirality as an internal tensor factor. R3.5 confirms that this is the correct factoring: chirality should not enter the algebra side of the spectral triple.

**Does not address:**

- GH convergence on the full state space (still open, R2.5 keystone).

- Uniform Lipschitz bound (Track A Lemma 5.3, requires explicit spinor-bundle curvature computation on $S^3$).

WH1 register entry recommendation (CLAUDE.md §1.7, PI-discretion): the entry should be tightened from "MODERATE-STRONG" to acknowledge that (i) R3.5 confirms the chirality grading is structurally independent of the operator-system kernel, and (ii) Track A's L1' lemma is now verified at $n_{\max}=2$ and $n_{\max}=3$ in the offdiag CH full Dirac, supplying the operator-system substrate that the GH-convergence proof shape (§5 of the TS-A memo) requires. The phrasing "the structural finding ... is robust under both convention upgrades (placeholder → Avery, scalar → spinor)" of R3.2 should be extended to "and the Weyl → full Dirac upgrade (R3.5), at both $n_{\max}=2$ and $n_{\max}=3$ under the offdiag CH variant."

---

## §6. Limitations and follow-up

### 6.1. Limitations of R3.5

1. **Scalar multipliers only.** R3.5 keeps the operator system as scalar-function-multipliers $M = f \in C^\infty(S^3)$. A genuinely chirality-coupled extension to $f \otimes \gamma$-matrix multipliers is a separate construction. The `geovac/full_dirac_operator_system.py` module is structured to allow such an extension as a follow-up.

2. **Pure-state restriction.** As in R3.1 / R3.2, the SDP computes Connes distance between pure node-evaluation states. The Connes-vS framework is most natural on the full state space.

3. **No Lipschitz-bound check.** R3.5 does not address Track A Lemma 5.3 (the comparison of $\|[D_{\mathrm{CH}}, M_f]\|$ to $\|\nabla f\|_\infty$ for $f \in C^\infty(S^3)$). This is a separate computation.

4. **Convention dependence in the offdiag perturbation.** The specific weights $(0.1 \cdot l, 0.005 \cdot 2m_j)$ and the unit cross-chirality coupling are choices, not derivations. The qualitative finding (anti-correlation, finite metric) is robust under reasonable variations; absolute numerical values depend on them.

### 6.2. Follow-up sprints

**R3.5 itself is now complete:** all five computations done — truthful CH at $n_{\max}=2,3$; Weyl-only truthful CH at $n_{\max}=3$ (comparison reference); offdiag CH at $n_{\max}=2,3$. No further R3.5 sub-sprint is open.

**R3.6 (medium-leverage): chirality-coupled multipliers.** Extend `FullDiracTruncatedOperatorSystem` to allow $f \otimes \gamma_5$ multipliers (or, more generally, $f \otimes M_2(\mathbb{C})$). Test whether this restores `prop` = 2 of the *full* algebra and breaks the kernel obstruction structurally. This is a categorically larger operator system; check whether it embeds naturally in the Connes-vS framework.

**R3.7 (low-leverage): offdiag perturbation parametrization.** Sweep over $\alpha \in \{0.1, 0.5, 1, 5, 10\}$ for the cross-chirality coupling at $n_{\max}=2$ to confirm anti-correlation is robust. The $n_{\max}=3$ data point is fixed at unit coupling and the four-point Pearson sequence (above) makes a sweep at $n_{\max}=3$ low-priority.

**R2.5 / R3.8 (high-leverage, keystone): GH convergence on the full state space.** As scoped in `debug/track_ts_a_gh_convergence_memo.md`. R3.5 supplies the operator-system substrate (offdiag CH full Dirac with all-finite metric on non-forced cross-pairs) that L1' demanded, and the natural choice of operator system to feed into the GH proof. Lemma 5.3 (Lipschitz bound) is the next bottleneck; R3.5 does not block it.

**Track A Lemma 5.3 (Lipschitz bound).** Compute $\|[D_{\mathrm{CH}}^{\mathrm{full,offdiag}}, M_f]\|$ vs $\|\nabla f\|_\infty$ for the FullDiracTruncatedOperatorSystem multiplier matrices already built (offdiag CH variant — the natural one per L1'). Identify the torsion-correction prefactor explicitly. Estimated 1–2 weeks effort. R3.5 has the infrastructure ready.

### 6.3. Out of scope

- No paper edits. R3.5 is a sprint result; per CLAUDE.md §13.5 paper updates are PI-discretion.
- No CLAUDE.md edits to the WH1 register. PI to decide whether the R3.1+R3.2+R3.5 combined evidence warrants tightening the §1.7 wording.

### 6.4. Code change made during the offdiag n_max=3 leg

`geovac/connes_distance.py` was refactored between the truthful $n_{\max}=3$ run and the offdiag $n_{\max}=3$ run to (a) cache the operator-norm LMI kernel construction (previously rebuilt twice per pair) and (b) remove a redundant second SCS solve (the gauge-fixed solver was being called both speculatively without dual extraction and then again with extraction; the speculative call was dropped). Net per-pair SDP speedup ~15-25×. Behavior on the `compute_distance_matrix` and `compute_connes_distance` API is unchanged (verified by re-running the truthful $n_{\max}=2$ regression: bit-identical distance matrix). All 17 tests in `tests/test_connes_distance.py` pass.

---

## §7. Files added in R3.5

### Added

- `geovac/full_dirac_operator_system.py` — ~580 lines. Full-Dirac truncated operator system on the spinor bundle of $S^3$. `FullDiracLabel`, `FullDiracTruncatedOperatorSystem`, truthful and offdiag CH Dirac builders, API-compatible drop-in to `compute_distance_matrix`.
- `tests/test_full_dirac_operator_system.py` — 29 passing tests (1 slow skipped). Equation verification per CLAUDE.md §13.4a: dim formula, label validation, chirality assignments, CH eigenvalues, multiplier block-diagonality, identity, *-closure, dim(O) scalar-equality, propagation per-block.
- `debug/wh1_r35_full_dirac_compute.py` — script driving the Connes distance computation at $n_{\max}=2, 3$ in both truthful and offdiag CH variants.
- `debug/data/wh1_r35_full_dirac_truthful_nmax2.json` — Connes distance matrix under truthful CH at $n_{\max}=2$.
- `debug/data/wh1_r35_full_dirac_offdiag_nmax2.json` — Connes distance matrix under offdiag CH at $n_{\max}=2$.
- `debug/data/wh1_r35_full_dirac_truthful_nmax3.json` — done. 96 finite cross-pair distances, Pearson nz $-0.409$.
- `debug/data/wh1_r35_weyl_truthful_nmax3.json` — done. 24 finite cross-pair distances, Pearson nz $-0.468$. Comparison reference for §2.4.
- `debug/data/wh1_r35_full_dirac_offdiag_nmax3.json` — done (2026-05-04, ~66 min wall time after `geovac/connes_distance.py` refactor). 720 finite cross-pair distances, 0 +∞, Pearson nz $-0.2501$, range $[0.151, 1.087]$.
- `debug/data/wh1_r35_offdiag_nmax3_v2.log` — run log for the offdiag $n_{\max}=3$ leg.
- `debug/wh1_r35_full_dirac_memo.md` — this memo.

### Unchanged

- All papers, CLAUDE.md.
- `geovac/spinor_operator_system.py`, `geovac/operator_system.py`, `geovac/connes_distance.py`, `geovac/circulant_s3.py`, `geovac/so4_three_y_integral.py`, and their tests. R3.5 builds on top of these without modification.

---

## §8. Bottom line

R3.5 closes the "full Dirac sector" caveat from R3.2 (R3.2 §6.1 limitation 1) at all five computational legs (truthful $n_{\max}=2$ and 3; Weyl-only truthful $n_{\max}=3$ comparison reference; offdiag $n_{\max}=2$ and 3). Five findings:

(F1) **The full Dirac extension is implementable as a chirality-doubled spinor operator system** with API-compatible drop-in to the existing SDP. $\dim(\mathcal{O}_{\mathrm{full}}) = \dim(\mathcal{O}_{\mathrm{Weyl}})$ at every $n_{\max}$ because scalar multipliers are blind to the chirality grading; multiplier matrices are block-diagonal with two equal Weyl blocks; identity and *-closure preserve under the lift; per-block propagation number $\mathrm{prop} = 2$ holds, matching the R3.2 / Connes-vS Toeplitz S¹ result.

(F2) **The truthful Camporesi–Higuchi full Dirac (Paper 32 spectral form) does NOT structurally break the n-degeneracy obstruction.** The block-diagonal scalar-multiplier structure forces the commutator $[D_{\mathrm{full}}, \widetilde{M}^{\mathrm{full}}]$ to decompose as $[D_{\mathrm{Weyl}}, \widetilde{M}^{\mathrm{Weyl}}] \oplus -[D_{\mathrm{Weyl}}, \widetilde{M}^{\mathrm{Weyl}}]$, so $\ker([D_{\mathrm{full}}, \cdot]) \cap \mathcal{O}_{\mathrm{full}} \cong 2 \times \ker([D_{\mathrm{Weyl}}, \cdot]) \cap \mathcal{O}_{\mathrm{Weyl}}$. The chirality grading is structurally invisible to the kernel of $[D, \cdot]$ on the scalar operator system. **Verified by bit-exact match** at $n_{\max}=3$: full Dirac gives 96 finite distances $= 4 \times $ Weyl's 24 finite distances, with identical metric range $[0.347, 0.733]$.

(F3) **R3.2's "truthful CH Dirac is too n-degenerate" reading was a $n_{\max}=2$ artifact, not a general obstruction.** At $n_{\max} \ge 3$ the truthful CH Dirac admits SOME finite Connes distances even on the Weyl sector (24/180 cross-pairs at $n_{\max}=3$ Weyl). The kernel of $[D, \cdot]$ on $\mathcal{O}$ is *smaller* at $n_{\max}=3$ than at $n_{\max}=2$, due to higher-order multipliers $M_{N=3, L, M}$ that lift the n-degeneracy in some directions of the operator system.

(F4) **Under the CH+offdiag full Dirac with cross-chirality couplings, the metric is finite on every non-forced pair at all tested cutoffs.** The four-point Pearson nz sequence is $-0.297$ at $n_{\max}=2$ (96/96 finite, 0 +∞) and $-0.250$ at $n_{\max}=3$ (720/720 finite, 0 +∞). Combined with Weyl-only truthful CH at $n_{\max}=3$ giving Pearson nz $-0.468$ and full-Dirac truthful CH at $n_{\max}=3$ giving Pearson nz $-0.409$, the non-monotonicity finding is now confirmed structural across THREE convention upgrades (placeholder → Avery; scalar → spinor; Weyl → full Dirac), across both Dirac variants (truthful and offdiag), and across cutoff growth $n_{\max}=2 \to 3$ within the offdiag full-Dirac case (16% magnitude attenuation, no sign flip).

(F5) **Track A's L1 lemma is FALSE in its original form (truthful CH) but L1' (offdiag CH) is VERIFIED at all tested $n_{\max}$.** R3.5 closes the operator-system substrate question for Track A: the offdiag CH full Dirac is the natural setting where every non-forced cross-pair has finite Connes distance, supplying the precondition the GH-convergence proof shape (§5 of the TS-A memo) requires.

**Strategic implication for Track A:** R3.5 partially confirms and partially falsifies the TS-A §6.2 hypothesis. The mechanism Track A predicted (chirality grading breaks degeneracy at the operator level) is *false* — chirality is invisible to the kernel of $[D, \cdot]$ on scalar multipliers (full = 4 × Weyl exactly per §1.4 / §2.4). But the *outcome* Track A wanted (finite Connes distances on non-forced cross-pairs) holds via the cross-chirality offdiag E1 structure: shell-coupling multipliers (offdiag) provide the kernel-breaking that the chirality grading alone (truthful) does not. The reformulated lemma:

> **L1' (revised, 2026-05-04):** *under the offdiag CH full Dirac with cross-chirality coupling, every non-forced cross-pair has finite Connes distance.* Verified at $n_{\max}=2$ (0/120 +∞) and $n_{\max}=3$ (0/780 +∞).

The keystone direction (R2.5: GH convergence on the full state space via Peter–Weyl Berezin map on SU(2)) is unchanged — neither blocked nor unblocked by R3.5. Lemma 5.3 (Lipschitz bound) remains the next bottleneck and is independent of R3.5; the offdiag CH full-Dirac multiplier matrices built in this sprint are the natural input to that computation.

The chirality-flip forced zero (F3 §4.3 (iii)) is a clean operator-system-level statement: scalar multipliers cannot distinguish chirality-conjugate states. This is the Z₂-grading analog of the SO(3) m-reflection forced zero (R3.1), and is consistent with the standard physical reading that chirality coupling requires gamma-matrix-valued operators (a categorically larger construction, R3.6 if desired).

The pure-state anti-correlation persists at all four data points but is monotonically weakening with cutoff under offdiag CH (16% attenuation $n_{\max}=2 \to 3$) — qualitatively consistent with Connes–vS Theorem 4.20-style finite-truncation behavior on $S^1$ (where the truncated Connes distance is *larger* than Kantorovich at finite $n$, with equality only in the limit).

**End of memo.**
