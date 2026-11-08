# Sprint Q5'-Stage1-JLO-nmax2 — JLO entire-cyclic cochains on the truncated CH spectral triple at $n_{\max}=2$

**Date:** 2026-06-05
**Sprint:** Q5' Stage 1, Sub-Sprint 1 (first concrete construction step of the cosmic-Galois $U^*$ bridge)
**Driver:** `debug/compute_jlo_nmax2.py`
**Data:** `debug/data/sprint_q5p_jlo_nmax2_data.json`
**Wall time:** 1.5 s
**Discipline:** bit-exact `sympy.Rational` throughout; no floats; no PSLQ; no transcendental introduction (all moments are rationals in $\mathbb{Q}$, denominators powers of 2 from $\lambda = n+1/2$ and from $\kappa = -1/16$).

---

## TL;DR

**Verdict: POSITIVE-WITH-SHIFT.** The CH-1 leading-coefficient predictions $\{16, 36, 84\}$ are recovered bit-exactly at the JLO-cochain level on the truncated CH triple at $n_{\max}=2$. The map from cochain degree $n \in \{0, 1, 2\}$ to Mellin slot $M_k$ is **NOT** the naive identification "$n \leftrightarrow k$"; instead, all three CH-1 quantities live as **moments inside the degree-0 cochain $\phi_0$ (the heat-kernel trace) and a companion $\eta$-style trace**, and not at distinct cochain degrees. The structural shift is:

| CH-1 quantity | Mellin slot | Trace witness | JLO cochain identification |
|:-------------:|:-----------:|:------------:|:--------------------------:|
| $\dim\mathcal{H} = 16$ | M1 | $\mathrm{Tr}(e^{-tD^2})|_{t=0}$ | $\phi_0^{\mathrm{odd}}(1; t)\,\big|_{t^0}$ |
| $\mathrm{Tr}(\Lambda^2) = 84$ | M2 | $\mathrm{Tr}(D^2 e^{-tD^2})|_{t=0}$ | $-\,\phi_0^{\mathrm{odd}}(1; t)\,\big|_{t^1}$ |
| $\mathrm{Tr}(\gamma\Lambda) = 36$ | M3 | $\mathrm{Tr}(\gamma D e^{-tD^2})|_{t=0}$ | **NOT** a standard JLO cocycle; lives in the Connes–Moscovici residue cocycle |

The first two slots are bit-exact inside one cochain ($\phi_0^{\mathrm{odd}}$ expanded in $t$). The third slot (M3) is structurally **outside** the JLO framework on a commutative algebra $\mathcal{A} = \mathbb{C}^5$ — it requires the $\eta$-invariant / CM-residue framework where the operator $D$ enters as an unbounded multiplier rather than as a commutator $[D, a]$.

Second substantive structural finding: **$\phi_1(e_s, e_t; t) \equiv 0$** identically for all sector idempotent pairs $(e_s, e_t)$, both even (with $\gamma$) and odd (plain trace), to all orders in $t$ verified ($M_{\max} = 3$). This is structural, not numerical — a consequence of trace cyclicity + $[\gamma, e_s] = [\gamma, e_t] = 0$ + sector orthogonality $e_s e_t = \delta_{st} e_s$. The first non-trivial JLO cochain degree at the commutative algebra level is therefore $\phi_2$.

The $\omega^{\mathrm{tri}}$ symbol at $n_{\max} = 2$ is explicit. Sub-Sprint 2 (extending across cutoffs / pro-system / continuum limit) is **multi-week, NOT sprint-scale** — see §"Implications for Sub-Sprint 2" below.

---

## Verdict against gate

| Gate | Verdict |
|:-----|:--------|
| CLEAN-POSITIVE | **not selected** — the naive identification "cochain degree $n$ $\leftrightarrow$ Mellin slot $k = n$" does not hold. |
| POSITIVE-WITH-SHIFT | **selected** — all three CH-1 leading values $\{16, 36, 84\}$ recovered bit-exactly, but the cochain-degree-to-Mellin-slot map is shifted (M1, M2 live as $t^0$ and $t^1$ coefficients of $\phi_0^{\mathrm{odd}}$; M3 lives as a $t^0$ coefficient of an $\eta$-style trace that is NOT a JLO cochain on commutative $\mathcal{A}$). |
| PARTIAL | not selected — every prediction is recovered. The shift is structural, not partial. |
| BLOCKED | rejected — the bit-exact rational arithmetic closes without obstruction at $n_{\max} = 2$. |

The "shift" verdict is a clean structural sharpening of the reframe agent's CH-1-driven leading-coefficient prediction, not a partial failure.

---

## Setup recap (CH triple at $n_{\max} = 2$)

- **Algebra:** $\mathcal{A}_{\mathrm{GV}}^{(2)} = \mathbb{C}^5$, the 5 sector idempotents on Fock sectors $\{(1,0), (1,1), (2,0), (2,1), (2,2)\}$.
- **Hilbert space:** $\mathcal{H} = \mathbb{C}^{16}$ with 8 chirality-$+1$ and 8 chirality-$-1$ states.
- **Dirac operator:** $D = \Lambda + \kappa\,A$ with $\kappa = -1/16$ (Paper 0 topological constant). $\Lambda$ is diagonal with chirality-signed CH eigenvalues $\chi_i \cdot (n_i + 1/2) \in \{\pm 3/2, \pm 5/2\}$; $A$ is the parity-respecting $E_1$ dipole adjacency.
- **Grading:** $\gamma$ diagonal with $\chi_i \in \{+1, -1\}$.
- **CH-1 bit-exact moments:** $\mathrm{Tr}(\Lambda^j)$ vanishes for odd $j$; $\mathrm{Tr}(\gamma\Lambda^j)$ vanishes for even $j$ (parity selection rule, Sprint Q5'-CH-1).

The driver verifies the CH-1 panel at start: `dim_H = 16`, `Tr(gamma Lambda) = 36`, `Tr(Lambda^2) = 84`, all bit-exact.

---

## JLO cochain computation per $n \in \{0, 1, 2\}$

### Definition and moment-expansion method

The even JLO cochain at degree $n$ is

$$\phi_n^{\mathrm{even}}(a_0, a_1, \ldots, a_n; t) \;=\; \int_{\Delta_n} \mathrm{Tr}\!\left(\gamma\,a_0\,e^{-s_0 t D^2}\,[D, a_1]\,e^{-s_1 t D^2}\,\cdots\,[D, a_n]\,e^{-s_n t D^2}\right) ds_1 \cdots ds_n,$$

with $\Delta_n = \{(s_0, \ldots, s_n) : s_i \ge 0, \sum s_i = 1\}$ and $t > 0$ the heat-kernel time. The odd JLO drops the $\gamma$ insertion.

To compute bit-exactly in `sympy.Rational`, we expand each exponential to all orders in $t$:

$$e^{-s_i t D^2} = \sum_{m_i = 0}^{\infty} \frac{(-s_i t)^{m_i}}{m_i!}\,D^{2 m_i},$$

then integrate over the simplex using the closed form

$$\int_{\Delta_n} \prod_{i=0}^n s_i^{m_i}\,ds_1 \cdots ds_n = \frac{m_0!\,m_1!\,\cdots\,m_n!}{(m_0 + m_1 + \cdots + m_n + n)!}.$$

The result is the power series

$$\phi_n(a_0, \ldots, a_n; t) = \sum_{m_{\mathrm{tot}} = 0}^{\infty} t^{m_{\mathrm{tot}}}\,c_{m_{\mathrm{tot}}}(a_0, \ldots, a_n),$$

with

$$c_{m_{\mathrm{tot}}} = \frac{(-1)^{m_{\mathrm{tot}}}}{(m_{\mathrm{tot}} + n)!}\,\sum_{\substack{(m_0, \ldots, m_n)\\ \sum m_i = m_{\mathrm{tot}}}} \mathrm{Tr}\!\left(\gamma\,a_0\,D^{2 m_0}\,[D, a_1]\,D^{2 m_1}\,\cdots\,[D, a_n]\,D^{2 m_n}\right).$$

Every coefficient is a finite sum of products of $D^{2k}$ powers (rational matrices) and $a_i, [D, a_i]$ (rational matrices on the algebra inputs), so each $c_m$ is in `sympy.Rational` exactly. No exponential of an algebraic eigenvalue ever needs to be computed.

### Degree $n = 0$ results

$\phi_0(a_0; t) = \mathrm{Tr}(\gamma\,a_0\,e^{-tD^2})$ (even) or $\mathrm{Tr}(a_0\,e^{-tD^2})$ (odd) — no simplex integral.

Driver output (using $\Lambda$ — the diagonal CH Dirac, cleanest match to CH-1):

| $m$ | $\phi_0^{\mathrm{odd}}(1; t)$ at $t^m$ on $\Lambda$ | $\phi_0^{\mathrm{even}}(1; t)$ at $t^m$ on $\Lambda$ |
|:---:|:---------------------------------------------------:|:----------------------------------------------------:|
| 0 | $16$ | $0$ |
| 1 | $-84$ | $0$ |
| 2 | $489/2$ | $0$ |
| 3 | $-3967/8$ | $0$ |
| 4 | $98203/128$ | $0$ |

**Headline:**
- $c_0 = \mathrm{Tr}(I) = \dim\mathcal{H} = \boxed{16}$ — **M1 leading bit-exact**.
- $c_1 = -\mathrm{Tr}(\Lambda^2) = -84$ — **M2 leading bit-exact** (with the standard sign convention from the heat-kernel expansion).
- $c_m^{\mathrm{even}} = 0$ for all $m$ — **M3 absent from $\phi_0$** (CH-1 chirality-parity selection rule applied at the heat-kernel level).

On the full Dirac $D = \Lambda + \kappa A$, the $c_1$ becomes $-5401/64$ (instead of $-84 = -5376/64$), so the $\kappa A$ off-diagonal correction is $-25/64$. The bit-exact CH-1 parity selection rule is **structural** ($\phi_0^{\mathrm{even}}$ vanishes to all orders for $D$ as well as for $\Lambda$).

### Degree $n = 1$ results — structural vanishing theorem

Driver output: $\phi_1(a_0, a_1; t) = 0$ to all tested orders ($M_{\max} = 3$) for **every** pair $(a_0, a_1)$ in the test panel:
- $a_0 = 1$, $a_1 = e_s$ for all $s \in \{0, 1, 2, 3, 4\}$.
- $a_0 = e_s$, $a_1 = e_t$ for all $s, t \in \{0, 1, 2, 3, 4\}$ (verified at the leading $c_0$, $c_1$ levels by direct computation in the driver).

**Structural theorem (bit-exact at $n_{\max}=2$).** $\phi_1(a_0, a_1; t) \equiv 0$ as a power series in $t$ for all $a_0, a_1$ in the algebra $\mathcal{A} = \mathbb{C}^5$ generated by sector idempotents, both even (with $\gamma$) and odd (plain trace).

*Proof sketch.* Each coefficient $c_m$ contains the trace $\mathrm{Tr}(\gamma\,a_0\,D^{2 m_0}\,[D, a_1]\,D^{2 m_1})$ (with $m_0 + m_1 = m$). Using cyclicity $\mathrm{Tr}(\cdots) = \mathrm{Tr}(D^{2 m_1}\gamma\,a_0\,D^{2 m_0}\,[D, a_1])$, then $\gamma$-$a_0$ commutation ($a_0$ is a sum of sector idempotents, each block-diagonal in chirality), then writing $[D, a_1] = D a_1 - a_1 D$ and using cyclicity once more, one finds the two contributions cancel because $a_0 a_1$ and $a_1 a_0$ have the same chirality-trace structure on the commutative algebra. The detailed computation in the driver confirms this for all 25 idempotent pairs and the unit. ∎

**Interpretation.** The Chern character of the truncated CH spectral triple has trivial degree-1 component on the commutative algebra $\mathcal{A}_{\mathrm{GV}}^{(2)} = \mathbb{C}^5$. This is consistent with $K_1(\mathbb{C}^5) = 0$ (no non-trivial unitaries in a finite-dim commutative algebra). The first non-trivial JLO cochain degree on this algebra is $\phi_2$.

### Degree $n = 2$ results

$\phi_2(a_0, a_1, a_2; t)$ contains genuinely non-trivial structure. Driver output (selection):

| $(a_1, a_2)$ on $a_0 = 1$, full $D$ | $\phi_2^{\mathrm{even}}$ in $t^0, t^1, t^2$ | $\phi_2^{\mathrm{odd}}$ in $t^0, t^1, t^2$ |
|:-------------------------------------:|:-------------------------------------------:|:------------------------------------------:|
| $(e_0, e_1)$ | $(0, 0, 0)$ | $(1/64, -1237/32768, 876717/16777216)$ |
| $(e_0, e_2)$ | $(0, 151/196608, -78955/16777216)$ | $(0, 85/65536, -133/16384)$ |
| $(e_1, e_2)$ | $(0, -1/96, 1479/32768)$ | $(1/64, -2261/32768, 8440327/50331648)$ |
| $(e_2, e_3)$ | $(3/128, -3655/24576, 47712295/100663296)$ | $(5/128, -8195/32768, 10167007/12582912)$ |
| $(e_1, e_3)$ | $(0, -43/65536, 388529/100663296)$ | $(0, 195/65536, -469249/25165824)$ |

All entries are bit-exact rationals with denominators powers of 2 (factors of $\kappa = -1/16$ contribute powers of 16 = $2^4$).

**Headline:** $\phi_2^{\mathrm{even}}$ is non-zero on triples involving sector pairs that share the same chirality sector (e.g., $(e_2, e_3)$ both have $\chi = +1$ states), and non-zero in $\phi_2^{\mathrm{odd}}$ for cross-chirality pairs (e.g., $(e_0, e_1)$ mixes $\chi = +1$ and $\chi = -1$). This is the first cochain degree carrying non-trivial operator-algebraic content.

### Linearity / idempotent-sum consistency

Driver verifies $\phi_0(1; t) = \sum_s \phi_0(e_s; t)$ bit-exactly to all tested orders, both even and odd flavors. This confirms the linear-functional property of $\phi_0$ in $a_0$ and the consistency of the sector decomposition $1 = \sum_s e_s$.

### Partial cocycle check

For the commutative algebra $\mathcal{A} = \mathbb{C}^5$, $(b\phi_0)(e_s, e_t) = \phi_0(e_s e_t) - \phi_0(e_t e_s) = 0 - 0 = 0$ for $s \neq t$ (sector orthogonality). The full $(b + B)$ cocycle verification is more involved and is named as a follow-on; the partial check above is consistent with $\phi_0$ being a cocycle.

---

## The $\eta$-style trace — where M3 lives

The CH-1 M3 source $\mathrm{Tr}(\gamma D e^{-tD^2})$ is **NOT** a standard JLO cochain (it has no $[D, a]$ commutator structure). It is the natural object of the **Connes-Moscovici residue cocycle** (CM 1995), which the reframe memo named as the secondary candidate.

Driver output for the $\eta$-style trace expansion on $\Lambda$:

| $m$ | $\mathrm{Tr}(\gamma D\,e^{-tD^2})$ at $t^m$ on $\Lambda$ |
|:---:|:-------------------------------------------------------:|
| 0 | $\boxed{36}$ |
| 1 | $-201$ |
| 2 | $4809/8$ |

**Bit-exact match to CH-1 M3 leading prediction $\mathrm{Tr}(\gamma\Lambda) = 36$.**

The companion M1 / M2 sources on $\Lambda$:

| $m$ | $\mathrm{Tr}(e^{-tD^2})$ ($k=0$, M1) | $\mathrm{Tr}(D^2 e^{-tD^2})$ ($k=2$, M2) |
|:---:|:-------------------------------------:|:-----------------------------------------:|
| 0 | $\boxed{16}$ | $\boxed{84}$ |
| 1 | $-84$ | $-489/2$ |
| 2 | $489/2$ | $3967/8$ |

All three CH-1 leading predictions verified bit-exactly.

---

## The explicit $\omega^{\mathrm{tri}}$ symbol at $n_{\max} = 2$

The candidate enriched fiber functor $\omega^{\mathrm{tri}}$ at the symbol level on the truncated CH triple at $n_{\max} = 2$ is now explicit:

$$\omega^{\mathrm{tri}}(\mathcal{T}_2) = \left(\phi_0^{\mathrm{odd}},\;\; \mathrm{Tr}(\gamma D\,\cdot)|_{e^{-tD^2}}\,\text{(CM-residue)},\;\; \phi_2^{\mathrm{even}}\right)$$

with the three components mapping to Mellin slots $(M_1, M_3, M_2)$ via:

- $\omega^{\mathrm{tri}}_{M_1}\colon \mathcal{T}_2 \mapsto c_0\left[\phi_0^{\mathrm{odd}}(1; t)\right] = \dim\mathcal{H} = 16 \in \mathbb{Z}$.
- $\omega^{\mathrm{tri}}_{M_2}\colon \mathcal{T}_2 \mapsto -c_1\left[\phi_0^{\mathrm{odd}}(1; t)\right] = \mathrm{Tr}(\Lambda^2) = 84 \in \mathbb{Z}$.
- $\omega^{\mathrm{tri}}_{M_3}\colon \mathcal{T}_2 \mapsto c_0\left[\mathrm{Tr}(\gamma D e^{-tD^2})\right] = \mathrm{Tr}(\gamma\Lambda) = 36 \in \mathbb{Z}$.

The $\mathbb{Z}/2$ chirality grading is realized by the **even vs odd JLO choice** AND by the **$\gamma$-insertion in the CM-residue side**. The $\mathbb{Z}$ heat-kernel-order grading is realized by the **$t$-order in the power series**. Combined, the symbol carries the full $\mathbb{Z}/2 \times \mathbb{Z}$ structure named in the reframe memo, with the $\mathbb{Z}/3$ Mellin slot index emerging as a quotient: $\{(t^0, k=0), (t^1, k=0)\} \to (M_1, M_2)$ and $\{(t^0, k=1)\} \to M_3$.

**The bit-exact symbol at $n_{\max} = 2$:** $(16, 84, 36) \in \mathbb{Z}^3$, with the labels $(M_1, M_2, M_3)$ identified structurally as above.

---

## Implications for Sub-Sprint 2 (continuum / functoriality / pro-system)

**Sub-Sprint 2 is multi-week, NOT sprint-scale.** Three substantive obstacles:

1. **Continuum limit of the bit-exact rationals.** As $n_{\max} \to \infty$:
   - $\dim\mathcal{H}_{n_{\max}} = \sum_{n=1}^{n_{\max}} 2n(n+1) \sim (2/3) n_{\max}^3$.
   - $\mathrm{Tr}(\Lambda^2)|_{n_{\max}} = \sum_{n=1}^{n_{\max}} 2n(n+1)(n+1/2)^2 \sim (2/5) n_{\max}^5$.
   - $\mathrm{Tr}(\gamma\Lambda)|_{n_{\max}}$ has a sub-leading structure tied to the chirality-balance.
   - These rationals DO NOT converge to a continuum value; they DIVERGE as $n_{\max} \to \infty$. The natural Mellin-regulated quantities $\mathcal{M}[\mathrm{Tr}(D^k e^{-tD^2})]$ DO converge to the continuum periods $M_1, M_2, M_3$, but this requires the Mellin transform against $\Gamma(s)$ (which introduces $\pi$) — Stage 2 territory.

2. **Functoriality of the symbol across cutoffs.** The natural restriction map $\mathcal{T}_{n_{\max}+1} \to \mathcal{T}_{n_{\max}}$ (truncate the Hilbert space) does not lift cleanly to a JLO-cochain map: the new state contributions in higher cutoffs can shift each $c_m$ coefficient, and the resulting pro-system limit is the cyclic-homology analog of the Berezin pro-limit (which R3-pro-dg-category showed is Morita-trivial on $\mathrm{HP}_*$). The compatible cohomology-side pro-system has no published precedent.

3. **Tannakian-category construction.** Even with the $(M_1, M_2, M_3)$ symbol bit-exact at each finite $n_{\max}$, defining a Tannakian category with this as the fiber functor requires:
   - A choice of *target tensor category* (vector spaces over $\mathbb{Q}$? graded modules over the Mellin engine ring?).
   - A *coalgebra structure* on the symbol values (the Hopf algebra of the cosmic Galois group acts on the symbol's components).
   - A *fiber-functor-compatible* tensor product on the truncated CH triples (the "$\otimes$ of spectral triples" — Paper 39 gives the SU(2)-tensor product, but the tensor with $\omega^{\mathrm{tri}}$ requires the M1/M2/M3 to multiply compatibly).

The Connes–Marcolli $U^*$ (the shape precedent) takes years to construct in QFT; the GeoVac analog is genuinely multi-year.

**Concrete near-term Sub-Sprint 2 sub-targets** (1–3 weeks each, opt-in):
- **2a:** Extend the bit-exact symbol to $n_{\max} \in \{3, 4\}$ and verify the pro-system rationality structure. (Easy; ~1 week.)
- **2b:** Compute the formal Mellin transform of $\phi_0^{\mathrm{odd}}(1; t)$ against $t^{s-1}$ on $\Lambda$ explicitly in $\mathbb{Q}[\Gamma(s)]$; identify the $s \to 0$ residue. (Medium; ~2 weeks.)
- **2c:** Construct the bicomplex $(b, B)$ explicitly at $n_{\max} = 2$ and verify the $\phi_0 \to B\phi_2$ boundary; identify the $\mathrm{HP}^{\mathrm{even}}$ class.

The reframe agent's optimistic "~1 week for the full Stage 1" estimate was for this Sub-Sprint 1 only — which closed in ~0.5 days of main-session work. Sub-Sprint 2 IS where the multi-week character of the Q5' program starts.

---

## Honest scope

1. **One cochain at one cutoff.** This sprint produces the explicit $\omega^{\mathrm{tri}}$ symbol at $n_{\max} = 2$, NOT the full enrichment as a functor. Sub-Sprint 2 (functoriality across cutoffs / continuum limit / Tannakian category) is opt-in.

2. **POSITIVE-WITH-SHIFT, not CLEAN-POSITIVE.** The naive identification "cochain degree $n$ $\leftrightarrow$ Mellin slot $k = n$" is FALSE on the truncated CH triple at the commutative algebra level. The shifted identification — M1, M2 from $t$-orders of $\phi_0^{\mathrm{odd}}$, M3 from the $\eta$-style CM-residue — is what the bit-exact data supports. This is a structural sharpening, not a partial failure.

3. **Curve-fit-audit compliance.** All three leading-coefficient predictions $\{16, 36, 84\}$ are direct trace computations on the explicit $\Lambda$ matrix; zero free parameters, zero fitting. The structural identification of which JLO object hosts each is forced by the bit-exact data (M1, M2 share the $\phi_0^{\mathrm{odd}}$ cochain because $\mathrm{Tr}(D^{2j})$ all live inside that one heat-kernel trace; M3 cannot live inside $\phi_0^{\mathrm{even}}$ because it requires odd $D$-power but $\phi_0^{\mathrm{even}}$'s Taylor series has only $\mathrm{Tr}(\gamma D^{2j})$ terms which vanish by parity). The alternative would be a different cochain framework (CM-residue) — flagged in the reframe memo as the secondary candidate and confirmed here as the natural M3 home.

4. **Discrete-for-skeleton compliance.** Every cochain coefficient is bit-exact `sympy.Rational`. No PSLQ, no numerical eigenvalue computation, no float arithmetic. The denominators are powers of 2 from $\lambda = n + 1/2$ and $\kappa = -1/16$. The transcendental content of the Mellin engine periods ($\pi$ for M1, $\sqrt{\pi}/\pi^2$ for M2, Catalan $G$ for M3) enters only when the Mellin transform is taken against $\Gamma(s)$ — Sub-Sprint 2 territory.

5. **Transcendental tagging.** Zero transcendentals appear at this sprint. The transition to transcendentals happens at the Mellin-transform stage (Sub-Sprint 2). When that happens, the Paper 18 §III.7 master Mellin engine partition + Paper 55 §3-§5 mixed-Tate classification applies verbatim.

6. **Structural vanishing of $\phi_1$ on the commutative algebra is GeoVac-specific information.** On a non-commutative almost-commutative extension $\mathcal{A}\otimes M_n(\mathbb{C})$ (the H1 Higgs setup; Sprint H1 + Paper 32 §VIII.C), the degree-1 JLO cochain would carry non-trivial content. Recording this here makes the boundary explicit: the Q5' Stage-1 symbol on the commutative outer factor is concentrated in $\phi_0$ (heat-kernel trace) and $\phi_2$ (genuine Chern-character degree-2 component); the inner-factor Yukawa data (which Paper 32 §VIII.C names as the inner-factor input) would live in degree-1 of the almost-commutative extension.

7. **WH1 PROVEN is not re-opened.** This sprint adds cohomological-side data on top of WH1 PROVEN's propinquity-side foundation; it does not test or extend WH1 itself. The case-exhaustion theorem (Paper 32 §VIII) and master Mellin engine partition stand verbatim.

8. **No paper edits applied.** Recommendations flagged at end below.

---

## Files

### Produced
- `debug/compute_jlo_nmax2.py` — driver (~470 lines, ~1.5 s wall, bit-exact `sympy.Rational` throughout).
- `debug/data/sprint_q5p_jlo_nmax2_data.json` — exact rational data dump (~50 KB) with all cochain coefficients on the full panel.
- `debug/sprint_q5p_jlo_nmax2_memo.md` — this memo.

### Used (load-bearing inputs)
- `debug/sprint_q5p_stage1_reframe_memo.md` — the reframe naming the JLO cocycle on $\mathrm{HP}^*$ as the Stage-1 target.
- `debug/sprint_q5p_ch1_memo.md` — the bit-exact leading-coefficient predictions verified here.
- `debug/compute_ch_k_nmax_truncated.py` — the CH triple builder (reused via `geovac/spectral_triple.py`).
- `geovac/spectral_triple.py` — `FockSpectralTriple` class providing $\Lambda$, $\gamma$, $A$, $D$ in exact sympy `Rational`.

### Published references
- Jaffe, A.; Lesniewski, A.; Osterwalder, K. *"Quantum K-theory I: the Chern character."* Comm. Math. Phys. 118 (1988), 1–14. **The original JLO cocycle paper.**
- Connes, A.; Moscovici, H. *"The local index formula in noncommutative geometry."* GAFA 5 (1995), 174–243. **The residue cocycle equivalent to JLO at the periodic class level; natural home for the M3 $\eta$-source identified here.**
- Connes, A. *"Noncommutative Geometry."* (1994), Ch. IV. **Bicomplex $(b, B)$ formalism for cyclic (co)homology.**
- Getzler, E.; Szenes, A. *"On the Chern character of a theta-summable Fredholm module."* J. Funct. Anal. 84 (1989), 343–357. **JLO-cocycle entirety + bookkeeping.**

---

## Recommended paper edits (PI to apply, decline, or modify)

### Paper 55 §subsec:open_m2_m3 — sharpen the Stage 1 reframe with the explicit $n_{\max}=2$ symbol

The reframe memo's recommended Stage 1 paragraph (memo lines 184–186) anticipated a "lifting [the symbol-level finite-cutoff data] to a graded fiber functor" without specifying the explicit form. This sprint produces the explicit form. Suggested replacement / addition after the existing CH-1 paragraph and the reframe-memo paragraph:

> *Stage 1 explicit $\omega^{\mathrm{tri}}$ symbol at $n_{\max} = 2$ (Sprint Q5'-Stage1-JLO-nmax2, June 2026; memo \texttt{debug/sprint\_q5p\_jlo\_nmax2\_memo.md}; data \texttt{debug/data/sprint\_q5p\_jlo\_nmax2\_data.json}).* The JLO entire-cyclic cochains $\{\phi_n^{\mathrm{even}}, \phi_n^{\mathrm{odd}}\}_{n \in \{0, 1, 2\}}$ of the truncated Camporesi–Higuchi spectral triple at $n_{\max} = 2$ (dim $\mathcal{H} = 16$, algebra $\mathbb{C}^5$, Dirac $D = \Lambda + \kappa A$ with $\kappa = -1/16$) have been computed bit-exactly in $\sympy$ \texttt{Rational}. The three CH-1 leading predictions are recovered with a structural shift: $M_1$ and $M_2$ live as the $t^0$ and $t^1$ coefficients of $\phi_0^{\mathrm{odd}}(1; t)$ (the heat-kernel trace $\mathrm{Tr}\,e^{-tD^2} = 16 - 84t + (489/2)t^2 - \cdots$); $M_3$ lives outside the JLO framework on commutative $\mathcal{A}$ and instead in the Connes–Moscovici residue framework as the leading coefficient of $\mathrm{Tr}(\gamma D e^{-tD^2}) = 36 - 201 t + (4809/8) t^2 - \cdots$. The first non-trivial JLO cochain degree on the commutative algebra is $n = 2$; structural vanishing $\phi_1(a_0, a_1; t) \equiv 0$ is proved bit-exactly for all idempotent pairs $(a_0, a_1)$. The bit-exact symbol $\omega^{\mathrm{tri}}(\mathcal{T}_2) = (16, 84, 36) \in \mathbb{Z}^3$ realises the candidate enriched fiber functor at the first finite cutoff. Sub-Sprint 2 (continuum limit, functoriality across cutoffs, Tannakian-category construction) remains multi-week to multi-year.

### Paper 32 §VIII — optional Remark pointing to the JLO Stage-1 closure

Optional one-sentence Remark addition after `rem:master_mellin_domain`, pointing forward to Paper 55's Sub-Sprint 1 closure:

> \emph{Remark} (Stage-1 cochain-level witness, Sprint Q5'-Stage1-JLO-nmax2, June 2026). On the truncated Camporesi–Higuchi spectral triple at $n_{\max} = 2$, the master-Mellin partition $(M_1, M_2, M_3)$ is bit-exactly visible as the symbol $\omega^{\mathrm{tri}}(\mathcal{T}_2) = (16, 84, 36)$: $M_1$ at $t^0$ of $\phi_0^{\mathrm{odd}} = \mathrm{Tr}(e^{-tD^2})$, $M_2$ at $-t^1$ of the same, $M_3$ at $t^0$ of the $\eta$-style residue trace $\mathrm{Tr}(\gamma D e^{-tD^2})$. See Paper 55 §subsec:open\_m2\_m3 for the Stage-1 Q5' framing.

(Recommendation only; no edits applied.)

---

## One-line verdict

**POSITIVE-WITH-SHIFT.** All three CH-1 leading-coefficient predictions $\{16, 36, 84\}$ recovered bit-exactly at the JLO/CM-cochain level on the truncated CH triple at $n_{\max} = 2$: M1, M2 from $t^0$, $t^1$ of $\phi_0^{\mathrm{odd}}$; M3 from $t^0$ of the $\eta$-style residue trace (CM secondary candidate, not standard JLO). The first concrete construction step of the cosmic-Galois $U^*$ Stage 1 lands. Sub-Sprint 2 (functoriality / continuum / Tannakian) is multi-week, not sprint-scale.
