# Prior art: symmetries, gauge groups and order-one defects of spectral truncations

Date 2026-10-31. Read-only web scan for WG (`docs/predictions/2026-10-31_window_gauge_group.md`) and the order-one edge-crossing formula. The items already credited in `gauge_invariant_truncation_memo.md` are not re-reported: Choi, Dong 2022, Connes–van Suijlekom 2021, Verhoeven 2026, Marcolli–van Suijlekom, and CCvS. The one exception is Verhoeven, re-checked below for symmetry content. Nothing was downloaded and read. One PDF was briefly copied into `dl_prior_truncation/` and deleted unread. Everything below comes from arXiv abstract pages and arXiv/ar5iv HTML full text.

## Verdict

- **Order-one failure on a Dirac spectral truncation: ANTICIPATED, qualitatively, on S².** Glaser and Stern measured it numerically on the spinor truncation of the round 2-sphere and called it a "boundary defect".
- **Explicit edge-crossing formula on S³: not found.** It is the spectral-window analogue of the classical Widom Toeplitz–Hankel identity, and should be credited as such.
- **Dong's gauge group G(E) on any truncation beyond the circle (torus, sphere, SU(2)/S³): not found.** WG extends Dong. Verhoeven's SU(2) paper treats no symmetries.

## Table

| Citation | What it shows | Relation | Verified |
|---|---|---|---|
| L. Glaser, A. B. Stern, "Understanding truncated non-commutative geometries through computer simulations", J. Math. Phys. 61, 033507 (2020), arXiv:1909.08054 | Truncation (P_Λ A P_Λ, P_Λ H, P_Λ D) with P_Λ = χ_[−Λ,Λ](D) on L²(S², S). For the round sphere the double commutator [[D, Y_Λ], Y_Λ], built from coordinate functions, is "a boundary defect of asymptotically constant norm" with trace norm O(Λ) (§3.2). The defect masks the first-order failure of perturbed Dirac operators (§3.1). | **ANTICIPATES** (qualitatively: order one fails on a spinor Dirac window, located at the cut). Our result adds the exact S³ formula. Note the growth: operator norm O(1), trace norm O(Λ). Check which norm Paper 32's "grows" refers to. | full text (ar5iv §3.1–3.2); journal ref from the abstract page |
| L. Glaser, A. B. Stern, "Reconstructing manifolds from truncated spectral triples", arXiv:1912.09227 | §6.3: failure of the first order condition "is not detectable by standard methods at the level of truncated spectral triples" (attributed to their companion paper). | RELATED (supports reading the defect as intrinsic to the cut) | full text (ar5iv §6.3); journal ref not checked |
| Widom identity T(ab) = T(a)T(b) + H(a)H(b̃) (Widom; also Basor–Ehrhardt, Böttcher; e.g. arXiv:math/0012200) | The failure of compression to be multiplicative is exactly a Hankel product P a (1−P) b P. | **SUPPORTS / credit**: our edge-crossing sum is the Dirac-window version, PaP·PbP − PabP = −P a Q b P with Q = 1−P, sandwiched with D. | secondary (search snippets of Böttcher's and Basor's texts; Widom 1976 primary not read) |
| D. Farenick, "The operator system of Toeplitz matrices", Trans. AMS Ser. B 8 (2021), arXiv:2103.16546 | The Connes–van Suijlekom identification is a unital complete order isomorphism. A linear unital isometry of Toep_n into M_n is a unitary similarity. | SUPPORTS (complete-isometry classification on the circle; matches Dong) | abstract; journal from AMS page title |
| D. Farenick, M. Mastnak, A. Popov, "Isometries of the Toeplitz matrix algebra", arXiv:1502.01573 | Isometries of upper-triangular Toeplitz algebras have the form A ↦ UAU* or UĀU*. | RELATED (the triangular algebra, not the operator system) | abstract |
| F. Latrémolière, "How to approximate the flat spectral triple of a quantum torus by fuzzy tori: a twisted tale", arXiv:2607.01681 (2026) | Discrete calculus is non-local, so "the Leibniz rule fails". Spectral truncations avoid this by "abandoning C*-algebras in favor of operator systems". The author instead uses twisted triples whose twist is a discretized Riesz transform. | RELATED (an alternative repair of the Leibniz defect on tori) | abstract |
| F. Latrémolière, "Continuity of the spectrum of Dirac operators of spectral triples for the spectral propinquity", Math. Ann. 389 (2024) 765–817, arXiv:2112.11000 | Propinquity convergence gives spectrum convergence, and action functionals are continuous for a wide class. | RELATED (spectral action in the limit, not gauge invariance at finite cut) | abstract |
| E.-M. Hekkelman, E. A. McDonald, "A noncommutative integral on spectrally truncated spectral triples…", J. Funct. Anal. (DOI 10.1016/j.jfa.2025.111154), arXiv:2412.00628 | Approximate NC integral on truncations; Szegő limit (Widom). | RELATED (Szegő/Widom machinery on truncations) | abstract |
| E.-M. Hekkelman, "Truncated geometry on the circle", Lett. Math. Phys. 112, 20 (2022), arXiv:2111.13865 | GH convergence of the pure state space of Toep_n. | RELATED | abstract (search snippet) |
| M. Leimbach, W. D. van Suijlekom, arXiv:2302.07877; Leimbach thesis (Radboud UP, 2025) | GH convergence of truncations of tori in every dimension; compact quantum groups with invariant Lip-norms. | RELATED (no gauge group) | abstract / snippets |
| A. Gaudillot-Estrada, W. D. van Suijlekom, arXiv:2310.14733; R. Toyota, arXiv:2309.13469 | Peter–Weyl truncations of compact groups. Dirac-type truncations are explicitly left "for future research". | RELATED | abstract |
| M. A. Rieffel, arXiv:2210.00387 (DOI 10.1016/j.geomphys.2023.104921); arXiv:2108.01136 | Fejér–Riesz systems for compact matrix quantum groups and ergodic actions; Dirac operators on matrix algebras converging to coadjoint orbits (S² has ≥ 3 inequivalent Dirac operators). | RELATED (equivariant setting, no G(E)) | abstract |
| Verhoeven, arXiv:2610.19791 (re-checked) | Full HTML read: no automorphism, isometry or gauge group, no order-one condition, no citation of Dong. §4.1 notes the SU(2) Dirac truncation is not a Peter–Weyl truncation. | **EXTENDS-OURS**: WG fills a gap in it | full text |
| J. W. Barrett, J. Burridge, "Fuzzy geometries with an internal space", arXiv:2604.19549 (2026) | Inner fluctuations on a matrix geometry produce a gauge-field analogue plus geometry fluctuations. | RELATED (matrix algebra, not operator system) | abstract |
| A. Connes, W. D. van Suijlekom, "Tolerance relations and operator systems", arXiv:2111.02903; Gielen–van Suijlekom arXiv:2207.07735 | C*-envelopes, propagation numbers and pure states for tolerance-relation operator systems. | RELATED (no symmetry group computed in the abstracts) | abstract |

## What we can build on

1. **Widom/Hankel form of the defect.** Write the order-one and pure-gauge edge sums as Hankel-type products P a Q · Q b P. This gives norm and trace-norm estimates from Toeplitz theory (Szegő/Widom asymptotics; the Hekkelman–McDonald route). It also gives a direct comparison with Glaser–Stern's O(1) operator norm and O(Λ) trace norm on S², which should be checked on S³ as the natural first test.
2. **Envelope argument for G(E).** Hamana's C*-envelope uniqueness, with Verhoeven's C*_env(S_J) = M_N, gives the following. Every unital complete order automorphism of the window is Ad U, so Dong's G(E) modulo phases is the full complete-order automorphism group. Farenick's isometry theorem is the circle precedent. WG's computed Lie algebra can be phrased as this automorphism group. The exact theorem number in Paulsen's book (ch. 15) is unverified, so cite Hamana generically.
3. **Orientation reversal.** The unitary that anticommutes with D is the S³ analogue of the Z₂ (reflection) factor in Dong's U(1) × (U(1) ⋊ Z₂). Isometries that reverse orientation sending D ↦ −D is standard spin geometry. Credit it as standard, not as new.

## Not found (searches empty)

- G(E), or complete-order automorphism groups, for any spectral truncation other than the circle: torus, S², SU(2)/S³.
- Any closed-form order-one defect at the spectral cut for any manifold. Glaser–Stern is numerical.
- Gauge invariance of the spectral action under a fixed spectral truncation, including Fejér/Toeplitz boundary terms (beyond the items in the earlier memo).
- Work by Kaad, Farsi, D'Andrea–Lizzi–Martinetti or Lizzi on truncation symmetry groups (only metric/state-space results surfaced).
- Any truncation result with spinors on S³ other than Verhoeven.
