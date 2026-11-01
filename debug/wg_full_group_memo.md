# WG follow-up: the full gauge group of the true CH window

Date 2026-10-31. Driver `debug/wg_full_group.py`. Data in `debug/data/wg_full_group.json` (full mode, n_max = 1, 2, 3) and `debug/data/wg_full_group_tonly.json` (torus part only, n_max = 4, 5). Builds on record WG (`docs/predictions/2026-10-31_window_gauge_group.md`). **This follow-up was not pre-registered.**

## Verdict

**FULL-GROUP, for unital complete order automorphisms, at n_max = 1, 2, 3.** Two results:

- Every unital complete order automorphism of O^CH is Ad U, with U in Dong's G. The theorem chain is cited below.
- Exactly as the record guessed, G = (U(1)·Spin(4)) ⋊ Z₂, so G/G₀ = Z₂. The Z₂ is the inversion. Consequence: the complete-order automorphism group G/U(1) is O(4), the isometry group of round S³. This is the S³ analogue of Dong's G(Toep_n)/U(1) = O(2).

## 1. Theorem chain: every complete order automorphism is Ad U

1. **C\*(O) = M_N.** Commutant dimension 1 measured at n_max = 1, 2, 3 (as Paper 38, Paper 32 prop:odd(c)); O unital and *-closed, so C\*(O) = O'' = M_N (double commutant, finite dimension). Verhoeven is not needed. His Thm 4.10 (arXiv:2610.19791, local full text) says the same, N = (2/3)(2J+1)(2J+2)(2J+3) = 4, 16, 40, 80, but his Def. 3.7 calls C\*(S) ⊂ M_n "the envelope", which is wrong in general; it is right here because M_N is simple (the Shilov boundary ideal is proper, hence 0).
2. **C\*-envelope universal property.** Connes–van Suijlekom, CMP 383 (2021), arXiv:2004.14115, Def. 2.24 (full text read): for every C\*-extension κ: E → B (complete order embedding, C\*(κ(E)) = B) there is a unique surjective \*-homomorphism ρ: B → C\*_env(E) with ρ∘κ = ι. Existence: Hamana, PRIMS 15 (1979) 773–785 (via CvS; also Paulsen Ch. 15, Blecher–Le Merdy §4.3); primary not read.
3. **Apply it.** Let φ be a unital complete order automorphism of O. Then κ = incl∘φ is a C\*-extension onto M_N, so some surjective ρ: M_N → M_N satisfies ρ∘φ = id_O. M_N is simple, so ρ is an automorphism. By Skolem–Noether, ρ = Ad W. Hence φ = Ad W\*|_O, with W\* ∈ G.
4. **Kernel.** Ad U is the identity on O exactly when U ∈ O' = U(1). So Aut_co(O) = G/U(1).

**Class covered.** Unital complete order automorphisms = unital surjective complete isometries (a unital complete contraction is CP; standard, Paulsen Ch. 3). CvS Prop. 2.8: a completely isometric CP isomorphism is automatically unital.

**Not covered:**
- (a) Non-unital complete isometries that are not completely positive. These would need the ternary-envelope Banach–Stone theorem, x ↦ U x W, which we did not check against a source.
- (b) Positive but not completely positive order automorphisms.
- (c) Antilinear maps, such as entrywise conjugation. O is closed under conjugation, to 1.9e-15. These maps are outside G by definition.

**Test case for (b).** Transposition (O closed under it, 1.6e-15) is not of type (b): U e = eᵀ U has a one-dimensional unitary solution (1e-14) at n_max = 1, 2, 3, commuting with D (4e-14) and fixing su(2)_L, so it lies in G₀. A random 2-positivity search found no violation (weak, consistent).

## 2. Component group: the reduction

G is a closed subgroup of U(N), so G/G₀ is finite. Let g ∈ G. Ad g preserves Lie(G₀) = u(1) ⊕ su(2)_L ⊕ su(2)_R (from WG1) and its derived ideal. The automorphism group of su(2) ⊕ su(2) is (SO(3) × SO(3)) ⋊ Z₂, so there is a homomorphism G → Z₂ that records whether the two factors are swapped.

**Kernel K.** Take g ∈ K. Ad g acts on each factor by an inner automorphism, Ad h with h ∈ G₀. The product g' = g h⁻¹ then commutes with the image of Spin(4). The Spin(4) representation on the window is multiplicity-free: the eigenspace ±(n + 3/2) carries ((n+1)/2, n/2) or (n/2, (n+1)/2), all distinct, which we measured from the Casimirs at n_max ≤ 5. By Schur, g' = e^{iθ(D)}. So K = G₀·T', where T' = {e^{iθ(D)} ∈ G}.

**T', computed exactly.** O contains exactly one copy of each harmonic type (k/2, k/2), k ≤ 2n_max − 1 (compression injective on band-limited functions). A function of D maps that copy into itself exactly when θ_λ − θ_λ' is the same for every block (λ, λ') in the copy's support. This gives integer constraints. Smith normal form of the constraint matrix:

| n_max | constraint rank | identity dim of T' | torsion | torsion mod Z(Spin4)·U(1) | residual of elements on O |
|:--:|:--:|:--:|:--:|:--:|:--|
| 1 | 3 | 1 | Z₂ | trivial | 2e-16 |
| 2 | 7 | 1 | Z₂ | trivial | 1.1e-15 |
| 3 | 11 | 1 | Z₂ | trivial | 1.9e-15 |
| 4 | 15 | 1 | Z₂ | trivial | 4.2e-15 |
| 5 | 19 | 1 | Z₂ | trivial | 6.7e-15 |

The single Z₂ is the central element (−1)_L of Spin(4). It already lies in G₀. Direct enumeration (orders 2, 3, 4), independent of the SNF, agrees.

Control: sign(D) is not in G. Its residual on O is 1.0, 0.99 and 1.0 at n_max = 2, 3, 4. At n_max = 1 sign(D) equals the central element, residual 2e-16, as it must. So **K = G₀**, and |G/G₀| ≤ 2.

**The swap is realised.** At n_max = 1, 2, 3 the inversion unitary ι is unique up to phase, unitary (1.4e-14), preserves O (2.6e-15), anticommutes with D (4.5e-14, so ι ∉ G₀), exchanges su(2)_L and su(2)_R (2.3e-15), and ι² is scalar (1.4e-14), so after a phase ι² = 1 and the extension splits.

**Result.** G = (U(1)·Spin(4)) ⋊ Z₂ at n_max = 1, 2, 3.

Spin(4) ∩ U(1) = {±(1, 1)}. The element (−1, −1) acts as −1 on every eigenspace, since a + b is half-odd. So G₀/U(1) = SO(4).

**Aut_co(O) ≅ O(4).** Pulling functions back along isometries of S³ gives a homomorphism O(4) → Aut(O). The map is well defined because compression is injective on band-limited functions (Paper 38 band ranks). Its kernel is trivial, because the degree-1 harmonics are the coordinates. It reaches both components. That gives Aut_co(O) = G/U(1) ≅ O(4). The circle analogue is Dong Cor. 4.3, where G(Toep_n)/U(1) = U(1) ⋊ Z₂ = O(2).

## 3. Status and scope

- **n_max = 1, 2, 3:** complete (commutant, T', ι here; dim g = 7 from WG).
- **n_max = 4:** dim g = 7 and T' trivial, so |G/G₀| ≤ 2; inversion not run. **n_max = 5:** T' only.
- **All n_max:** not proved; dim g = 7 and the T' pattern are measured, not derived.

**Prior art.** Isometries survive truncation (CvS §3); D ↦ −D under orientation reversal is standard; circle: Dong, Farenick (arXiv:2103.16546). New: the S³ computation and Aut_co ≅ O(4).

## 4. Before citing

- The component result is post-hoc: WG3 was "reported, not predicted", and this follow-up was not registered. For a paper statement at n_max ≤ 3, label it [MEASURED + reduction] and disclose that it was not pre-registered.
- Any general-n claim needs a pre-registration first. Suggested falsifier: G/G₀ = Z₂ at n_max = 4, meaning the inversion has null dimension 1 and anticommutes with D; T' torsion equals the centre at n_max = 6; and transposition lies in G₀.
- Guard is a separate step (§9); it should reject sign(D) in T', an extra T' element, ι² non-scalar.
- Cite Hamana via CvS Def. 2.24; no Paulsen theorem number.
