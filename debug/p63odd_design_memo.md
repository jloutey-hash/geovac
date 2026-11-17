# Paper 63, owed items after P63L1: odd-λ pieces and l > 0 attainment. Design pass (2026-11-17)

Status: design only; nothing scored; no paper, test or CLAUDE.md edits; no new shape evaluated at eps < 0.1. Drivers: `debug/p63odd_{factor,tail,engine,controls,toy,mpspot}.py`; data:
`debug/data/p63odd_*`.

Owner text: `prop:tree_rung_l` ends "pieces with odd λ are bounded but not identified with the operators of
F_k, and attainment for ℓ > 0 is shown only for the p_z-leaf chain". Notation: P_λ := φ(L_m(θ − λ)).

---

## (A) The odd-λ pieces

### A.1 Where odd λ comes from (parity lemma). Proof sketch, bookkeeping only

Track λ mod 2 through the three steps of `prop:tree_rung_l`.
- **Densities.** A term ∂_a^i ∇^{2j} h_ℓ(∂) has radial part r^{-ℓ}·poly(θ)·K0 (Hobson). ∇² = r^{-2}θ(θ+1) and
  a∂_a = θ are even. So λ ≡ ℓ (mod 2).
- **Vertices.** Channel L of a product of harmonics has L ≡ Σℓ (Gaunt parity), and radial twists add. So in every
  channel, λ ≡ L. Λ²B is trivial (constant Wronskian), so the Clebsch–Gordan step also keeps parity.
- **Potentials.** Φ_L W lies in r²N_W (parity kept) plus two order-one pieces, r^{L} (parity kept) and
  r^{-L-1} (parity **flipped**). Each order-one piece replaces the whole subtree's Bessel content with a pure power,
  so it drops at least one Bessel factor.
- **Root.** It pairs channel L (λ ≡ L) with the root's r^{-ℓ0} part, ℓ0 ≡ L, and H's r² factor. The total is even.

**Lemma (parity).** A piece with odd λ arises only through an odd number of exterior tails r^{-L-1}. In
particular it has m ≤ n − 2 Bessel factors. Every piece with m = n − 1 has even λ.

Examples: in the s-chain of three, χ = π/2 − R gives −πR ∈ rB (λ odd, m = 1); in P63L1, g = −1/r² + E, the
r^{-2} is the L = 1 exterior tail and the cross term −2E has m = 1. Not separately guarded.

### A.2 What the odd pieces are. Proposition, proved for m ≤ 5

**Proposition (odd pieces = 2-D banana).** For odd λ, the composition factors of P_λ = φ(L_m(θ−λ)) of order > 1
are those of V_m(s = t²). V_m is Vanhove's operator for the equal-mass two-dimensional banana with m propagators
(arXiv:1401.6438, Table 1), of order m − 1. All other factors have order one. For even λ, the factors of order > 1
are those of φ(L_m), the operator of F_{m+1}, **for every m**. The paper states this only for m + 1 ≤ 8, via
irreducibility.

*Proof.*

(i) **λ = 1.** The transfer rules of `prop:tree_rung` hold with K0(tr) replaced by I0(tr), which obeys the same
θ-equation. The boundary terms vanish for |t| < m. So P_1 annihilates b_m(t) = ∫ x I0(tx) K0(x)^m dx. Up to the
factor 2^{m−1}, this is Vanhove's banana I_m(s) at s = t² (his eq. 9.1).

(ii) **Exact division** (`p63odd_factor.py`, E1/E2). P_1 = c_m t^{2a} (θ − r_m)² V_m exactly, for m = 2, 3, 4, 5,
with (c_m, r_m) = (−2, 0), (4, 2), (−8, 2), (16, 4).

The quotient has constant coefficients, so its two factors are order one. Control E4: V_m does *not* right-divide
the even operator P_0 (m = 2..5). The right-factor fact was already in the archived banana-dictionary sprint
(`debug/archive/misc/banana_dict_vanhove.py`); the quotient is new.

(iii) **All odd λ: parity transport.** The intertwiner P_{λ+2}U = U P_λ holds, with U = t^{-2}θ². E3 checked it
exactly for m = 2..5 and λ = −3..5 (odd); the guard checks four more cases. It gives a module map
D/DP_{λ+2} → D/DP_λ, [1] ↦ [U]. Its cokernel is a quotient of D/Dθ², so it has order-one factors only. Hence the
factors of P_λ of order > 1 are a sub-multiset of those of P_{λ+2}.

The adjoint identity U*P*_{λ+2} = P*_λ U* gives the map D/DP*_λ → D/DP*_{λ+2}. Its cokernel also has order-one
factors only. Dualizing (D/DP)^∨ ≅ D/DP* over C(t) gives the reverse inclusion. So the multiset of factors of order
> 1 depends only on λ mod 2. For odd λ it is that of P_1, which is V_m by (ii). For even λ it is that of P_0. ∎

**Numerical checks** (`p63odd_factor.json`, N1/N2, mp30):
- P_λ annihilates b_m^{(λ)} = ∫ x^λ I0(tx) K0^m for λ = 1, 3 to ~1e-31 relative at two t.
- The even operator applied to the odd function leaves 0.65 relative (the control).
- V_2 b_2 = −1, constant to 1e-31 at t = 0.5 and 0.9: the inhomogeneous banana equation.
- m = 3, 4: P_λ annihilation ≤ 2e-31; controls 0.14 and 0.022; V_m b_m = −1.5 and −3, i.e. −m!/2^{m−1} as
  Vanhove's normalization gives, constant to ≤ 1e-31.

### A.3 What this does to the bound

Combining A.1 and A.2, every composition factor of an l > 0 tree fibre (degree ≤ 2) is one of:
- of order one;
- a factor of F_{m+1}'s operator, m ≤ n − 1 (even λ); or
- a factor of the banana operator V_m, m ≤ n − 2, of order m − 1 ≤ n − 3 (odd λ).

So odd pieces sit at least three orders below the bound n (A.2 alone gives n − 2).

For n ≤ 4, V_m with m ≤ 2 has order ≤ 1. This explains the paper's observation that "none [of the factors computed
at n = 3, 4] comes from an odd moment." The first possible non-trivial odd piece is the sunrise V_3 (elliptic,
singular at t² = 1, 9) in a fibre of n ≥ 5. V_3 is not equivalent to F_2's Legendre operator: it has a genuine
singularity at t = 3 (exponents 0, 0, so a log), where F_2 is regular. Likewise V_4 has a singularity at t² = 16,
which F_3 lacks.

### A.4 What remains open (precise)

1. **m ≥ 6** (needed for n = 8 trees). That V_m has order m − 1 and is the right factor of P_1 is literature
   (Vanhove's Table 1 and the Borwein–Salvy recursion). It is not checked here for m = 6. The parity transport in
   (iii) holds for every m.
2. **Presence.** Whether any fibre actually carries a V_m piece, m ≥ 3, is not shown. Attainment of odd pieces is a
   separate question; the lemma says it needs n ≥ 5.
3. **Irreducibility of V_m** (m = 3..5) is cited (sunrise elliptic; Sym² for m = 4; Calabi–Yau 3 for m = 5), not
   checked here. The bound statement does not need it.
4. The parity lemma is a sketch. It needs to be written into the proof of `prop:tree_rung_l` and fire-tested by a
   separately written guard. The candidate wrong answer is "an odd piece with m = n − 1".

---

## (B) Attainment for l > 0 beyond the p_z-leaf chain

### B.0 The mechanism (local rule). It is stronger than a shape test

**Local rule.** For a source s = e^{-kr} r^a Σ s_j r^{-j} in channel L, the e^{-kr} part of its Coulomb potential
is the unique formal solution of (θ−L)(θ+L+1)Φ = −4πr²s. Its leading term is −4πs_0/k², independent of L.
L first enters at relative order r^{-2}.

Analytic controls (`p63odd_tail.py`, all pass, exact): T1 P63L1's leaf tail √(2/π)r^{-1/2}e^{-r}[1 + 11/8r +
97/128r² − 463/1024r³] (its "r^{1/2}e^{-r} cancellation" is just −4πf_p); T1b the s-leaf tail −(2/π)R/r; T2
P63L1's A = 1/(3π²), κ = 47/8; T3 SA's A = 1/π², κ = 31/8; T4 star of five, c_5Γ(1/2) = π²/(8√2).

**Consequence (proof sketch; a candidate proposition, PI call).** Take a tree on the slice with t at the root v0
and n ≤ 8. The e^{-(n−1)x} part of H has leading coefficient equal to a product of non-zero radial constants
times ∫_{S²} ∏_v σ_v(−ω) dω. Here σ_v is the principal symbol of P_v: ∂ → −ω, ∇² → 1, the top ∂_a power. The
coefficient is the "all electrons at one far point" limit, because every potential's exponential part is local.

Then the star proof goes through verbatim:
- The other exponential orders of H lie in modules with ≤ n − 2 Bessel factors. Each tail-vs-moment choice at an
  edge drops the exponential and the Bessel factor together. So they are holomorphic at −(n−1).
- The I-part of the root kernel then gives the non-holomorphic term x^{(4−n)/2+Σi_v}e^{-εx}, with
  i_v the top ∂_a power of P_v. With no ∂_a, these are the same exponents as `prop:star_attain`.
- By the parity lemma, the surviving c = n − 1 piece has even λ, so it is ≅ F_n.

So:
- **(a)** every **s-type tree**, n ≤ 8, on any one-scale slice (root at the t density) reaches rung n. Here σ ≡ 1
  and the angular integral is 4π. This would close the paper's "for other trees with n ≥ 5 it is not shown".
- **(b)** an l > 0 tree with ∫∏σ_v ≠ 0 reaches rung n.

For a star of pure harmonics, ∫∏σ_v = 0 forces the fibre to vanish identically. For trees with internal vertices
it does not (B2). The step needing care in a write-up: the exponential-order split of H on [1, ∞) into module
pieces, with boundary terms at x = 1 as in the star proof.


### B.1 Shape B1: p_z leaf – p_z middle – d root, t on the root (exercises all three P63L1 caps)

The guard's reduction chain (`tests/test_paper63_tree_rung_l.py`): an internal p_z density, an L = 2 channel
(cos² = 1/3 + (2/3)P2), a d root h_2(∂)t³K0(tr)/(2π²) = t⁵K2(tr)P2/(2π²).

The fibre is (4π/5) t⁵/(2π²) ∫ K2(tr) r² Φ_2[(2/3)f_p g] dr.

**Derivation.**
- At t = −(2−ε)+i0, K2(ye^{iπ}) = K2(y) − iπI2(y). Only the e^{-2r} part of Φ_2 reaches −2.
- The source is w = (2/3)f_pE ~ −(1/(3π²))x^{-1}e^{-2x}(1 + 7/4x).
- The local rule with k = 2, a = −1 gives Φ ≈ −πw through first order (φ_1 = −πw_1 + (a+1)φ_0, and a + 1 = 0).
- With I2's −15/(8y) and (2−ε)^{-1/2}, the branch integral S1(ε) = ∫_1^∞ I2(y) x²Φ_{2,exp} dx satisfies

  **ε^{3/2}S1 → A1 = 1/(12π) = 0.0265258,  κ1 = 15/8** (sympy, `p63odd_tail.json`; the hand derivation agrees).

**Predictions.**
- FB1: last-5 log-slope in [−1.55, −1.45] (p = 0.93).
- FB2: |A_fit/A1 − 1| ≤ 1e-5, smallest 6 points (p = 0.92).
- FB3: u² coefficient/(A1κ1) ∈ [0.95, 1.05] (p = 0.75; non-gating).

### B.2 Shape B2: the same chain with a cancelling root (genuine risk)

**Object.** The root is (∇² − (5/2)h_2(∂))t³K0(tr)/(2π²). Its radial parts are t⁵[(K0 − K1/ρ) − (5/2)K2 P2].
∫cos²θ(1 − (5/2)P2) dΩ = 0 (checked), so the leading e^{-2x} amplitude **cancels exactly**. The fibre itself is
non-zero: −3.07e-3 at t = 0.5, root swap 4e-15 (control C3). This mixed root has degree two, so it is inside
`prop:tree_rung_l`. The proposition of B.0 does not apply, because ∫∏σ = 0.

**Derivation.** Φ_0 and Φ_2 agree through relative order 1/x (a = −1), so the surviving term comes only from
the root kernels: I0 + I1/y ~ 1 + 9/(8y) against I2 ~ 1 − 15/(8y), difference 3/y.

So the exponent shifts by one (−3/2 → −1/2), and S2 = ∫_1^∞ x²[(I0 + I1/y)Φ_{0,exp} − I2Φ_{2,exp}] dx satisfies

**ε^{1/2}S2 → A2 = 3/(8π) = 0.1193662,  κ2 = −3/8.** The ε^{−3/2} coefficient is exactly 0. The expansion carries
odd powers of u = √ε from the analytic part, so the fit basis is {1, u, u², u³, u⁴}.

Attainment for B2 rests on a sub-leading coefficient being non-zero: the P63L1 situation one level up. The
coefficient is derived non-zero, so B2 also attains by proof, given A.1. The record checks the exponent shift and
the root-kernel amplitude.

**Predictions.**
- FB4: last-5 slope in [−0.75, −0.25] (p = 0.9; the toy gives −0.513, slowed by the u term).
- FB5: |A_fit/A2 − 1| ≤ 1e-4 (p = 0.85).
- FB6: u² coefficient/(A2κ2) ∈ [0.9, 1.1] (p = 0.6; non-gating).

### B.3 Controls (measured; `p63odd_controls.json`, `p63odd_toy.json`, `p63odd_mpspot.json`)

| id | what | result |
|:--|:--|:--|
| P | positive control: frozen P63L1 instrument re-run (sha256 match) | 9 S_mp values identical to the committed run (max rel diff 0.0); A ratio 0.9999999883, slope −1.49775 |
| C1 | generic tail route (L = 1, k = 1) vs P63L1 closed form e^xE, x = 1..80 | 7.8e-15 |
| C2 | P63L1 branch rebuilt from the generic route vs committed mp30, 9 eps | 1.1e-13 (42 s) |
| C3 | reduction root vs far leaf, t = 0.5, 1.3 | B1 2e-16 (value = guard's MC-pinned 3.31319e-4); B2 4e-15, 4e-15 |
| C4 | integrand e^{2x}Φ_Lx: Laguerre vs cumulative-grid route, L = 0, 2 | 4.4e-10, 5.4e-10 (grid-limited) |
| M | integrand: float Laguerre vs mp30 quadrature, L = 0, 2, x = 1, 3, 10 | ≤ 7.1e-15; 29–49 s per mp node |
| C5 | integrand vs the derived 8-term local-rule series, x = 20..60 | 3e-8 → 6e-12 (falls as x^{-8}) |
| T1 | TOY1 (I2 kernel, H = xe^{-2x}; A = 1/4, κ = −13/8), float vs mp30, fit | ≤ 4e-14; A ratio −2.9e-9 (smallest 6); slope −1.5006; u²/(Aκ) 0.99995 |
| T2 | TOY2 (B2's cancelling kernel, A = 3/4, κ = 9/8) | ≤ 1.8e-12; A ratio −3.2e-9; slope −0.5126; u²/(Aκ) 0.99992 |
| C6 | B1, B2 branch at eps = 0.5, 0.2, 0.1 only (values not read): 10k vs 20k trapezoid nodes; wall time | 1.5e-9 to 5.6e-9 (B1), 3e-10 to 2.5e-9 (B2), i.e. O(h²); 0.4 ms (B1), 0.9 ms (B2) per node |

### B.4 Gate, caps, cost (for record "P63L2", two shapes, scored independently)

**Gate.**
- **VOID:** FB0 fails, i.e. the float route and the referee disagree by more than 1e-9 relative at any eps, or a
  frozen file changed, or a step has not finished 48 h after launch. The referee is mp30 for x < 32 plus the
  derived series (≥ 40 terms) for x ≥ 32. In that region the float integrand is compared with the derivation's input, as with P63L1's R series.
- **GO:** for B1, FB1 and FB2 hold. For B2, FB4 and FB5 hold.
- **STOP** (per shape):
  - B1: slope outside [−2, −1].
  - B2: slope > −0.25. The branch is absent and the cancellation reaches the next order: attainment fails at this
    order.
- **BORDERLINE:** otherwise. This includes B2 slope < −1, which would mean the angular cancellation failed: a design
  error, reported as such.

**Caps.**
- One slice; t on the root; degree ≤ 2.
- B2's cancellation is engineered.
- Both attainments are provable from the derivations plus A.1. The record tests the derivation and the
  full-integrand instrument, not an open conjecture (the same honest weight as P63L1).

**Cost** (measured):
- Float branch at the P63L1 grid (400k nodes): about 160 s (B1) and 380 s (B2) per eps, so 9 eps is about
  25 min + 60 min.
- mp referee: about 200 cached nodes in [1, 32] (the panel edges do not depend on eps for eps ≤ 0.06) × 30–50 s ×
  number of channels. That is about 2.5 CPU-h for B1 and 5 CPU-h for B2; 4-way parallel, about 2 h.

**Recommended order.**
1. Write A.1 + A.2 into the paper as a proposition. It is proved for m ≤ 5, and the guard is a separate pass (E1–E4
   are its natural core).
2. Register P63L2 (B1 + B2).
3. PI call: write B.0 as a proposition (it changes an OPEN item, "other trees with n ≥ 5"), proof first, guard
   in a separate pass. B1 is then a check of its l > 0 case; B2 lies outside it.
