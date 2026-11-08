# Scoping memo: the deuteron's size from its own binding (candidate door 27)

Date: 2026-10-20. PI-directed scoping. Nothing was computed that scores a prediction, and no existing
file was modified. Companion draft: `debug/deuteron_sturmian_prereg_draft.md` (not committed).

**Verdict in one paragraph.** The door is buildable in a single session and is worth building, but it is
an **accounting door**, like J14 and J16, not a finding. Every candidate outcome is KNOWN or PROBABLY
KNOWN except one low-weight framing item. Use the **Malfliet–Tjon MT-III** triplet potential. It is
s-wave only, a sum of two Yukawas, and has published parameters and a published deuteron binding, which
we saw in primary text. It also has a structural advantage the other candidates lack: in an
exponential-tail (Laguerre/Coulomb-Sturmian) basis its matrix elements are **rational numbers**, so every
truncation is exact and π-free, while the same potential in the HO basis brings in `erfc`. The atlas
reading holds, in sharpened form. The radius stops being consumed and becomes an output, and the free
numbers move to the potential's parameters. At low energy they move further, to two numbers (the binding
B and the triplet scattering length a_t), because effective-range universality fixes r_m from those two
to about 1%. The charge radius also consumes the nucleon radii. One side finding needs a PI decision: the
ledger row that excludes Minnesota appears to rest on two code defects, not physics (Sec. 6).

## 1. What exists

- **Paper 23 `sec:deuteron`** uses Minnesota in an HO basis (N_shells = 2, ħω = 10 MeV) as an
  *encoding-validation* benchmark only. The "~2.1 fm" radius it quotes is experiment, not computed
  (l. 465–471). No computed deuteron radius exists in the corpus.
- **Paper 23 `sec:finite_nucleus`, seam J14/J16, atlas doors 24/26:** the nuclear radius R and the Zemach
  moment are **consumed** (verdict C). This door would supply R for A = 2 from the nuclear side.
- **`geovac/nuclear/`** has HO machinery (Moshinsky, `minnesota.py` with Gaussian HO matrix elements,
  `nuclear_hamiltonian.py`). There is **no** Laguerre/Sturmian code and no Yukawa-in-HO matrix element.
  Everything needed is new and small.
- **Paper 24 `sec:ks_frame` (J15):** the Coulomb Sturmians are the Kustaanheimo–Stiefel pushdown of the
  4D oscillator, so the exponential-tail basis carries the hydrogen skeleton's own (n, l) labels. It is used
  here only as a basis; by Fock rigidity (Paper 23 `sec:rigidity`) the Yukawa Hamiltonian is not
  diagonal in it.

## 2. Choice of potential

| Candidate | Channels | Matrix elements in Laguerre basis | Published deuteron numbers we could verify | Verdict |
|:--|:--|:--|:--|:--|
| **MT-III** (revised MT I–III triplet) | ³S₁ only, central | **rational** (Yukawa × exponential) | B = 2.23069 MeV, a_t = 5.5132 fm (Barletta–Kievsky, text seen) | **chosen** |
| MT-V | spin-averaged s-wave | rational | E_d = −0.3500 MeV (Elster, Schadow, Nogga, Glöckle, nucl-th/9805018, text seen) | rejected: not the deuteron |
| Reid soft core | ³S₁–³D₁ + tensor | Yukawa sums, but singular tensor pieces | not verified | later door |
| Argonne v18 | coupled, 18 operators | not rational: (1−e^{−cr²}) form factors | WSS 1995 table, from memory only | later door (D-state, Q) |
| Minnesota | central Gaussians | `erfc` (algebraic in HO) | — | excluded by PI (Sec. 6) |

**MT-III, as seen in primary text.** Barletta & Kievsky, arXiv:0807.2333, Sec. IV:
V(r) = (1438.72 e^{−3.11 r} − 626.885 e^{−1.55 r})/r, with r in fm and V in MeV, ħ²/m = 41.47 MeV fm²,
acting on l = 0 only, "E2b = −2.23069 MeV … as = 5.5132 fm". Garrido, Romero-Redondo, Kievsky and
Viviani, arXiv:1211.1802, Sec. III, give the same triplet form (with singlet 513.968), call it the
"revised Malfliet–Tjon I-III" of Friar et al., PRC 42, 1838 (1990), and give "binding energy for the
deuteron of 2.2307 MeV". Original: Malfliet & Tjon, Nucl. Phys. A127, 161 (1969) (not consulted).
**Published MT-III deuteron rms radius: NOT FOUND** in three web searches. It is probably tabulated in
the 1980s Friar–Payne–Gibson benchmark papers, which we did not consult. A_S and P_D: none published
that we saw. P_D = 0 by construction.

**Caveat on the forced/free reading.** B and a_t were consumed **upstream** by Malfliet and Tjon's fit to
NN data. MT-III's B is within 0.3% of experiment; whether B was a fit target is not verified. We introduce
**zero** fitted parameters. The potential carries five decimal numbers (two strengths, two ranges,
ħ²/m), and they encode NN low-energy data.

## 3. Literature (what was seen, and how)

- **Caprio, Maris, Vary, PRC 86, 034312 (2012)**, arXiv:1208.4156. Abstract seen: Gaussian HO falloff is
  a poor description of the tail, and Coulomb-Sturmian functions give "a complete, discrete set" with
  exponential falloff. Example ⁶Li (details via search summary only).
- **Caprio, Maris, Vary, PRC 90, 034305 (2014)**, arXiv:1409.0877, ar5iv text read through a fetch
  summary. Their radial function is S_nl(b; r) ∝ (2r/b)^{l+1} L_n^{2l+2}(2r/b) e^{−r/b}, orthonormal in
  dr. **For ⁴He the CS basis
  does not converge the radius better than HO**: r_p is 1.4358 (CS) vs 1.4361 fm (HO), and CS energy
  convergence lags HO by about two N_max steps. For ⁶He/⁸He they report "improved ℏΩ independence" of
  radii, but r_m is "not possible to give a definitive value". **No deuteron.** Our gain therefore rests on
  the deuteron being a halo-like, weakly bound system, which is the case where CS should help most. That
  is an expectation, not something they show.
- **More, Ekström, Furnstahl, Hagen, Papenbrock, PRC 87, 044326 (2013)**, arXiv:1302.3815, text read.
  An HO basis acts as a hard wall at L₂ = √(2(N + 3/2 + 2)) b. Their Fig. 5 shows the deuteron r² (for
  the Entem–Machleidt N3LO potential) falling on one curve versus L₂, converging to 3.9006 fm². **The slow,
  IR-limited HO radius convergence for the deuteron is KNOWN.**
- **Buchmann, Henning, Sauer, nucl-th/9610018**, text read. Their Eq. (3) is
  r_d² = (r_ch² − r_p²) − r_n² − r_DF², with r_DF² = 3/(4M_N²) = 0.0331 fm². Their Eq. (5) adds the
  spin-orbit term and the two-body (meson-exchange) term r²_[2] ≈ 0.02 fm². They state that a_t, A_S and
  r_m are empirically linearly related, with r_m = 1.969(4) fm at a_t = 5.419(7) fm and A_S = 0.885(2)
  fm^{−1/2}, and that scattering analyses give r_m ≈ 1.950 fm (a 1% tension).
- **Effective-range theory.** Bethe, Phys. Rev. 76, 38 (1949), not consulted. Phillips, Rupak, Savage,
  nucl-th/9908054, seen via a search summary: at NLO the Z-parametrization "exactly reproduces …
  effective range theory … of a deuteron containing point-like nucleons". Klarsfeld et al., NPA 456, 373
  (1986), cited by Buchmann, not consulted.
- **Measured values.** B = 2.224566 MeV (PI-supplied; primary source not checked). r_ch = 2.12799(74) fm,
  CODATA 2018 (seen via a search summary of arXiv:2307.09680). r_p = 0.8414(19) fm (CODATA 2018) and
  r_n² ≈ −0.116 fm² (PDG) are from memory and must be re-read from CODATA 2022/PDG at scoring. The gates
  below are insensitive to them at the 0.1% level.

## 4. Design

**Radius is a ground-state expectation value.** ⟨r²⟩ = ⟨ψ₀|r²|ψ₀⟩ needs one eigenvector, with no sum
over excited or continuum states. The failure mode of the ledger row ("sum-over-states fragile for
continuum response", polarizability unstable from 0.4 to 19,000 fm³ across α) therefore **does not
apply**. The variational principle bounds B. ⟨r²⟩ has no bound, but its convergence follows the
eigenvector's.

**Basis.** φ_n(r) = N_n x L_n^{(2)}(x) e^{−x/2} with x = 2λr (l = 0; Caprio's L^{2l+2} form,
orthonormal in dr, S = I). The non-orthogonal Coulomb-Sturmian set with L^{(1)} spans the same space at
each N, so the choice between them affects only conditioning. Every integral reduces to
∫x^k e^{−a x} dx = k!/a^{k+1}, which gives:
- overlap I; kinetic (ħ²/m)λ² × (rational matrix);
- ⟨e^{−μr}/r⟩ = λ × (rational function of μ/λ); ⟨r²⟩ = λ^{−2} × (rational matrix).
With λ rational, H_N ∈ ℚ^{N×N} in MeV. Plan: exact `Fraction` assembly, then an `mpmath` eigensolve at
dps 50.

**Scale λ, with no fitting.** Primary rule: λ_N minimizes E_N at each N (the variational rule; it uses no
physics input). We also report a scan λ ∈ {0.5, 1, 2, 3, 4, 6, 8} × κ_N and the λ-plateau width. A
self-consistent tail rule (λ = κ_N) is run as well, but **excluded from the tail test**, because there it
is tautological (Sec. 5).

**Comparison basis.** The HO relative s-wave on the same MT-III potential, with the same function count
N_r (N_quanta = 2(N_r − 1)), at ħω ∈ {2, 4, 6, 8, 10, 15, 20, 30, 40, 60} MeV. The Yukawa HO matrix
elements are evaluated in closed form, using e^{y²}erfc(y) with y = μb/2, at dps 30. They are
cross-checked by quadrature.

**Independent reference route.** Direct radial ODE: series start u = r + c₂r² at the 1/r core, Numerov
on a fine grid with Richardson extrapolation in h, r_max = 60 fm, shooting on B. This route gives B,
r_m, A_S (from u e^{κr} at large r) and a_t (from the zero-energy solution).

**Computed quantities.** B; r_m = ½√⟨r²⟩ (point nucleons); κ_N = √(B_N/41.47) fm^{−1}; the plateau
A_S(r) = u(r) e^{κ_N r} on r ∈ [8, 16] fm; r_ch from Buchmann's Eq. (5) with only r_p², r_n² and r_DF²
(MT-III has no model for the spin-orbit or meson-exchange terms; their size, about 0.02 fm², is
disclosed); and r_ERE from the computed κ and a_t (shape-independent form,
r_m² ≈ 1/(8κ²(1 − κρ)) with ρ = 2(κ − 1/a_t)/κ²).

**Hand estimates (disclosed, not computed).** For MT-III: κ = 0.23193 fm^{−1}, ρ = 1.880 fm,
A_S ≈ 0.907 fm^{−1/2}, r_ERE ≈ 2.030 fm. The same formula on experiment (κ = 0.23161, a_t = 5.419)
gives 1.982 fm against Buchmann's 1.969, an overshoot of +0.6%. Applying that overshoot gives
**r_m(MT-III) ≈ 2.017 fm** and r_ch ≈ 2.166 fm, about +1.8% above CODATA.

## 5. Predictions (summary; full text in the draft)

| ID | Statement | Gate | Novelty |
|:--|:--|:--|:--|
| DR1 | At N_r = 30 the Laguerre r_m is converged to ≤1e-6 relative (vs ODE). HO at the same function count is ≥1e-4 off at every ħω in the scan | both halves | HO half KNOWN (More 2013); deuteron Laguerre-vs-HO NOT FOUND (low weight) |
| DR2 | At λ = 2κ_N and 3κ_N (not κ, to avoid the tautology), N_r = 40, u e^{κ_N r} is flat to 1e-4 on [8, 16] fm, and its value equals the ODE A_S to 1e-4. The HO local decay rate at N_r = 20 deviates from κ by >10% somewhere in the window | flatness, A_S match, HO control | PROBABLY KNOWN |
| DR3 | Number class. Laguerre-basis H_N and ⟨r²⟩ are rational for rational λ. The HO Yukawa element has a nonzero e^{y²}erfc(y) part. Minnesota-in-HO is algebraic and Minnesota-in-Laguerre has erfc | exact symbolic | KNOWN (math); the pairing as an atlas statement NOT FOUND |
| DR4 | r_m(MT-III) ∈ [1.99, 2.05] fm. r_ch exceeds CODATA by +0.8% to +3.5%. Basis error ≤1e-5 relative, at least 100× below the potential's residual | value window | PROBABLY KNOWN (value likely in the Friar–Payne–Gibson benchmarks) |
| DR5 | Forced/free split: 0 < (r_ERE − r_m)/r_m < 2%, with r_ERE from the **computed** κ and a_t. So the radius's free content is two low-energy numbers | window | KNOWN in form (ERE; Klarsfeld 1986; Phillips–Rupak–Savage) |

Controls: B = 2.23069 MeV and a_t = 5.5132 fm reproduced to 1e-5 (positive, KNOWN, Barletta–Kievsky); a
hydrogen limit (V = −1/r, λ = 1, atomic units) exact at N = 1 (positive, exact); the MT-I singlet
(513.968) has no bound state, and its lowest eigenvalue stays positive as N grows (negative, KNOWN); the
ODE route versus the basis route (independent-route rule).

**Expected outcome.** At best four CHECKs and one low-weight FINDING. The likeliest MISS is DR1's HO
half (HO reaching 1e-4 at a good ħω). Thresholds come from hand estimates; none was tuned.

## 6. The Minnesota exclusion: the ledger's lesson looks wrong (PI decision)

Reading the v2 code (`debug/deuteron_sturmian_v2.py`) against `geovac/nuclear/minnesota.py` shows two
defects:
1. **Wrong channel.** Lines 33–34 label `V_T, K_T = -91.4, 0.465` as "triplet". Those are Minnesota's
   **singlet** values (`minnesota.py` l. 56–59: V_S = −91.4, κ_S = 0.465; V_T = −178.0, κ_T = 0.639).
   The refit (v3) inherits κ = 0.465.
2. **Non-orthonormal basis treated as orthonormal.** The basis is x L_n^{(1)}(x) e^{−x/2} in dr, with
   S = I assumed. By hand: ∫x² L_0^{(1)} L_1^{(1)} e^{−x} dx = ∫x²(2 − x)e^{−x} dx = 4 − 6 = −2 ≠ 0, and
   ∫u_0² dr = 2/α ≠ 1. The kinetic matrix used is for a different basis, so the variational bound is
   lost and arbitrary over-binding follows.

The row's stated cause ("Minnesota is effective interaction for HO") is therefore **unsupported**.
Minnesota was fit to free-space NN data (Thompson–LeMere–Tang 1977, not consulted). Its deuteron binding
in a correct basis is not established in the corpus. We did not compute it, because doing so would
score an outcome. **Recommendation (PI's call):** append a RE-INTERPRETATION row to §3 and the ledger,
in the 2026-09-25 CUSP-2 pattern (append only, §13.5). Keep MT-III as the primary for this door, and
optionally add "Minnesota in a correct Laguerre basis" as a separate, separately registered control.

## 7. Atlas reading (proposed door 27, for the PI)

- **Releases:** the consumed radius of J14/J16, which becomes an output for A = 2.
- **Survives:** l = 0 and the exponential-tail labels as a basis. The Hamiltonian is not diagonal in
  them (Fock rigidity).
- **Price:** the potential's five numbers, which collapse at about 1% to (B, a_t); plus r_p and r_n²
  (consumed) for r_ch.
- **Number class:** rational at every truncation in the matched basis, `erfc` in the mismatched one.
  Whether the N → ∞ limit is transcendental is unknown and not claimed.
- **Verdict:** FREE (moved, not removed), with the convergence structure FORCED if DR1–DR3 hold.
- **Scope:** it does not reach A ≥ 3 (Faddeev/NCCI), the tensor force, the D-state, or Q. AV18 would be
  the next door.

## 8. Cost and risks

**Cost.** About 400 lines in `debug/` (Laguerre exact assembly, HO-Yukawa closed form, ODE route,
driver), with minutes of runtime at dps 50. Guards come afterwards, as separate work (§9 guard rule).
The paper home, if it lands, is a new Paper 23 subsection after `sec:deuteron`, with seam row J17.

**Risks.**
1. **Tautology in the tail test.** At λ = κ the finite expansion decays as r^N e^{−κr} by
   construction, so DR2 excludes λ = κ.
2. **Laguerre conditioning at N ≥ 40:** handled by exact rationals plus mpmath.
3. **The 1/r core.** The net 811.8/r MeV repulsion at the origin gives u ∝ r + c₂r². The HO basis
   carries only odd powers of r and so converges algebraically in the UV; this is part of why DR1
   expects HO to lag. If HO converges fast anyway, DR1 misses.
4. **Overclaim risk.** The radius depends on NN data through the potential, so it is not "derived" in any
   absolute sense. The papers stay under §1.5.
5. **The published MT-III radius is unverified.** DR4's window rests on hand ERE plus a 0.6% transfer;
   if the window misses, the transfer was wrong, not the method.
6. **Novelty.** Accounting-grade; Caprio et al. saw no CS radius gain for ⁴He.
