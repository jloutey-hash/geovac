# Sprint memo — Does the TC three-body operator L3 collapse to two-body under GeoVac's angular selection rules?

**Date:** 2026-08-23 · **Type:** Decomposer structural probe (angular-only, single-center)
**Verdict: STOP** — L3 is genuinely three-body under Gaunt/6j selection rules. No collapse. xTC (1-RDM contraction) is the only reduction route; it inherits GeoVac's angular sparsity.
**Driver:** `debug/tc_threebody_collapse_angular.py` · **Data:** `debug/data/tc_threebody_collapse.json`

## 1. The plane-wave mechanism (why PW collapses)
Transcorrelation similarity-transforms `H → e^{-τ}He^{τ}`, `τ=Σ_{i<j}u(r_ij)`. The `-½(∇_i u)²`
term splits into a two-body piece (`j=k`) and the genuine three-body operator
`L3 = -½ Σ_i Σ_{j≠k≠i} ∇_i u(r_ij)·∇_i u(r_ik)`. Verbatim abstracts: Luo & Alavi
(arXiv:1712.07524, HEG/plane waves) — the effective Hamiltonian "takes a simple form for
plane wave bases, **containing up to two-body operators only**"; Cohen–Luo–Guther–
Dobrautz–Tew–Alavi (arXiv:1908.02882, Gaussian basis, first-row atoms) — it "**contains
three-body interactions**." The PW mechanism: with `u(r)=Σ_q ũ_q e^{iq·r}`, the two
correlator lines meeting at the shared vertex `i` deposit momenta `q` (i–j line) and `q'`
(i–k line); the plane-wave product `e^{iq·r_i}e^{iq'·r_i}=e^{i(q+q')·r_i}` **adds the two
transfers into a single resultant `q+q'`**. Translation symmetry is abelian → the resultant
is one vector (multiplicity 1) → the three-body integral is the analytic number
`ũ_q ũ_{q'}(q·q')` fully pinned by momentum conservation, and in the homogeneous PW basis
it reduces to two-body-and-lower effective operators.

## 2. L3's GeoVac angular structure
Same Legendre/Gaunt machinery as Paper 22: `u(r_ij)=Σ_{L,M}(4π/(2L+1))u_L(r_i,r_j)
Y*_{LM}(r̂_j)Y_{LM}(r̂_i)`. Particles j, k each see ONE correlator harmonic (ordinary 2-body
Gaunt). The shared vertex i carries **two** correlator harmonics plus its bra/ket — the
four-harmonic integral `W(a;L,L';c)=∫ Y*_{l_a m_a}Y_{LM}Y_{L'M'}Y_{l_c m_c}dΩ_i`, which
decomposes exactly through the SO(3) tensor product `Y_{LM}Y_{L'M'}=Σ_Λ G_pair(L,L',Λ)
Y_{Λ,M+M'}` as `W=Σ_Λ G_pair(L,L',Λ)·G_ext(a,Λ,c)`. **Λ is the SO(3) analog of `q+q'`.**
Crucially SO(3) is *non-abelian*: Λ runs over `|L−L'|…L+L'` intersected with the external
triangle `|l_a−l_c|…l_a+l_c` — a Clebsch–Gordan **sum over multiple channels**, not one
resultant. (Scope: gradient's radial·radial part gives exactly this four-Y; the angular·
angular part only redistributes within the same L,L' via 6j — measuring W's Λ-rank is the
faithful, conservative collapse test. Verified against spherical quadrature, max err 1.9e-15.)

## 3. Measured rank / sparsity (l_max ≤ 2, single center)
| Object | l=1 | l=2 (L_corr=2) | l=2 (L_corr=3) |
|---|---|---|---|
| Shared-vertex bilinear-form max rank (1 = collapse) | **4** | **9** | **16** |
| Frac. external pairs with rank ≥ 2 | 0.875 | 0.975 | **1.000** |
| Frac. four-Y with SO(3) multiplicity n_Λ ≥ 2 | 0.19 | 0.49 | 0.49 |
| Full L3 angular-tensor density | 6.1% | **5.0%** | — |
| Full L3 Schmidt rank (i : jk bipartition) | 9 | 25 | — |
| **PW-mimic residual** (norm in higher SO(3) channels) | **14.5%** | **28.4%** | — |

The shared vertex is **rank ≥ 2 for 87.5%–100%** of external pairs; max rank grows 4→9→16
with the basis. Keeping only the lowest resultant Λ (the abelian PW mimic) discards **14.5%
of the operator norm at l=1, rising to 28.4% at l=2** — the higher SO(3) channels are
load-bearing and grow with l_max. The angular selection rules DO make L3 sparse (≈5% dense,
same single-digit-% order as Paper 22's 8.5% two-body density), but sparsity ≠ collapse.

## 4. Verdict — STOP, with the honest reason
Momentum conservation collapses PW because the abelian translation group has a **single**
resultant `q+q'`. GeoVac's SO(4)/Gaunt/6j rules are **non-abelian**: the two correlator
lines at the shared vertex couple through a genuine Clebsch–Gordan multiplet (multiplicity
up to `2l_max+1`), which does not reduce to one channel at any l_max — it grows. L3 stays
irreducibly three-body. There is no free selection-rule cusp fix. **xTC (normal-order /
1-RDM contraction) is the required route** — and the constructive silver lining is that the
xTC-reduced effective two-body operator inherits the ≈5% Gaunt angular sparsity measured
here (Paper 22 applies to the contracted operator). A full multi-center calculation is **not
warranted** for the collapse question: two centers only enlarge the resultant multiplet
(add m-mixing), never shrink it to one — the single-center measurement is decisive.
Cross-ref: CLAUDE.md §3 TC rows were solver bolt-ons; this is the structural operator-rank
statement they lacked. Paper 59 §f12 (geminals native in Fock momentum) is consistent —
the correlator is cheap, but its *three-body* transform is not abelian-separable off plane waves.
