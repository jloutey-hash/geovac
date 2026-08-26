"""Backing tests for the Paper 32 remark: the multi-center composition reducible->irreducible
transition (LiH 2-center REDUCIBLE / Halmos vs BeH2 3-center IRREDUCIBLE full M6).

The overlap matrices are HARDCODED VALIDATED values -- LiH: exact Slater two-center overlap to
12 digits (Topos-3 engine, debug/compute_topos3_two_center_meet.py); BeH2: exact prolate-
spheroidal evaluator validated vs mpmath to 2.3e-14 (debug/beh2_ci_exact_landscape.py). The
load-bearing logic tested here is the ALGEBRA STRUCTURE (commutant dimension => reducible vs
irreducible), which is fast exact linear algebra on the validated overlaps; the overlap values
themselves are validated in the debug engines (their provenance is cited above, permanent record
in CHANGELOG v5.1.0).

Claims pinned:
  (1) commutant routine self-validates on 3 cases with known answers;
  (2) LiH (2 centers): the pair of center-projections is REDUCIBLE (commutant dim > 1), the three
      principal angles are 7.6/44.7/67.3 deg, and ||[P_A,P_B]|| = 1/2 (the 45-deg ceiling);
  (3) BeH2 (3 centers): the three center-projections generate M6 IRREDUCIBLY (commutant dim = 1
      <=> full matrix algebra, by double-commutant), while EVERY pair is REDUCIBLE.
  (4) BASIS ROBUST: enriching the sigma menu {1s,2p0}->{1s,2s,2p0} keeps the triple IRREDUCIBLE
      (M6->M9) with pairs reducible -- the categorical jump is not a minimal-basis artifact.
  (5) BENDING WELDS: in-plane {s,px,pz}x3, EXACTLY-linear is REDUCIBLE (sigma/pi decouple,
      commutant dim 2); ANY bend (H2O 104.5 deg) welds it IRREDUCIBLE M9 -- the operator-algebra
      face of Renner-Teller.
  (6) TWO INVARIANTS: irreducibility (commutant=1) is a KNIFE-EDGE generic to >=3 coupled
      projections (holds at any nonzero coupling), while the nested-commutator MAGNITUDE is the
      bonding-specific three-body signal (large at coupling, ->0 as coupling->0). Scaling the
      off-diagonal overlaps by eps demonstrates both on the SAME projectors.

Validated overlaps for the new tests are hardcoded from the same prolate-spheroidal evaluator
(fast float64, cross-checked vs mpmath topos3 to <=4e-13 incl. the 2s states; permanent record
in CHANGELOG v5.1.0, drivers debug/diag_three_center_generality.py).
"""
import numpy as np
import pytest

TOL = 1e-6


def commutant_dim(Ps):
    """dim of A' = {X : [P_i, X] = 0 for all i}, X in C^{n x n} (via the stacked commutator map)."""
    n = Ps[0].shape[0]
    rows = [np.kron(P, np.eye(n)) - np.kron(np.eye(n), P.T) for P in Ps]
    s = np.linalg.svd(np.vstack(rows), compute_uv=False)
    return int((s < TOL * max(1.0, s[0])).sum()) + (n * n - len(s))


def projectors_from_G(G, block_sizes):
    """Metric-orthogonal projectors onto consecutive orbital blocks (diagonal blocks of G = I).
    Commutant dim is whitening-invariant, so Cholesky whitening is used (results identical to
    Loewdin; the whitening only affects even/odd labelling, not the algebra structure)."""
    Xh = np.linalg.cholesky(G).T
    Ps, k = [], 0
    for b in block_sizes:
        A = Xh[:, k:k + b]
        Ps.append(A @ np.linalg.pinv(A))
        k += b
    return Ps


def cn(X):
    return float(np.linalg.norm(X, 2))


# --- validated overlap matrices ---
# LiH bond block, m=0 {1s,2s,2p0}, R_eq=3.015, both Z=1 (exact Slater, 12 digits)
S_AB_LIH = np.array([[0.3455, -0.2518, -0.4950],
                     [-0.2518, 0.7993, 0.1268],
                     [0.4950, -0.1268, 0.4786]])
# BeH2 sigma pair {1s,2p0}, Be Z=2 / H Z=1: Be-H at R=2.5, H-H raw block at R=5.0 (exact, 2.3e-14)
PAR = np.diag([1.0, -1.0])
S_BEH = np.array([[0.227626, -0.278504],
                  [0.572885, 0.262679]])
S_HH_RAW = np.array([[0.096577, -0.415672],
                     [0.415672, 0.00513]])


def _lih_pair():
    I = np.eye(3)
    G = np.block([[I, S_AB_LIH], [S_AB_LIH.T, I]])
    return projectors_from_G(G, [3, 3])


def _beh2_triple():
    S1 = S_BEH
    S2 = PAR @ S_BEH @ PAR
    SHH = PAR @ S_HH_RAW @ PAR
    I = np.eye(2)
    G = np.block([[I, S1, S2], [S1.T, I, SHH], [S2.T, SHH.T, I]])
    assert np.linalg.eigvalsh(G).min() > 1e-6, "overlap metric not positive-definite"
    return projectors_from_G(G, [2, 2, 2])


# ============================== tests ==============================
def test_commutant_selfcheck():
    """The commutant routine on three cases with known answers (self-validation)."""
    rng = np.random.default_rng(1)
    Q, _ = np.linalg.qr(rng.standard_normal((6, 2)))
    assert commutant_dim([Q @ Q.T]) == 20                       # single rank-2 in 6-dim: 4+16
    assert commutant_dim([np.diag([1., 1, 0, 0, 0, 0]),
                          np.diag([1., 0, 1, 0, 0, 0])]) > 1     # commuting pair: reducible
    c, s = np.cos(0.7), np.sin(0.7)
    assert commutant_dim([np.array([[1., 0], [0, 0]]),
                          np.array([[c * c, c * s], [c * s, s * s]])]) == 1  # generic pair: irreducible


def test_lih_pair_reducible_halmos():
    """LiH: the two center-projections are REDUCIBLE (Halmos, dim A' > 1); 3 principal angles; ceiling."""
    PA, PB = _lih_pair()
    assert commutant_dim([PA, PB]) > 1                          # reducible (decomposes into <=2x2)
    ang = np.degrees(np.arccos(np.clip(np.linalg.svd(S_AB_LIH, compute_uv=False), 0, 1)))
    assert np.allclose(np.sort(ang), [7.6, 44.7, 67.3], atol=0.2)   # the three principal angles
    assert abs(cn(PA @ PB - PB @ PA) - 0.5) < 1e-2             # ||[P,P]|| at the 45-deg ceiling (max)


def test_beh2_triple_irreducible_pairs_reducible():
    """BeH2: the three center-projections are IRREDUCIBLE (dim A' = 1 <=> full M6); every pair reducible."""
    PBe, PH1, PH2 = _beh2_triple()
    assert commutant_dim([PBe, PH1, PH2]) == 1                 # IRREDUCIBLE => A = full M6
    for pair in [(PBe, PH1), (PBe, PH2), (PH1, PH2)]:
        assert commutant_dim(list(pair)) > 1                   # each PAIR reducible


def test_reducible_irreducible_is_a_real_contrast():
    """Non-tautology guard: the SAME routine gives >1 for pairs and exactly 1 for the triple, on the
    SAME BeH2 projections -- so 'irreducible' is a measured property of adding the third center."""
    PBe, PH1, PH2 = _beh2_triple()
    pair_dims = [commutant_dim([PBe, PH1]), commutant_dim([PBe, PH2]), commutant_dim([PH1, PH2])]
    triple_dim = commutant_dim([PBe, PH1, PH2])
    assert min(pair_dims) > 1 and triple_dim == 1 and triple_dim < min(pair_dims)


# ============================== (4) basis robustness ==============================
# enriched sigma {1s,2s,2p0}, Be Z=2 / H Z=1: Be-H at R=2.5, H-H raw block at R=5.0
S_BEH_9 = np.array([[0.227626, -0.066271, -0.278504],
                    [-0.528422, 0.474605, 0.448796],
                    [0.572885, 0.151987, 0.262679]])
S_HH_9 = np.array([[0.096577, -0.270245, -0.415672],
                   [-0.270245, 0.672071, 0.299268],
                   [0.415672, -0.299268, 0.00513]])
PAR9 = np.diag([1.0, 1.0, -1.0])   # z->-z parity of {1s(even),2s(even),2p0(odd)}


def _beh2_triple_enriched():
    S1 = S_BEH_9
    S2 = PAR9 @ S_BEH_9 @ PAR9
    SHH = PAR9 @ S_HH_9 @ PAR9
    I = np.eye(3)
    G = np.block([[I, S1, S2], [S1.T, I, SHH], [S2.T, SHH.T, I]])
    assert np.linalg.eigvalsh(G).min() > 1e-6
    return projectors_from_G(G, [3, 3, 3])


def test_beh2_basis_robust_irreducible():
    """Enriching {1s,2p0}->{1s,2s,2p0} keeps the triple IRREDUCIBLE (M6->M9), pairs reducible."""
    PBe, PH1, PH2 = _beh2_triple_enriched()
    assert PBe.shape[0] == 9
    assert commutant_dim([PBe, PH1, PH2]) == 1                 # still IRREDUCIBLE (full M9)
    for pair in [(PBe, PH1), (PBe, PH2), (PH1, PH2)]:
        assert commutant_dim(list(pair)) > 1                   # pairs still reducible


# ============================== (5) bending welds ==============================
# Slater-Koster fundamentals (ss, SP=<s|pz'>, PS=<pz'|s>, ppsigma, pppi), validated evaluator.
FUNDS_BEH_25 = (0.227626, -0.278504, 0.572885, 0.262679, 0.559775)   # apex-H, R=2.5, Z2/1
FUNDS_HH_50 = (0.096577, -0.415672, 0.415672, 0.00513, 0.578015)     # H-H, R=5.0   (linear 180)
FUNDS_HH_3954 = (0.19495, -0.479592, 0.479592, 0.236805, 0.700114)   # H-H, R=3.954 (bent 104.5)


def _sk_block(funds, u):
    """lab-frame {s,px,pz} overlap block for an in-plane unit bond vector u=(ux,uz)."""
    ss, SP, PS, pps, ppp = funds
    Mloc = np.array([[ss, 0.0, SP], [0.0, ppp, 0.0], [PS, 0.0, pps]])
    ux, uz = u
    C = np.array([[1.0, 0.0, 0.0], [0.0, uz, ux], [0.0, -ux, uz]])
    return C @ Mloc @ C.T


def _triatomic_inplane(funds_AH, funds_HH, half_deg):
    a = np.radians(half_deg)
    S_A1 = _sk_block(funds_AH, (np.sin(a), np.cos(a)))
    S_A2 = _sk_block(funds_AH, (np.sin(a), -np.cos(a)))
    S_HH = _sk_block(funds_HH, (0.0, -1.0))
    I = np.eye(3)
    G = np.block([[I, S_A1, S_A2], [S_A1.T, I, S_HH], [S_A2.T, S_HH.T, I]])
    assert np.linalg.eigvalsh(G).min() > 1e-6
    return projectors_from_G(G, [3, 3, 3])


def test_bending_welds_irreducible():
    """In-plane {s,px,pz}x3: EXACTLY-linear REDUCIBLE (sigma/pi split, dimA'=2); a BEND welds M9."""
    lin = _triatomic_inplane(FUNDS_BEH_25, FUNDS_HH_50, 0.0)         # 180 deg
    bent = _triatomic_inplane(FUNDS_BEH_25, FUNDS_HH_3954, 37.75)    # 104.5 deg (H2O)
    assert commutant_dim(lin) == 2                                   # sigma + pi decouple => reducible
    assert commutant_dim(bent) == 1                                  # bend welds => IRREDUCIBLE M9


# ============================== (6) two invariants ==============================
def _nested_max(Ps):
    def nst(A, B, C):
        c = A @ B - B @ A
        return cn(c @ C - C @ c)
    return max(nst(Ps[0], Ps[1], Ps[2]), nst(Ps[1], Ps[2], Ps[0]), nst(Ps[0], Ps[2], Ps[1]))


def _beh2_triple_scaled(eps):
    """Minimal BeH2 triple with off-diagonal (inter-center) overlaps scaled by eps in [0,1]."""
    S1, S2, SHH = S_BEH, PAR @ S_BEH @ PAR, PAR @ S_HH_RAW @ PAR
    I = np.eye(2)
    G = np.block([[I, eps * S1, eps * S2], [eps * S1.T, I, eps * SHH], [eps * S2.T, eps * SHH.T, I]])
    return projectors_from_G(G, [2, 2, 2])


def test_two_invariants_coupling_vs_category():
    """Irreducibility is generic to coupled triples (knife-edge); the nested-commutator MAGNITUDE
    is the bonding signal. On the SAME family, full coupling => irreducible + large 3-body norm;
    coupling below resolution => reducible + ~0 norm."""
    Ps_full = _beh2_triple_scaled(1.0)
    Ps_off = _beh2_triple_scaled(1e-9)
    assert commutant_dim(Ps_full) == 1 and _nested_max(Ps_full) > 0.2     # irreducible + strong 3-body
    assert commutant_dim(Ps_off) > 1 and _nested_max(Ps_off) < 1e-6       # reducible + no 3-body


# ------------- non-abelian != wild: compact-group coupling is tame (Peter-Weyl) -------------
def _su2_generators(js):
    """Block-diagonal J_x, J_y, J_z for the direct sum of spin-j irreps listed in `js`."""
    from scipy.linalg import block_diag
    bx, by, bz = [], [], []
    for j in js:
        d = int(round(2 * j)) + 1
        m = np.array([j - k for k in range(d)])
        Jp = np.zeros((d, d), complex)
        off = np.sqrt(j * (j + 1) - m[1:] * (m[1:] + 1))
        for k in range(d - 1):
            Jp[k, k + 1] = off[k]
        Jm = Jp.conj().T
        bx.append((Jp + Jm) / 2); by.append((Jp - Jm) / (2j)); bz.append(np.diag(m).astype(complex))
    return [block_diag(*b) for b in (bx, by, bz)]


def _max_noncommutativity(As):
    return max(cn(As[i] @ As[j] - As[j] @ As[i])
               for i in range(len(As)) for j in range(i + 1, len(As)))


@pytest.mark.parametrize("js,expected", [([0, 1, 2], 3),      # three irreps, mult 1 each -> 1+1+1
                                         ([1, 1], 4),          # one irrep at mult 2 -> 2^2
                                         ([1, 2, 2], 5),       # 1^2 + 2^2
                                         ([1], 1)])            # single irrep -> irreducible
def test_compact_group_coupling_is_tame(js, expected):
    """SU(2)/Gaunt coupling is compact-group rep theory => type I (Peter-Weyl) => REDUCIBLE
    whenever more than one irrep is present, however non-abelian. Commutant dim = sum m_lambda^2,
    which is a self-validating prediction (Schur)."""
    As = _su2_generators(js)
    assert _max_noncommutativity(As) > 0.9, "generators must genuinely not commute"
    assert commutant_dim(As) == expected


def test_nonabelian_and_wild_are_independent():
    """The discriminator: SU(2) coupling is MORE non-abelian than the molecular projections yet
    REDUCIBLE, while the three center-projections are IRREDUCIBLE at smaller non-commutativity.
    So the TC/Gaunt non-abelianness is not the multi-center wildness (different axes)."""
    As = _su2_generators([0, 1, 2])
    Ps = _beh2_triple()
    na_group, na_proj = _max_noncommutativity(As), _max_noncommutativity(Ps)

    assert na_group > na_proj, "group case should be the MORE non-abelian one"
    assert commutant_dim(As) > 1, "compact-group algebra must stay reducible"
    assert commutant_dim(Ps) == 1, "three center-projections must be irreducible"
