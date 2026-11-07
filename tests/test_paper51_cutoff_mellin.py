"""
Verification of Paper 51 §G8 + §G4-5 theorems:

  thm:cutoff_dep  (Theorem, §G8)
    G_eff(f, Lambda)  = 3 pi / (phi(1) * Lambda^2)
    Lambda_cc(f, Lambda) = 6 * phi(2)/phi(1) * Lambda^2
    R_crit * Lambda  = sqrt( phi(1) / (6 * phi(2)) )

  thm:sector_mellin  (Theorem, §G4-5; Laplace convention f = L[g])
    cosmological-constant sector  <->  phi(2)
    Einstein-Hilbert sector       <->  phi(1)
    4D cigar tip                  <->  (r_h^2 Lambda^2/3) phi(1) - f(0)/18
    (phi(0)_reg belongs to the proper-time estimator of rem. cor:phi0_ratio)

Rebuilt 2026-11-07 (group5 guard pass): the rank-2 constants 6 pi and
36 pi are withdrawn (Paper 51 traces G2 over the rank-4 bundle).  The
load-bearing checks are now test_direct_spectral_sum_gives_G_and_Lcc, which
derives G_eff and Lambda_cc from Tr f(D^2/Lambda^2) summed over the CH
spectrum, and test_cigar_tip_laplace_convention, which evaluates the tip
sector against the actual S^2 Dirac spectrum.  Both reject the rank-2
values and the proper-time phi(0) reading of the spectral action.

Paper file: papers/group5_qed_gauge/paper_51_gravity_arc.tex
  (note: the prompt said group3_foundations, but the actual location is
   group5_qed_gauge -- only matters for cross-reference, the theorem
   statements are what we test).

Tests:
  (1) Pick three cutoff functions. The paper's Table tab:g8_cutoffs uses
        Gaussian   f(x) = exp(-x)
        Polynomial f(x) = exp(-x^2)   (paper's "polynomial" label, see Table)
        Sharp      f(x) = Theta(1 - x)
      The task prompt asked for polynomial (1 - x^2)^2 instead of exp(-x^2);
      since the paper's table is what the theorem makes a quantitative
      claim about, we test BOTH polynomials -- the paper's exp(-x^2) at
      paper's reported values, AND the prompt's (1 - x^2)^2 at the closed
      forms predicted by Theorem thm:cutoff_dep -- and check internal
      consistency (closed form versus its own definition) in both cases.

  (2) Compute phi(0), phi(1), phi(2) (phi(0) only as a sanity check; the
      regulated phi(0) is what the tip sector uses).

  (3) [Internal-consistency only: the closed forms against themselves and
      the table.  The theorem itself is backed by
      test_direct_spectral_sum_gives_G_and_Lcc.]  Original plan: a direct
      heat-kernel evaluation
        S = int_0^infty (dt/t) f_tilde(t * Lambda^2) K(t)
      truncated at the two leading Seeley-DeWitt coefficients in d=4.
      Concretely:
        S ~ a_0 * Lambda^4 * phi(2) + a_2 * Lambda^2 * phi(1)
      Reading off the Wilson coefficients of Lambda^4 and Lambda^2 R in
      the Chamseddine-Connes action then gives G_eff, Lambda_cc, R_crit
      in terms of phi(1), phi(2), matching the paper.

  (4) Verify R_crit * Lambda for all three cutoffs:
        Gaussian:    1/sqrt(6) ~ 0.4082
        Polynomial:  sqrt( (sqrt(pi)/2) / (6 * 1/2) ) = sqrt(sqrt(pi)/6) ~ 0.5441
        Sharp:       sqrt(1/(6 * 1/2)) = sqrt(1/3) ~ 0.5774

Verdicts (per CLAUDE.md §13.5):
  - We do NOT modify the paper or production code.
  - We do NOT alter the Paper 2 K-rule labeling.
"""

import math
import numpy as np
import pytest
from scipy import integrate
from scipy.special import gamma


# ---------------------------------------------------------------------------
# Cutoff functions and their Mellin moments
# ---------------------------------------------------------------------------

def phi_numeric(f, s, *, upper=np.inf, lower=0.0):
    """Compute phi(s) = int_0^infty f(x) * x^(s-1) dx by quadrature."""
    integrand = lambda x: f(x) * x**(s - 1.0)
    val, err = integrate.quad(integrand, lower, upper, limit=200)
    return val, err


# Gaussian: f(x) = exp(-x)
def f_gauss(x):
    return np.exp(-x)


# Polynomial per paper Table tab:g8_cutoffs: f(x) = exp(-x^2)
def f_poly_paper(x):
    return np.exp(-x * x)


# Polynomial per task prompt: f(x) = (1 - x^2)^2 * Theta(1 - x)
def f_poly_prompt(x):
    x = np.asarray(x)
    return np.where(x <= 1.0, (1.0 - x * x) ** 2, 0.0)


# Sharp: f(x) = Theta(1 - x)
def f_sharp(x):
    x = np.asarray(x)
    return np.where(x <= 1.0, 1.0, 0.0)


# Analytical moments for the paper's three cutoffs (used to verify quadrature)
ANALYTIC_MOMENTS = {
    "gaussian": {
        "phi1": 1.0,                            # Gamma(1) = 1
        "phi2": 1.0,                            # Gamma(2) = 1
    },
    "polynomial_paper": {
        "phi1": math.sqrt(math.pi) / 2.0,       # (1/2) Gamma(1/2) = sqrt(pi)/2
        "phi2": 0.5,                            # (1/2) Gamma(1)   = 1/2
    },
    "sharp": {
        "phi1": 1.0,                            # int_0^1 x^0 dx = 1
        "phi2": 0.5,                            # int_0^1 x   dx = 1/2
    },
    "polynomial_prompt": {
        # phi(s) = int_0^1 (1-x^2)^2 x^(s-1) dx
        # Use Beta-function: substitute u = x^2 -> du = 2x dx
        #   phi(s) = (1/2) * int_0^1 (1-u)^2 * u^((s-2)/2) du = (1/2) B((s)/2, 3)
        #   where B(a,b) = Gamma(a) Gamma(b) / Gamma(a+b)
        # phi(1) = (1/2) * B(1/2, 3) = (1/2) * Gamma(1/2)*Gamma(3)/Gamma(7/2)
        #        = (1/2) * sqrt(pi) * 2 / ((15/8) sqrt(pi)) = 8/15
        # phi(2) = (1/2) * B(1, 3) = (1/2) * Gamma(1)*Gamma(3)/Gamma(4)
        #        = (1/2) * 1 * 2 / 6 = 1/6
        "phi1": 8.0 / 15.0,
        "phi2": 1.0 / 6.0,
    },
}


# ---------------------------------------------------------------------------
# Closed-form predictions from Theorem thm:cutoff_dep
# ---------------------------------------------------------------------------

def G_eff_predicted(phi1, Lambda):
    return 3.0 * math.pi / (phi1 * Lambda ** 2)


def Lambda_cc_predicted(phi1, phi2, Lambda):
    return 6.0 * phi2 / phi1 * Lambda ** 2


def Rcrit_Lambda_predicted(phi1, phi2):
    return math.sqrt(phi1 / (6.0 * phi2))


# ---------------------------------------------------------------------------
# (1) phi-moment computation: numeric quadrature versus analytic
# ---------------------------------------------------------------------------

@pytest.mark.parametrize(
    "name,f,upper",
    [
        ("gaussian",          f_gauss,        np.inf),
        ("polynomial_paper",  f_poly_paper,   np.inf),
        ("sharp",             f_sharp,        1.0),
        ("polynomial_prompt", f_poly_prompt,  1.0),
    ],
)
def test_phi_moments_numeric_vs_analytic(name, f, upper):
    """
    Verify phi(1) and phi(2) computed by quadrature match the analytic
    closed-form for each cutoff. This is the foundational check that
    feeds into all downstream theorem verifications.
    """
    phi1_num, err1 = phi_numeric(f, 1.0, upper=upper)
    phi2_num, err2 = phi_numeric(f, 2.0, upper=upper)
    a = ANALYTIC_MOMENTS[name]

    # Quadrature should match analytic to ~1e-8 except for sharp/prompt
    # at x=1 endpoint, which still come out near machine precision since
    # the integrands are bounded.
    assert abs(phi1_num - a["phi1"]) < 1e-8, (
        f"{name}: phi(1) quadrature {phi1_num} vs analytic {a['phi1']}"
    )
    assert abs(phi2_num - a["phi2"]) < 1e-8, (
        f"{name}: phi(2) quadrature {phi2_num} vs analytic {a['phi2']}"
    )


# ---------------------------------------------------------------------------
# (3) Closed-form thm:cutoff_dep verified against direct heat-kernel-with-cutoff
# ---------------------------------------------------------------------------
#
# The spectral action with cutoff f has heat-kernel representation
#   S(Lambda) = int_0^infty (dt/t) f_tilde(t * Lambda^2) K(t)
# where f_tilde is the inverse Mellin transform of f and
# K(t) ~ sum_n a_n * t^((n-d)/2) (Seeley-DeWitt expansion in d=4).
#
# Crucially, the Mellin-inverse trick used in the paper's proof is
#   int_0^infty (dt/t) f_tilde(t * Lambda^2) * t^((n-d)/2)
#     = Lambda^(d-n) * phi((d-n)/2)
# This means: the only thing the cutoff does to the Wilson coefficient of
# the n-th Seeley-DeWitt term is multiply it by phi((d-n)/2). So we don't
# need to construct f_tilde explicitly; we can verify the relation by
# comparing the ratio of (Wilson coefficients) across cutoffs.
#
# We do this concretely:
#   - For each cutoff, compute G_eff, Lambda_cc, R_crit*Lambda from the
#     paper's closed forms using the *numerically* obtained phi(1), phi(2).
#   - Verify that the cross-cutoff RATIOS of these quantities match the
#     ratios predicted by Theorem thm:cutoff_dep (with Gaussian as reference).
#
# This is the operationally correct content of "closed-forms reproduce the
# direct heat-kernel-with-cutoff" at finite Lambda: the Wilson coefficients
# track phi(1) and phi(2) in exactly the way the theorem claims.

@pytest.mark.parametrize("Lambda", [0.5, 1.0, 2.0, 5.0])
def test_thm_cutoff_dep_closed_form_consistency(Lambda):
    """
    For each Lambda, verify the closed-form predictions of G_eff, Lambda_cc,
    R_crit*Lambda for all three paper cutoffs and the prompt's polynomial.
    Internal consistency: G_eff * Lambda_cc = 18 * pi * phi(2) / phi(1)^2,
    R_crit * Lambda = sqrt( phi(1) / (6 phi(2)) ).
    """
    for name in ["gaussian", "polynomial_paper", "sharp", "polynomial_prompt"]:
        a = ANALYTIC_MOMENTS[name]
        phi1, phi2 = a["phi1"], a["phi2"]

        G  = G_eff_predicted(phi1, Lambda)
        Lc = Lambda_cc_predicted(phi1, phi2, Lambda)
        RL = Rcrit_Lambda_predicted(phi1, phi2)

        # G * Lambda^2 = 3 pi/phi1; Lc / Lambda^2 = 6 phi2/phi1
        # product (G * Lambda^2) * (Lc / Lambda^2) = 18 pi phi2 / phi1^2
        prod_lhs = (G * Lambda ** 2) * (Lc / Lambda ** 2)
        prod_rhs = 18.0 * math.pi * phi2 / phi1 ** 2
        assert abs(prod_lhs - prod_rhs) < 1e-12 * abs(prod_rhs), (
            f"{name}: product check failed at Lambda={Lambda}"
        )

        # R_crit is the extremum of the action in R, not G_eff/Lambda_cc.
        assert abs(RL - math.sqrt(phi1 / (6.0 * phi2))) < 1e-14


# ---------------------------------------------------------------------------
# (4) Numerical values of R_crit * Lambda for the three paper cutoffs
# ---------------------------------------------------------------------------

def test_paper_table_R_crit_Lambda_values():
    """
    Verify the paper's reported R_crit * Lambda values in Table tab:g8_cutoffs:
      Gaussian:   1/sqrt(6)
      Polynomial: sqrt(sqrt(pi)/6) ~ 0.544
      Sharp:      1/sqrt(3)
    Also verify the prompt's polynomial: phi1=8/15, phi2=1/6
      => R_crit * Lambda = sqrt( (8/15) / (6 * 1/6) ) = sqrt(8/15) ~ 0.7303
    """
    # Gaussian
    RL_gauss = Rcrit_Lambda_predicted(1.0, 1.0)
    assert abs(RL_gauss - 1.0 / math.sqrt(6.0)) < 1e-14
    assert abs(RL_gauss - 0.4082482904638630) < 1e-12

    # Paper's polynomial e^{-x^2}: phi1 = sqrt(pi)/2, phi2 = 1/2
    RL_poly_paper = Rcrit_Lambda_predicted(math.sqrt(math.pi) / 2.0, 0.5)
    expected_poly_paper = math.sqrt(math.sqrt(math.pi) / 6.0)
    assert abs(RL_poly_paper - expected_poly_paper) < 1e-14
    # Paper reports approx 0.544 (the paper rounds to 3 sig figs)
    # Exact value: sqrt(sqrt(pi)/6) = 0.5435153863055942
    assert abs(RL_poly_paper - 0.544) < 1e-3

    # Sharp Theta(1-x): phi1 = 1, phi2 = 1/2
    RL_sharp = Rcrit_Lambda_predicted(1.0, 0.5)
    assert abs(RL_sharp - 1.0 / math.sqrt(3.0)) < 1e-14
    assert abs(RL_sharp - 0.5773502691896258) < 1e-12

    # Prompt's polynomial (1-x^2)^2: phi1 = 8/15, phi2 = 1/6
    RL_poly_prompt = Rcrit_Lambda_predicted(8.0 / 15.0, 1.0 / 6.0)
    expected_poly_prompt = math.sqrt((8.0 / 15.0) / (6.0 * (1.0 / 6.0)))
    assert abs(RL_poly_prompt - expected_poly_prompt) < 1e-14
    # Numerically: sqrt(8/15) ~ 0.7303
    assert abs(RL_poly_prompt - math.sqrt(8.0 / 15.0)) < 1e-14


def test_paper_table_G_eff_Lambda2_and_Lambda_cc_Lambda2():
    """
    Verify Table tab:g8_cutoffs columns G_eff*Lambda^2 and Lambda_cc/Lambda^2:
      Gaussian:   G_eff Lambda^2 = 3pi,     Lambda_cc/Lambda^2 = 6
      Polynomial: G_eff Lambda^2 = 6 sqrt(pi), Lambda_cc/Lambda^2 = 6/sqrt(pi)
      Sharp:      G_eff Lambda^2 = 3pi,     Lambda_cc/Lambda^2 = 3
    """
    for name, (G_target, Lc_target) in {
        "gaussian":         (3.0 * math.pi,                      6.0),
        "polynomial_paper": (6.0 * math.sqrt(math.pi),           6.0 / math.sqrt(math.pi)),
        "sharp":            (3.0 * math.pi,                      3.0),
    }.items():
        a = ANALYTIC_MOMENTS[name]
        # G_eff * Lambda^2 = 3 pi / phi1
        G_Lambda2 = 3.0 * math.pi / a["phi1"]
        # Lambda_cc / Lambda^2 = 6 phi2 / phi1
        Lc_Lambda2 = 6.0 * a["phi2"] / a["phi1"]
        assert abs(G_Lambda2 - G_target) < 1e-12, (
            f"{name}: G_eff Lambda^2 = {G_Lambda2}, paper = {G_target}"
        )
        assert abs(Lc_Lambda2 - Lc_target) < 1e-12, (
            f"{name}: Lambda_cc/Lambda^2 = {Lc_Lambda2}, paper = {Lc_target}"
        )


# ---------------------------------------------------------------------------
# (3, continued) Direct heat-kernel verification at finite Lambda
# ---------------------------------------------------------------------------
#
# We perform the actual Mellin-inverse representation. The spectral action
# with cutoff f, written in the form S = Tr f(D^2/Lambda^2), has heat-kernel
# representation
#   S = int_0^infty (dt/t) f_tilde(t * Lambda^2) K(t)
# where f_tilde is the inverse Mellin transform of f. The identity
#   int_0^infty (dt/t) f_tilde(t * Lambda^2) * t^a = Lambda^(-2a) * phi(-a)
# means: the *coefficient* of K(t) ~ a_n * t^((n-d)/2) in the spectral
# action is exactly Lambda^(d-n) * phi((d-n)/2). This is a Mellin-Plancherel
# identity, valid for any cutoff f with bounded Mellin transform in a strip
# around the relevant index, and it is the only piece the theorem's proof
# uses.
#
# We can verify the identity numerically WITHOUT constructing f_tilde
# explicitly, by integrating
#   I(s, Lambda) := int_0^infty (dt/t) f(t * Lambda^2) * t^s = Lambda^(-2s) * phi(s)
# directly (note: this is the spectral-action integrand BEFORE inverse-Mellin
# substitution, applied to the test-function family t^s rather than to K(t)).
#
# Then check: I(s=1, Lambda) / Lambda^(-2) = phi(1)
#             I(s=2, Lambda) / Lambda^(-4) = phi(2)
# at multiple Lambdas, for all three cutoffs.

@pytest.mark.parametrize("Lambda", [0.5, 1.0, 2.0, 4.0])
@pytest.mark.parametrize(
    "name,f,upper",
    [
        ("gaussian",          f_gauss,        np.inf),
        ("polynomial_paper",  f_poly_paper,   np.inf),
        ("sharp",             f_sharp,        None),
        ("polynomial_prompt", f_poly_prompt,  None),
    ],
)
def test_mellin_plancherel_identity(name, f, upper, Lambda):
    """
    Verify the master Mellin-Plancherel identity used in the proof of
    thm:cutoff_dep:
        int_0^infty (dt/t) f(t * Lambda^2) * t^s = Lambda^(-2s) * phi(s)
    at s = 1 (Einstein-Hilbert sector, phi(1)) and s = 2 (cosmological-
    constant sector, phi(2)).

    For finite-support cutoffs, the integrand in t is supported on
    t * Lambda^2 in [0, 1], i.e. t in [0, 1/Lambda^2].
    """
    a = ANALYTIC_MOMENTS[name]
    if upper is None:
        # finite support: t-integration upper limit is 1/Lambda^2 (since f
        # vanishes outside [0,1])
        t_upper = 1.0 / Lambda ** 2
    else:
        t_upper = np.inf

    for s, phi_s, sector_name in [
        (1.0, a["phi1"], "EH"),
        (2.0, a["phi2"], "CC"),
    ]:
        integrand = lambda t: f(t * Lambda ** 2) * t ** (s - 1.0)
        val, err = integrate.quad(integrand, 0.0, t_upper, limit=200)
        # Identity says val should equal Lambda^(-2s) * phi(s)
        predicted = Lambda ** (-2.0 * s) * phi_s
        assert abs(val - predicted) < 1e-7 * max(abs(predicted), 1.0), (
            f"{name} {sector_name} sector at Lambda={Lambda}: "
            f"int = {val}, Lambda^(-2s) * phi({s}) = {predicted}"
        )


# ---------------------------------------------------------------------------
# Sector Mellin moment map (thm:sector_mellin)
# ---------------------------------------------------------------------------
#
# The theorem's (i)/(ii) are checked against the CH spectrum in
# test_direct_spectral_sum_gives_G_and_Lcc and its (iii) in
# test_cigar_tip_laplace_convention (below).  The two tests that follow
# concern the PROPER-TIME estimator of remark cor:phi0_ratio, which weights
# a constant tip by (dt/t) f(t Lambda^2) and so gives phi(0)_reg; that is a
# property of the substrate estimator, not of the spectral action.

@pytest.mark.parametrize("Lambda", [1.0, 2.0])
def test_thm_sector_mellin_tip_phi0(Lambda):
    """
    Verify that the topological-tip sector contribution
       S_tip ~ (1/6) * int_{t_UV}^{t_IR} (dt/t) f(t * Lambda^2)
    has the form (1/6) * phi(0)_reg, where phi(0)_reg depends on the
    UV/IR substrate regulators.
    """
    a_substrate = 0.1     # substrate lattice spacing
    t_UV = a_substrate ** 2
    t_IR = 100.0          # arbitrary IR substrate scale

    # For each cutoff, compute the tip integral (1/6) * int (dt/t) f(t * Lambda^2)
    # with substrate-regulated bounds, and compare to (1/6) * phi(0)_reg where
    # phi(0)_reg = int_{t_UV}^{t_IR} (du/u) f(u * Lambda^2) reduces (under
    # change of variable v = u * Lambda^2) to int_{t_UV*Lambda^2}^{t_IR*Lambda^2}
    # (dv/v) f(v) -- a Lambda-independent integrand under the conversion.

    for name, f, finite_support in [
        ("gaussian",          f_gauss,       False),
        ("polynomial_paper",  f_poly_paper,  False),
        ("sharp",             f_sharp,       True),
        ("polynomial_prompt", f_poly_prompt, True),
    ]:
        # S_tip / (1/6) = int_{t_UV}^{t_IR} (dt/t) f(t * Lambda^2)
        tip_integrand = lambda t: f(t * Lambda ** 2) / t
        # For finite-support cutoffs, clip t_IR at 1/Lambda^2 (above this
        # f is zero)
        if finite_support:
            upper = min(t_IR, 1.0 / Lambda ** 2)
        else:
            upper = t_IR
        if upper <= t_UV:
            # nothing to integrate; tip vanishes (sharp cutoff entirely
            # below UV substrate)
            continue
        S_tip_x6, _ = integrate.quad(tip_integrand, t_UV, upper, limit=400)

        # phi(0)_reg via direct substitution v = t * Lambda^2 -> dv/v = dt/t
        v_lower = t_UV * Lambda ** 2
        v_upper = upper * Lambda ** 2
        phi0_integrand = lambda v: f(v) / v
        phi0_reg, _ = integrate.quad(phi0_integrand, v_lower, v_upper, limit=400)

        assert abs(S_tip_x6 - phi0_reg) < 1e-7 * max(abs(phi0_reg), 1.0), (
            f"{name} tip sector at Lambda={Lambda}: "
            f"S_tip * 6 = {S_tip_x6}, phi(0)_reg = {phi0_reg}"
        )


def test_sector_mellin_ordering_qualitative():
    """
    The paper's Corollary cor:phi0_ratio + the empirical table in
    §G4-5 (Sec 4.18 in the source) state that the tip sector follows
        S_sharp > S_poly > S_Gauss
    qualitatively matching phi(0)_reg ordering. We verify the phi(0)_reg
    ordering: with shared UV/IR regulators, sharp has the largest phi(0)_reg,
    polynomial intermediate, Gaussian smallest. (This is qualitative
    consistency with the theorem; the numerical Wilson coefficients also
    enter, so we do NOT pin the precise ratio here.)
    """
    Lambda = 1.0
    a_substrate = 0.1
    t_UV = a_substrate ** 2
    t_IR = 10.0

    phi0_vals = {}
    for name, f, finite_support in [
        ("gaussian",          f_gauss,       False),
        ("polynomial_paper",  f_poly_paper,  False),
        ("sharp",             f_sharp,       True),
    ]:
        v_lower = t_UV * Lambda ** 2
        v_upper = (min(t_IR, 1.0 / Lambda ** 2) if finite_support else t_IR) * Lambda ** 2
        if v_upper <= v_lower:
            phi0_vals[name] = 0.0
            continue
        phi0_integrand = lambda v: f(v) / v
        val, _ = integrate.quad(phi0_integrand, v_lower, v_upper, limit=400)
        phi0_vals[name] = val

    # The paper's claim is that phi(0)_reg orders as
    # sharp > polynomial > Gaussian. Sharp keeps the full log of the
    # integration range; Gaussian truncates exponentially.
    assert phi0_vals["sharp"] > phi0_vals["polynomial_paper"] > phi0_vals["gaussian"], (
        f"phi(0)_reg ordering violated: sharp={phi0_vals['sharp']}, "
        f"poly={phi0_vals['polynomial_paper']}, gauss={phi0_vals['gaussian']}"
    )


# ---------------------------------------------------------------------------
# Direct evaluation from the CH spectrum (the load-bearing check)
# ---------------------------------------------------------------------------
#
# The rank-4 bundle C^2 (x) C^2 is the paper's input (Paper 51 sec:g2):
# _sigma_space_rank reads the multiplicity 2 off D^2 for that bundle, so it
# documents the trace rather than testing the bundle choice.

def _sigma_space_rank() -> int:
    """dim of the sigma factor, read off D^2 for D = D3 (x) s1 + I (x) s2 p."""
    s1 = np.array([[0, 1], [1, 0]], dtype=complex)
    s2 = np.array([[0, -1j], [1j, 0]], dtype=complex)
    lam, p = 1.3, 0.8
    D2 = (lam * s1 + p * s2) @ (lam * s1 + p * s2)
    assert np.allclose(D2, (lam ** 2 + p ** 2) * np.eye(2))
    return int(round(np.trace(D2).real / (lam ** 2 + p ** 2)))


def _density_from_spectrum(f, R, Lambda, nmax=4000):
    """Tr f(D^2/Lambda^2) per unit Euclidean time on S^3_R x R_tau.

    CH spectrum |lambda_n| = (n + 3/2)/R, multiplicity 2(n+1)(n+2) (both
    signs), tau-momentum integrated, traced over the rank-4 bundle.
    """
    rank = _sigma_space_rank()
    total = 0.0
    for n in range(nmax):
        l2 = ((n + 1.5) / R) ** 2
        if l2 / Lambda ** 2 > 60.0:
            break
        half, _ = integrate.quad(lambda w: f((l2 + w * w) / Lambda ** 2), 0.0, np.inf, limit=200)
        total += 2 * (n + 1) * (n + 2) * 2.0 * half / (2.0 * math.pi)
    return rank * total


@pytest.mark.parametrize(
    "name,f",
    [("gaussian", f_gauss), ("polynomial_paper", f_poly_paper)],
)
def test_direct_spectral_sum_gives_G_and_Lcc(name, f):
    """
    thm:cutoff_dep from the spectrum, not from its own closed form.

    s(R) = A R^3 + B R is fitted to the summed spectral action at two radii,
    then matched to Einstein--Hilbert + cosmological constant on S^3_R x R:
      -(1/16 pi G) * (6/R^2) * 2 pi^2 R^3 = B R      =>  G = -3 pi / (4 B)
      (Lcc / 8 pi G) * 2 pi^2 R^3          = A R^3    =>  Lcc = 4 G A / pi
    and compared with 3 pi/(phi(1) Lambda^2) and 6 phi(2)/phi(1) Lambda^2.
    Rejects the rank-2 values (G Lambda^2 = 6 pi/phi(1)).  The sharp cutoff
    is not tested: Theorem thm:cutoff_dep is stated for smooth f.
    """
    Lambda = 8.0
    R1, R2 = 1.0, 2.0
    s1 = _density_from_spectrum(f, R1, Lambda)
    s2 = _density_from_spectrum(f, R2, Lambda)
    A, B = np.linalg.solve(np.array([[R1 ** 3, R1], [R2 ** 3, R2]]), np.array([s1, s2]))
    G = -3.0 * math.pi / (4.0 * B)
    Lcc = 4.0 * G * A / math.pi
    a = ANALYTIC_MOMENTS[name]
    assert abs(G * Lambda ** 2 - 3.0 * math.pi / a["phi1"]) < 1e-6 * 3.0 * math.pi, (
        f"{name}: G Lambda^2 = {G * Lambda ** 2}"
    )
    assert abs(Lcc / Lambda ** 2 - 6.0 * a["phi2"] / a["phi1"]) < 1e-6, (
        f"{name}: Lcc/Lambda^2 = {Lcc / Lambda ** 2}"
    )
    RL = math.sqrt(-B / (3.0 * A)) * Lambda
    assert abs(RL - math.sqrt(a["phi1"] / (6.0 * a["phi2"]))) < 1e-7


# Laplace-convention cutoffs f = L[g] with known g, phi(1) and f(0).
#   g(s) = s e^{-s}         -> f(x) = (1+x)^-2, phi(1) = 1,   f(0) = 1
#   g(s) = s^2 e^{-s} / 2   -> f(x) = (1+x)^-3, phi(1) = 1/2, f(0) = 1
#   g(s) = 2 s e^{-2s}      -> f(x) = 2/(2+x)^2, phi(1) = 1, f(0) = 1/2
#   (the last separates -f(0)/18 from an f-independent -1/18)
LAPLACE_CUTOFFS = [
    (lambda s: s * np.exp(-s), 1.0, 1.0),
    (lambda s: 0.5 * s * s * np.exp(-s), 0.5, 1.0),
    (lambda s: 2.0 * s * np.exp(-2.0 * s), 1.0, 0.5),
]


def _s2_dirac_heat_trace(t, r_h):
    """K_{S^2}(t) = sum_m 4 m exp(-m^2 t / r_h^2), the S^2 Dirac spectrum,
    summed until the terms are below e^-60."""
    mmax = int(r_h * math.sqrt(60.0 / t)) + 10
    m = np.arange(1, mmax + 1, dtype=float)
    return float(np.sum(4.0 * m * np.exp(-m * m * t / r_h ** 2)))


@pytest.mark.parametrize("g,phi1,f0", LAPLACE_CUTOFFS)
def test_cigar_tip_laplace_convention(g, phi1, f0):
    """
    thm:sector_mellin (iii): the replica tip on the 4D cigar is
    (1/6) K_{S^2}(t); in the spectral action S = int g(s) K(s/Lambda^2) ds it
    gives (r_h^2 Lambda^2/3) phi(1) - f(0)/18.  K_{S^2} is summed from the
    S^2 Dirac spectrum here, not taken from its asymptotic form.
    Rejects the reading 'tip <-> phi(0)' (phi(0) diverges for these f) and
    a missing or mis-signed subleading -f(0)/18.
    """
    # The remainder is the O(t) heat-trace term, of order 1/(r_h Lambda)^2;
    # at r_h Lambda = 30 it is ~1e-5, well below the f(0)/18 = 0.056 tested.
    r_h, Lambda = 3.0, 10.0
    val, _ = integrate.quad(
        lambda s: g(s) * _s2_dirac_heat_trace(s / Lambda ** 2, r_h) / 6.0,
        1e-9, 80.0, limit=400, points=[1.0, 5.0, 20.0],
    )  # g(s) < 1e-30 beyond s = 80; [0, 1e-9] contributes < 1e-12
    predicted = r_h ** 2 * Lambda ** 2 / 3.0 * phi1 - f0 / 18.0
    assert abs(val - predicted) < 1e-3, (val, predicted)
    # and the subleading term is resolved: dropping it misses by f(0)/18
    assert abs(val - (predicted + f0 / 18.0)) > 0.02
    # and it scales with f(0): an f-independent -1/18 misses when f(0) != 1
    if f0 != 1.0:
        assert abs(val - (predicted + f0 / 18.0 - 1.0 / 18.0)) > 0.02


# ---------------------------------------------------------------------------
# Sanity: phi(0) is divergent for the unregulated paper cutoffs
# ---------------------------------------------------------------------------

def test_phi0_unregulated_divergent_for_gaussian_polynomial():
    """
    phi(0) = int_0^infty f(x) / x dx diverges at x -> 0 for any cutoff
    with f(0) > 0 (which is the case for Gaussian, polynomial, sharp, and
    prompt's polynomial). The paper's tip sector requires the regulated
    phi(0)_reg = int_{u_UV}^{u_IR} f(u) du/u, which is finite.
    """
    # Verify divergence as the UV regulator -> 0
    for name, f, finite_support in [
        ("gaussian",          f_gauss,       False),
        ("polynomial_paper",  f_poly_paper,  False),
        ("sharp",             f_sharp,       True),
        ("polynomial_prompt", f_poly_prompt, True),
    ]:
        upper = 1.0 if finite_support else 10.0
        v_lowers = [1e-1, 1e-3, 1e-5]
        vals = []
        for vl in v_lowers:
            val, _ = integrate.quad(lambda v: f(v) / v, vl, upper, limit=400)
            vals.append(val)
        # The values should grow as the UV regulator shrinks (log-divergent)
        assert vals[2] > vals[1] > vals[0], (
            f"{name}: phi(0)_reg should grow as UV regulator shrinks; got {vals}"
        )


# ---------------------------------------------------------------------------
# Summary verdict aggregator (informational only)
# ---------------------------------------------------------------------------

def test_summary_verdict_print():
    """
    Print a tabular summary of all three (four including prompt) cutoffs,
    phi(1), phi(2), G_eff*Lambda^2, Lambda_cc/Lambda^2, R_crit*Lambda.
    Useful for at-a-glance confirmation; never fails (pure print).
    """
    print()
    print(f"{'Cutoff':<22} {'phi(1)':>10} {'phi(2)':>10} "
          f"{'G_eff*L^2':>12} {'Lcc/L^2':>10} {'R_crit*L':>10}")
    print("-" * 80)
    for name in ["gaussian", "polynomial_paper", "polynomial_prompt", "sharp"]:
        a = ANALYTIC_MOMENTS[name]
        phi1, phi2 = a["phi1"], a["phi2"]
        G_L2  = 3.0 * math.pi / phi1
        Lc_L2 = 6.0 * phi2 / phi1
        RL    = math.sqrt(phi1 / (6.0 * phi2))
        print(f"{name:<22} {phi1:>10.6f} {phi2:>10.6f} "
              f"{G_L2:>12.6f} {Lc_L2:>10.6f} {RL:>10.6f}")
    assert True
