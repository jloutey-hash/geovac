"""
Paper 60 -- can the isoenergetic-secular-matrix 1-norm sublinearity (atoms, eq:sublinear)
be measured for a MOLECULE?  This driver DIAGNOSES why the naive assembly fails and states
precisely what a faithful molecular measurement requires.  (Result: OBSTRUCTION identified,
not a number -- the atomic sublinearity is a single-center property.)

WHAT THE ATOMIC SUBLINEARITY RIDES ON (eq:secular / eq:sublinear).
  The atomic secular matrix is  M = diag(Z R_nu) + T',  M B = p_kappa B,  E = -1/2 p_kappa^2.
  Two single-center facts make it a clean LINEAR, p_kappa-INDEPENDENT eigenproblem whose
  1-norm is sublinear:
    (i)  the Goscinskian orbitals scale WITH p_kappa (charge Q_nu = p_kappa/R_nu), so the
         p_kappa in <1/r> ~ Q_nu cancels the -1/p_kappa prefactor -> M is p_kappa-independent;
    (ii) on ONE center <1/r> = Q/n^2 collapses the collective scale to the root-sum-of-squares
         T0_bare = sqrt(sum_p k_p^2) = Z*sqrt(sum 1/n^2) = Z R_nu  (a clean diagonal).

THE OBSTRUCTION FOR MOLECULES.
  * Collective scale is root-sum-of-squares, NOT a plain sum.  Verified below: for He 1s^2,
    sqrt(k^2+k^2) = 2 sqrt2 = Z R_nu gives the exact bare E = -4; a plain sum gives -8.
  * A fixed-scale molecular-Sturmian shortcut (all orbitals at one k, T0 = sum_p lambda_p)
    therefore FAILS the non-interacting limit: it returns p_kappa = 2 k_sg instead of the
    correct sqrt2 k_sg (E off by 2x).  Demonstrated below.
  * On two centers the nuclear attraction <1/r_A>+<1/r_B> acquires OFF-CENTER terms that
    depend on Q_nu * R (dimensionless), so T0 is neither diagonal nor p_kappa-independent:
    the matrix is p_kappa-dependent (via Q*R) and must be solved self-consistently -- exactly
    the k*R nonlinearity already seen in the one-electron Shibuya-Wulfman case.

FAITHFUL MEASUREMENT REQUIRES the p_kappa-scaled molecular Goscinskian: 2-electron configs of
atomic Sturmians on the two centers at config charge Q_nu = p_kappa/R_nu, with two-center
nuclear-attraction T0 and two-center ERIs T' at those config-specific scales, solved
self-consistently.  VALIDATION GATE: R->infinity must give 2x the H-atom energy (-1.0 Ha);
R->0 must give the He atom.  That build is the defined next step; it is NOT the fixed-scale
shortcut, which this driver shows is structurally wrong.

Diagnostic only.
"""
import warnings
warnings.filterwarnings("ignore")
import numpy as np

if __name__ == "__main__":
    np.set_printoptions(precision=5, suppress=True)

    print("=" * 78)
    print("(1) The collective isoenergetic scale is ROOT-sum-of-squares, = Z R_nu (single center)")
    print("=" * 78)
    # He 1s^2: k_p = Z/n = 2 for each 1s.  Bare (no ee) E must be -4 (two He+ 1s), paper sec:atomic.
    Z = 2.0
    kp = np.array([Z / 1, Z / 1])                 # two 1s at charge Z
    Rnu = np.sqrt(np.sum(1.0 / np.array([1, 1]) ** 2))
    pk_rss = np.sqrt(np.sum(kp ** 2))             # correct collective scale
    pk_sum = np.sum(kp)                           # naive sum (the bug)
    print(f"  He 1s^2:  Z R_nu = {Z*Rnu:.5f}   sqrt(sum k^2) = {pk_rss:.5f}   (equal -> the diagonal T0)")
    print(f"    bare E:  -1/2 (sqrt sum k^2)^2 = {-0.5*pk_rss**2:.4f}  (paper: -4.0)  [correct]")
    print(f"    bare E:  -1/2 (sum k)^2        = {-0.5*pk_sum**2:.4f}  (=-8, WRONG)  [the naive-sum bug]")

    print("\n" + "=" * 78)
    print("(2) A fixed-scale molecular-Sturmian shortcut FAILS the non-interacting limit")
    print("=" * 78)
    # H2 sigma_g^2, sigma_g one-electron SW scale k_sg (~1.75 at k_basis=1, R=1.4).
    ksg = 1.751
    print(f"  H2 sigma_g^2 (one-electron scale k_sg={ksg}):")
    print(f"    correct  p_kappa = sqrt2 k_sg = {np.sqrt(2)*ksg:.4f}   E = {-0.5*(np.sqrt(2)*ksg)**2:.4f}"
          f"  (= 2 x one-electron)")
    print(f"    shortcut p_kappa = sum lam    = {2*ksg:.4f}   E = {-0.5*(2*ksg)**2:.4f}"
          f"  (2x too deep -> NOT isoenergetic)")
    print("  => T0 = sum_p lambda_p is structurally wrong; the fixed-scale shortcut is not the")
    print("     isoenergetic secular matrix.  (This is why the earlier self-consistency scan")
    print("     returned p_kappa~3.5, E~-6 Ha -- garbage from the doubled diagonal.)")

    print("\n" + "=" * 78)
    print("(3) VERDICT")
    print("=" * 78)
    print("""  The atomic secular-matrix sublinearity is a SINGLE-CENTER property: it rides on the
  clean diagonal T0 = Z R_nu, which exists only because (i) the Goscinskian orbitals scale
  with p_kappa and (ii) on one center <1/r>=Q/n^2 collapses the collective scale to a root-
  sum-of-squares.  On two centers the nuclear attraction gains off-center, Q*R-dependent
  terms: T0 is non-diagonal and p_kappa-dependent, and a fixed-scale molecular-Sturmian
  shortcut is not the isoenergetic matrix (it fails the non-interacting limit above).

  A faithful molecular 1-norm therefore needs the p_kappa-scaled molecular Goscinskian
  (config-specific charge Q_nu=p_kappa/R_nu, two-center nuclear T0 + two-center ERIs, solved
  self-consistently; validated at R->inf = 2 H atoms and R->0 = He).  That is the defined
  next step.  The config-space sublinearity for molecules is OPEN -- with the obstruction now
  identified, not hand-waved: it is not that the integrals are hard (they close, Papers 58/59)
  but that the atomic diagonal-T0 structure that PRODUCED the sublinearity does not transfer.""")
    print("\nDONE.")
