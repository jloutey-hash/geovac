"""LiH R12-CI: the analytic 2x2 -> E_R12 (capstone).  Closes the loop: assemble the fully
ANALYTIC, RI-free {Phi0, (F-Fbar)Phi0} 2x2 and diagonalize, comparing to the VMC ground truth.

  H = [[E0, h],[h, g]] ,  S = [[1, 0],[0, sigma2]]   (S01=<F-Fbar>=0 exactly),
  E_R12 = 1/2[(E0 + g/sigma2) - sqrt((E0 - g/sigma2)^2 + 4 h^2/sigma2)]   (lowest root).
V_NN adds V_NN*S to H -> shifts every eigenvalue by V_NN (cancels in dE); we assemble the
electronic 2x2 and add V_NN at the end.

ALL matrix elements are analytic RI-free, geminal f = exp(-GAM r) (the Stage-1 f-tensor geminal
that every validated piece uses).  Provenance / reproduce (each validated vs its VMC target):
  E0     = stage2_E0()                        lih_r12ci_energy.py         (-7.887822, control-validated)
  sigma2 = ANALYTIC (block-separable)         lih_r12ci_sigma2_analytic.py ( 0.137228 vs VMC 0.13699)
  h_T    = PartA(grad dens) + PartB(IBP->Yuk) lih_r12ci_hT_analytic.py    (+0.749453 vs +0.7485)
  h_Vne  = <V_ne F> - Fbar<V_ne>              lih_r12ci_sigma2_analytic.py (-0.930200 vs -0.9299)
  h_Vee  = cov_FA_FB(f, coul, Yukawa)         lih_r12ci_hVee_analytic.py  (+0.309790 vs +0.3082)
  g_Vne  = block-decomp <F^2 V_ne> trilinear  lih_r12ci_gVne_analytic.py  (-2.851580 vs -2.8454)
  g_Vee  = <F^2 V_ee> (216-triple enum, incl. the RI-free 3-body TRIANGLE) - 2 Fbar<FV> + Fbar^2<V>
                                              lih_r12ci_gVee_analytic.py  (+0.538342 vs +0.53605)
  g_T    = gT1 - gam^2 sigma^2 + 2 gam Cov[F,Y_sum]   (IBP: no vector nabla-f needed)
                                              lih_r12ci_gT_analytic.py    (+1.394191 vs +1.38640)

NOTE ON GEMINAL.  These pieces (and their VMC targets) use f = exp(-GAM r).  The headline PoC
energy E_R12 = -7.932 (dE = -29.5 mHa) in the build plan is the CUSP-CORRECT geminal f = r e^{-GAM r}
(a different trial function; larger |lowering| here is not "better" -- exp has the wrong cusp sign
but at single-zeta gives a larger variational coupling).  The apples-to-apples validation is the
exp-geminal analytic 2x2 vs the exp-geminal VMC 2x2 (both below).  Extending the analytic
machinery to r e^{-GAM r} is mechanical (the enumerator/triangle/IBP framework is geminal-agnostic;
only build_kernel, the f^2 kernel, the same-pair h_Vee kernel f/r=e^{-gam r}, and grad^2 f change).

Run from debug/:  python lih_r12ci_assemble.py
"""
import numpy as np

GAM = 0.5
V_NN = 0.995025             # Z_A Z_B / R = 3/3.015

# --- analytic RI-free matrix elements (exp geminal), each VMC-validated (see docstring) --- #
E0_tot = -7.887822
sigma2 = 0.137228
h_T, h_Vne, h_Vee = 0.749453, -0.930200, 0.309790
g_T, g_Vne, g_Vee = 1.394191, -2.851580, 0.538342

# --- VMC exp-geminal targets (independent ground truth, same geminal) --- #
VMC = dict(sigma2=0.13699, h=0.1268, g=-0.92300)   # g,sigma2 from _g_targets.py ; h from vmc.py


def two_by_two(E0_elec, h, g_elec, sig2):
    """lowest generalized eigenvalue of [[E0,h],[h,g]] c = E [[1,0],[0,sig2]] c (electronic)."""
    a = E0_elec; c = g_elec / sig2; b = h / np.sqrt(sig2)
    return 0.5 * (a + c - np.sqrt((a - c) ** 2 + 4 * b ** 2))


if __name__ == "__main__":
    print("=" * 78)
    print("LiH R12-CI: analytic RI-free 2x2 -> E_R12  (geminal f = exp(-%.1f r))" % GAM)
    print("=" * 78)

    h = h_T + h_Vne + h_Vee
    g_elec = g_T + g_Vne + g_Vee
    E0_elec = E0_tot - V_NN

    print(f"\n  off-diagonal  h = h_T + h_Vne + h_Vee = {h_T:+.5f} {h_Vne:+.5f} {h_Vee:+.5f} = {h:+.6f}")
    print(f"                  (VMC h_total = +{VMC['h']:.4f})")
    print(f"  diagonal      g = g_T + g_Vne + g_Vee = {g_T:+.5f} {g_Vne:+.5f} {g_Vee:+.5f} = {g_elec:+.6f}")
    print(f"                  (VMC g = {VMC['g']:+.5f})")
    print(f"  overlap       sigma^2 = {sigma2:.6f}   (VMC {VMC['sigma2']:.5f})")
    print(f"  reference     E0 = {E0_tot:.6f} (incl V_NN={V_NN:.5f});  E0_elec = {E0_elec:.6f}")

    # analytic
    ER_elec = two_by_two(E0_elec, h, g_elec, sigma2)
    E_R12 = ER_elec + V_NN
    dE = E_R12 - E0_tot
    print(f"\n  >>> ANALYTIC  E_R12 = {E_R12:.6f} Ha   (dE = E_R12 - E0 = {dE*1e3:+.2f} mHa)")

    # VMC exp-geminal 2x2 (apples-to-apples: same geminal, independent MC pieces)
    ER_vmc = two_by_two(E0_elec, VMC['h'], VMC['g'], VMC['sigma2']) + V_NN
    dE_vmc = ER_vmc - E0_tot
    print(f"      VMC (exp)  E_R12 = {ER_vmc:.6f} Ha   (dE = {dE_vmc*1e3:+.2f} mHa)")
    print(f"      agreement: {abs(E_R12-ER_vmc)*1e3:.2f} mHa")

    print(f"\n  variational checks:  E_R12 < E0 : {E_R12 < E0_tot}   "
          f"E_R12 > exact -8.070 : {E_R12 > -8.070}")
    print(f"  correlation captured: {abs(dE)*1e3:.1f} mHa of LiH's ~83 mHa "
          f"(= {abs(dE)/0.083*100:.0f}%), variational (ionic single-zeta ref, PoC).")
    print(f"\n  (Headline PoC -7.932/-29.5 mHa uses the cusp-correct geminal f=r*exp(-{GAM} r);")
    print(f"   this closes the exp-geminal loop analytically, RI-free, incl. the 3-body triangle.)")
