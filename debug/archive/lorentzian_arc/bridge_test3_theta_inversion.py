"""Bridge Test 3: is the discrete<->continuum bridge the functional equation?

The wedge KMS partition function (state-side, large-beta, q-series face) and
the spectral-action heat trace (spectral-side, small-beta, Mellin/Bernoulli
face) are two analytic faces of one spectral function built from the SAME
odd-integer spectrum found in Test 1. Test whether the bridge between the faces
is the theta / Bernoulli functional equation zeta(2k) = pi^{2k} * Bernoulli
(= GB-4 two-window duality = Bernoulli-ladder zeta(s)<->zeta(1-s)).

Mechanism demonstrator: constant-degeneracy wedge Z(b) = sum_{m odd>0} e^{-b m}
= 1/(2 sinh b).
   large-b face: q-series e^{-b}(1+e^{-2b}+...)            [DISCRETE / state-side]
   small-b face: 1/(2b) - b/12 + 7 b^3/720 - ...           [CONTINUUM / spectral]
The small-b coefficients are Bernoulli B_2k/(2k)! -> tie to zeta(2k) ~ pi^{2k}
by the functional equation. That tie IS the bridge.

HONEST SCOPE: this confirms the MECHANISM on the idealized wedge and verifies
the zeta(2k)=pi^{2k}*Bernoulli identity exactly. The GeoVac-specific version
(real growing degeneracies) is the Eisenstein-class generalization and rests on
the already-verified Sprint MR-B Jacobi theta inversion of the Dirac heat trace.
"""
import json
from mpmath import mp, mpf, sinh, exp, pi, zeta, bernoulli, factorial, nsum, inf
mp.dps = 50

beta = mpf("0.37")
Z_exact = 1/(2*sinh(beta))
Z_qseries = nsum(lambda k: exp(-(2*k+1)*beta), [0, inf])           # discrete face

def csch_over2_series(b, N=16):                                     # continuum face
    return sum(2*(mpf(2)**(1-2*n)-1)*bernoulli(2*n)/factorial(2*n)*b**(2*n-1)
               for n in range(N)) / 2
Z_bern = csch_over2_series(beta)

fe = []                                                            # the bridge content
for k in (1, 2, 3, 4, 5):
    lhs = zeta(2*k)
    rhs = (-1)**(k+1)*bernoulli(2*k)*(2*pi)**(2*k)/(2*factorial(2*k))
    fe.append(mp.nstr(abs(lhs-rhs), 3))

out = {"beta": str(beta), "Z_exact": str(Z_exact),
       "discrete_qseries_face_resid": mp.nstr(abs(Z_qseries-Z_exact), 4),
       "continuum_bernoulli_face_resid": mp.nstr(abs(Z_bern-Z_exact), 4),
       "zeta2k_eq_pi2k_bernoulli_resid_k1to5": fe,
       "leading_continuum_term": "1/(2 beta)  (Weyl/Mellin pole)",
       "first_correction": "-beta/12 = -zeta(2)/(4 pi^2) * beta  (M2-ring via FE)"}
json.dump(out, open("debug/data/bridge_test3.json", "w"), indent=1)
print("Z_exact                  :", mp.nstr(Z_exact, 16))
print("discrete (q-series) face : resid", out["discrete_qseries_face_resid"])
print("continuum (Bernoulli)face: resid", out["continuum_bernoulli_face_resid"])
print("zeta(2k)=pi^2k*Bernoulli : resid k=1..5", fe)
