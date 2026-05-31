"""Bridge Test 3 (v2): the two faces of the wedge partition function meet,
and the bridge between them IS the functional equation.

v1 bug: hand-coded csch Laurent coeffs (wrong factor) evaluated outside radius
-> continuum face garbage. v2 lets mpmath.taylor generate the continuum
coefficients (no hand formula) and evaluates at small beta inside radius pi.

Idealized constant-degeneracy wedge (odd-integer spectrum from Test 1):
    Z(b) = sum_{m odd>0} e^{-b m} = 1/(2 sinh b).
  DISCRETE / state-side face (large b): q-series e^{-b}(1+e^{-2b}+...).
  CONTINUUM / spectral face (small b): 1/(2b) - b/12 + ...; coeffs are Bernoulli,
     tied to zeta(2k) ~ pi^{2k} by the functional equation = THE BRIDGE.
"""
import json
from mpmath import mp, mpf, sinh, exp, pi, zeta, bernoulli, factorial, nsum, inf, taylor
mp.dps = 60

b = mpf("0.5")  # inside radius pi
Z_exact = 1/(2*sinh(b))
Z_discrete = nsum(lambda k: exp(-(2*k+1)*b), [0, inf])

g = lambda x: x/(2*sinh(x)) if x != 0 else mpf("0.5")   # regular part b*Z
N = 20
coeffs = taylor(g, 0, N)
Z_continuum = sum(coeffs[k]*b**(k-1) for k in range(N+1))

bern_check = []
for k in range(0, 6):
    predicted = mpf(1)/2*(2 - mpf(2)**(2*k))*bernoulli(2*k)/factorial(2*k)
    bern_check.append(mp.nstr(abs(coeffs[2*k] - predicted), 3))

fe = []
for k in (1, 2, 3, 4, 5):
    lhs = zeta(2*k)
    rhs = (-1)**(k+1)*bernoulli(2*k)*(2*pi)**(2*k)/(2*factorial(2*k))
    fe.append(mp.nstr(abs(lhs-rhs), 3))

out = {"beta": str(b), "Z_exact": mp.nstr(Z_exact, 25),
       "discrete_qseries_face": mp.nstr(Z_discrete, 25),
       "continuum_bernoulli_face": mp.nstr(Z_continuum, 25),
       "discrete_resid": mp.nstr(abs(Z_discrete-Z_exact), 4),
       "continuum_resid": mp.nstr(abs(Z_continuum-Z_exact), 4),
       "continuum_coeffs_are_Bernoulli_resid_k0to5": bern_check,
       "functional_eqn_zeta2k_eq_pi2k_Bernoulli_resid_k1to5": fe}
json.dump(out, open("debug/data/bridge_test3_v2.json", "w"), indent=1)
print("Z_exact              :", out["Z_exact"])
print("discrete (q-series)  : resid", out["discrete_resid"])
print("continuum (Bernoulli): resid", out["continuum_resid"])
print("continuum coeffs = Bernoulli? resid k=0..5:", bern_check)
print("functional eqn resid k=1..5:", fe)
