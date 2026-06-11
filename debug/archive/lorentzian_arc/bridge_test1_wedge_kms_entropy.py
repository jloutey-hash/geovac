"""Bridge Test 1: is the wedge KMS entropy a modular (odd-power) q-series?

Pull the TRUE wedge spectrum of K_alpha = J_polar, build rho_W = e^{-beta K_W}/Z
at the BW point beta = 2*pi, compute S(rho_W) = -Tr(rho_W log rho_W), and check
whether S is an odd-power q-series in q = e^{-2*pi}.

PASS  = clean odd-power q-series (modular / theta arithmetic).
FAIL  = carries pi-powers or non-modular content.
"""
import math
import numpy as np

import geovac.modular_hamiltonian as M

BETA = 2.0 * math.pi
Q = math.exp(-BETA)  # q = e^{-2pi} ~ 1.867e-3

def wedge_K_eigs(n_max):
    """Return the wedge-restricted K_alpha eigenvalues (the integers it keeps)."""
    mh = M.for_bisognano_wichmann(n_max)
    # Preferred API
    for name in ("restrict_K_alpha_to_wedge", "restrict_K_to_wedge"):
        fn = getattr(mh, name, None)
        if fn is not None:
            Kw = np.asarray(fn())
            eigs = np.linalg.eigvalsh(Kw) if Kw.ndim == 2 else np.asarray(Kw)
            return np.real(eigs), name
    raise RuntimeError("no wedge-restriction method found: " +
                       ",".join(a for a in dir(mh) if "wedge" in a.lower()))

def entropy_from_eigs(eigs, beta=BETA):
    w = np.exp(-beta * np.asarray(eigs, float))
    Z = w.sum()
    p = w / Z
    p = p[p > 0]
    S = float(-(p * np.log(p)).sum())
    return S, float(Z)

print(f"q = e^(-2pi) = {Q:.6e}")
print(f"{'n_max':>5} {'wedge_dim':>9} {'distinct_eigs':>28} {'Z':>12} {'S(rho_W)':>16}")
results = {}
for n_max in (2, 3, 4, 5):
    try:
        eigs, api = wedge_K_eigs(n_max)
    except Exception as e:
        print(f"{n_max:>5}  ERROR: {e}")
        continue
    # group eigenvalues -> (value, degeneracy)
    vals, counts = np.unique(np.round(eigs, 9), return_counts=True)
    S, Z = entropy_from_eigs(eigs)
    results[n_max] = (vals, counts, S, Z)
    disp = ", ".join(f"{v:g}x{c}" for v, c in zip(vals, counts))
    print(f"{n_max:>5} {len(eigs):>9} {disp:>28} {Z:>12.8f} {S:>16.12f}")

# --- q-series / modular structure check -------------------------------------
# If eigenvalues are positive integers m with degeneracies g_m, then
#   Z = sum_m g_m q^m,   S = beta*<K> + log Z,   <K> = sum_m m g_m q^m / Z.
# Verify S reconstructs from a pure q-power series (q = e^{-2pi}); report the
# leading exponent (odd vs even) and how fast it saturates.
print("\n--- q-series reconstruction (q = e^{-2pi}) ---")
for n_max, (vals, counts, S, Z) in results.items():
    ints = np.allclose(vals, np.round(vals))
    all_pos = np.all(vals > 0)
    leading = vals[0] if len(vals) else None
    parity = ("odd" if abs(leading - round(leading)) < 1e-9 and int(round(leading)) % 2 == 1
              else "even/other")
    # closed-form q-series rebuild
    g = {int(round(v)): int(c) for v, c in zip(vals, counts)} if ints else {}
    Zq = sum(c * Q**m for m, c in g.items()) if g else float("nan")
    Kq = sum(m * c * Q**m for m, c in g.items()) / Zq if g else float("nan")
    Sq = BETA * Kq + math.log(Zq) if g else float("nan")
    print(f"n_max={n_max}: integer_spectrum={ints}  all_positive={all_pos}  "
          f"leading_exp={leading:g} ({parity})  "
          f"S_qseries={Sq:.12f}  |S-S_qseries|={abs(S-Sq):.2e}")

print("\nLeading-order S ~ g1*q^(m1)*(beta*m1 - log(g1*q^m1)) ...")
print("(small q => S is exponentially small and dominated by the lowest gap)")
