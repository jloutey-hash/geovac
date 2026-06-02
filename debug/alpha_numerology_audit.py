"""
alpha_numerology_audit.py
Reality-check on the GeoVac alpha observation K = pi*(B + F - Delta) ~ 1/alpha.

Question: is the numerical match to 1/alpha a genuinely sparse coincidence,
or is it dense (i.e. comparably-simple expressions hit 137.036 easily)?

This tests ONLY the numerical-surprise component. The framework's structural
stories for B=42 (Casimir), F=zeta(2) (Fock-Dirichlet), Delta=1/40 (Dirac
degeneracy) are post-hoc and labeled conjectural by the project itself; this
script does not adjudicate those stories. It asks the narrower, decidable
question: how surprising is the NUMBER match, on its own?
"""
import math
from math import gcd

try:
    from mpmath import mp, mpf, pi as MPI, zeta as MZETA
    mp.dps = 50
    HIGH = True
except Exception:
    HIGH = False

INV_ALPHA = 137.035999177  # CODATA ~2018 value of 1/alpha
PI = math.pi
ZETA2 = PI * PI / 6.0

# ---- 1. Verify the claimed expression ----
K_claim = PI * (42 + ZETA2 - 1.0 / 40.0)
rel_claim = abs(K_claim - INV_ALPHA) / INV_ALPHA
print("=" * 66)
print("1. The claimed expression")
print(f"   1/alpha                = {INV_ALPHA:.9f}")
print(f"   pi*(42 + zeta2 - 1/40) = {K_claim:.9f}")
print(f"   absolute diff          = {K_claim - INV_ALPHA:+.3e}")
print(f"   relative diff          = {rel_claim:.3e}")
print(f"   alpha^2                = {(1/INV_ALPHA)**2:.3e}   (is the diff ~ alpha^2?)")
if HIGH:
    Kh = MPI * (42 + MZETA(2) - mpf(1) / 40)
    ia = mpf('137.035999177')
    print(f"   [mpmath 50dps] relative diff = {abs(Kh - ia) / ia}")
print()

# ---- 2. Build comparable-simplicity menus ----
const_menu = {
    "0": 0.0, "1": 1.0, "3/2": 1.5,
    "zeta2=pi^2/6": ZETA2, "zeta3": 1.2020569031595943,
    "pi^2/12": PI * PI / 12, "pi/2": PI / 2, "pi": PI,
    "sqrt2": math.sqrt(2), "sqrt3": math.sqrt(3), "sqrt5": math.sqrt(5),
    "ln2": math.log(2), "ln3": math.log(3),
    "e": math.e, "e/2": math.e / 2,
    "gamma": 0.5772156649015329, "Catalan": 0.915965594177219,
    "phi": (1 + math.sqrt(5)) / 2,
}
delta_menu = {"0": 0.0}
for n in range(1, 51):
    delta_menu[f"1/{n}"] = 1.0 / n
for m in (1, 2, 3):
    for n in range(2, 13):
        if gcd(m, n) == 1:
            delta_menu[f"{m}/{n}"] = m / n
B_range = range(0, 61)
total = len(B_range) * len(const_menu) * len(delta_menu)

# ---- 3. Count matches at / below the claim tolerance ----
matches = []
for B in B_range:
    for fn, F in const_menu.items():
        for dn, D in delta_menu.items():
            K = PI * (B + F - D)
            rel = abs(K - INV_ALPHA) / INV_ALPHA
            if rel <= rel_claim:
                is_gv = (B == 42 and fn.startswith("zeta2") and dn == "1/40")
                matches.append((rel, B, fn, dn, K, is_gv))
matches.sort()
print("=" * 66)
print("2. Comparable-simplicity search:  K = pi*(B + F - Delta)")
print(f"   menu: B in 0..60, F in {len(const_menu)} constants, Delta in {len(delta_menu)} small rationals")
print(f"   total expressions tried = {total}")
print(f"   # matching 1/alpha AT LEAST AS WELL as GeoVac (rel <= {rel_claim:.2e}) = {len(matches)}")
print("   best 15:")
for rel, B, fn, dn, K, is_gv in matches[:15]:
    print(f"     rel={rel:.2e}  pi*({B:>2} + {fn:<12} - {dn:<5}) = {K:.9f}{'   <== GeoVac' if is_gv else ''}")
print()

# ---- 4. Density curve ----
print("=" * 66)
print("3. Density vs tolerance (how the count grows as precision is relaxed)")
for tol in (1e-3, 1e-4, 1e-5, 1e-6, 1e-7):
    c = 0
    for B in B_range:
        for F in const_menu.values():
            for D in delta_menu.values():
                if abs(PI * (B + F - D) - INV_ALPHA) / INV_ALPHA <= tol:
                    c += 1
    print(f"   rel <= {tol:.0e} : {c:>5} expressions")
print("=" * 66)
