"""
Sprint GB: the Bernoulli ladder and the functional-equation reading of the
two-layer split. Diagnostic-only. Exact sympy (skeleton quantities).

Question: are the gravity-sector rational coefficients and the alpha-conjecture
F = pi^2/6 shadows of a single Bernoulli ladder B_{2n}, glued by the zeta
functional equation -- i.e. is the Layer-1/Layer-2 split the ratio
zeta(1-s) <-> zeta(s)?
"""
import sympy as sp

pi = sp.pi

def show(label, val):
    print(f"  {label:<46} = {val}")

print("=" * 78)
print("1. Bernoulli numbers and the two zeta windows")
print("=" * 78)
for n in (1, 2, 3):                       # B_2, B_4, B_6
    s_pos = 2 * n                          # zeta(2), zeta(4), zeta(6)
    s_neg = 1 - 2 * n                      # zeta(-1), zeta(-3), zeta(-5)
    B = sp.bernoulli(2 * n)
    z_pos = sp.zeta(s_pos)
    z_neg = sp.zeta(s_neg)
    # standard closed forms
    z_pos_cf = (-1)**(n + 1) * B * (2 * pi)**(2 * n) / (2 * sp.factorial(2 * n))
    z_neg_cf = -B / (2 * n)
    print(f"\n B_{2*n} = {B}")
    show(f"zeta({s_pos})  [M2 / transcendental window]", z_pos)
    show(f"   closed form  (-1)^(n+1) B (2pi)^2n / 2(2n)!", sp.simplify(z_pos_cf))
    show(f"   matches?", sp.simplify(z_pos - z_pos_cf) == 0)
    show(f"zeta({s_neg}) [skeleton / rational window]", z_neg)
    show(f"   closed form  -B/(2n)", z_neg_cf)
    show(f"   matches?", sp.simplify(z_neg - z_neg_cf) == 0)

print("\n" + "=" * 78)
print("2. GeoVac gravity-sector coefficients, expressed in zeta/Bernoulli")
print("=" * 78)

z_m1 = sp.zeta(-1)     # -1/12
z_m3 = sp.zeta(-3)     #  1/120
z_2  = sp.zeta(2)      # pi^2/6
z_4  = sp.zeta(4)      # pi^4/90

coeffs = {
    "conical tip  (1/12)(1/a - a), coeff 1/12 (G4-2)":  sp.Rational(1, 12),
    "replica entropy derivative  1/6 (G4-4f)":          sp.Rational(1, 6),
    "per-t UV target  1/(24 pi t), coeff 1/24 (v3.20)": sp.Rational(1, 24),
    "scalar Casimir S^3  1/240 (Paper 35)":             sp.Rational(1, 240),
    "Dirac Casimir S^3  17/480 (Paper 35)":             sp.Rational(17, 480),
    "alpha-conj F = pi^2/6 = zeta(2) (Paper 2)":         z_2,
    "kappa = -1/16 (Paper 0, control: should NOT fit)": sp.Rational(-1, 16),
}

# candidate "skeleton-window" anchors
print("\n  anchors:")
show("zeta(-1)", z_m1); show("zeta(-3)", z_m3)
show("zeta(2)", z_2);   show("zeta(4)", z_4)
show("B_2", sp.bernoulli(2)); show("B_4", sp.bernoulli(4))

print("\n  coefficient -> simplest zeta/Bernoulli expression (exact):")
for name, c in coeffs.items():
    rels = []
    if c == -z_m1:            rels.append("= -zeta(-1) = B_2/2")
    if c == -2 * z_m1:        rels.append("= -2 zeta(-1) = B_2")
    if c == sp.Rational(1,2)*(-z_m1)/sp.Rational(1,2):  pass
    if c == -z_m1 / 2:        rels.append("= -zeta(-1)/2 = B_2/4")
    if c == -z_m3 / 2:        rels.append("= -zeta(-3)/2 = B_4-shadow")
    if c == sp.Rational(17,480): rels.append("= 17/480 (NOT a single -zeta(-3) multiple; see Hurwitz check)")
    if c == z_2:              rels.append("= zeta(2) = pi^2 B_2  [M2 window of B_2]")
    if c == sp.Rational(-1,16): rels.append("NO clean zeta/Bernoulli form (control)")
    show(name, " ; ".join(rels) if rels else "(no match found)")

print("\n" + "=" * 78)
print("3. Internal-consistency check: tip -> replica derivative -> entropy")
print("=" * 78)
a = sp.symbols('alpha', positive=True)
tip = sp.Rational(1, 12) * (1/a - a)
dtip = sp.diff(tip, a)
show("tip coefficient a_tip(a) = (1/12)(1/a - a)", tip)
show("d/da a_tip at a=1  (replica derivative)", dtip.subs(a, 1))
show("  equals -1/6 = -B_2 ?", dtip.subs(a, 1) == sp.Rational(-1, 6))
show("  equals 2*zeta(-1) ?", dtip.subs(a, 1) == 2 * z_m1)

print("\n" + "=" * 78)
print("4. Dirac Casimir 17/480: half-integer Hurwitz (M3-flavored) check")
print("=" * 78)
# scalar S^3 Casimir uses zeta(-3); Dirac uses half-integer-shifted Hurwitz.
# Hurwitz zeta(-3, 1/2) carries the half-integer shift -> different rational.
hz = sp.zeta(-3, sp.Rational(1, 2))
show("Hurwitz zeta(-3, 1/2)", hz)
show("scalar S^3 Casimir 1/240 vs -zeta(-3)/2", sp.Rational(1,240) - (-z_m3/2))
show("  (zero => scalar Casimir IS -zeta(-3)/2)", sp.Rational(1,240) == -z_m3/2)
# is 17/480 expressible via integer + half-integer zeta(-3)?
show("17/480 - (-zeta(-3)/2)", sp.simplify(sp.Rational(17,480) - (-z_m3/2)))
show("  -> residual is rational; 17 needs the half-integer shift (M3), not pure M2")
