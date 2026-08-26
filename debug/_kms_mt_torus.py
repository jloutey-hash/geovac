"""T1: KMS/BW modular circle vs Kramers/Hodge Mumford-Tate torus on the j=1/2 doublet.
All checks EXACT (sympy). Objects taken verbatim from the corpus:
  K = diag(2*m_j) = diag(1,-1)  -- BW modular generator (WH7 Lorentzian closure)
  J = [[0,-1],[1,0]]            -- Kramers/Hodge complex structure (P56 prop:hodge_cm_point)
  Q = [[0,1],[-1,0]]            -- SL2-invariant symplectic form (P56)
"""
import sympy as sp

t = sp.symbols('t', real=True)
x = sp.symbols('x')
I2 = sp.eye(2); Z2 = sp.zeros(2)

K = sp.diag(1, -1)
J = sp.Matrix([[0, -1], [1, 0]])
Q = sp.Matrix([[0, 1], [-1, 0]])

ok = []

# (1) KMS beta=2pi closure: e^{2 pi i K} = I   [WH7: compactness => discreteness]
ok.append(("e^{2 pi i K} = I (KMS closure)",
           sp.simplify((2*sp.pi*sp.I*K).exp() - I2) == Z2))

# (2) HALF-period: e^{i pi K} = -I  -- the spin double cover INSIDE the thermal circle.
#     (Scalar sector: 2*m_l even => e^{i pi K_scalar} = +I. The spinor/scalar discriminant
#      lives HERE, in the flow, not in the ERI.)
ok.append(("e^{i pi K} = -I (spin double cover at beta/2)",
           sp.simplify((sp.pi*sp.I*K).exp() + I2) == Z2))
K_scalar = sp.diag(2, 0, -2)   # 2*m_l for l=1, integer l
ok.append(("scalar sector: e^{i pi K_l} = +I (no double cover)",
           sp.simplify((sp.pi*sp.I*K_scalar).exp() - sp.eye(3)) == sp.zeros(3)))

# (3) QUARTER-period: e^{i (pi/2) K} = diag(i,-i), which transported to the Kramers
#     frame IS J. So J = the beta/4 point of the modular flow.
U = sp.Matrix([[1, 1], [-sp.I, sp.I]]) / sp.sqrt(2)   # Kramers frame (cols = J eigvecs for +i,-i)
ok.append(("U unitary", sp.simplify(U*U.H - I2) == Z2))
ok.append(("U diag(i,-i) U^dag = J  (J = flow at beta/4)",
           sp.simplify(U*(sp.pi*sp.I*K/2).exp()*U.H - J) == Z2))

# (4) The WHOLE flow transports to the Hodge circle: U e^{itK} U^dag = cos t I + sin t J
flowK = (sp.I*t*K).exp()
hodge = sp.cos(t)*I2 + sp.sin(t)*J
ok.append(("U e^{itK} U^dag = cos t I + sin t J (flow = Hodge circle)",
           sp.simplify(sp.expand_complex(U*flowK*U.H - hodge)) == Z2))

# (5) Q[J] = Q(i): minimal polynomial x^2+1, irreducible over Q; J^2 = -I
ok.append(("charpoly(J) = x^2+1", sp.expand(J.charpoly(x).as_expr() - (x**2+1)) == 0))
ok.append(("J^2 = -I", J*J == -I2))

# (6) Norm-1 torus of Q(i): (aI+bJ)(cI+dJ) = complex multiplication; det = Norm
a,b,c,d = sp.symbols('a b c d', real=True)
ok.append(("aI+bJ multiplies as Q(i)",
           sp.simplify((a*I2+b*J)*(c*I2+d*J) - ((a*c-b*d)*I2+(a*d+b*c)*J)) == Z2))
ok.append(("det(aI+bJ) = a^2+b^2 = N_{Q(i)/Q}",
           sp.simplify((a*I2+b*J).det() - (a**2+b**2)) == 0))

# (7) Riemann form QJ = I positive definite (polarization, P56)
ok.append(("QJ = I", Q*J == I2))

# (8) mu_4 torsion: the four quarter-points of the flow in the Kramers frame are {I,J,-I,-J}
quarters = [sp.simplify(U*(sp.I*sp.Rational(k,2)*sp.pi*K).exp()*U.H) for k in range(4)]
ok.append(("quarter points = {I, J, -I, -J} = mu_4",
           quarters == [I2, J, -I2, -J]))

for name, passed in ok:
    print(("PASS  " if passed else "FAIL  ") + name)
print("\nALL PASS" if all(p for _,p in ok) else "\nSOME FAILED")
