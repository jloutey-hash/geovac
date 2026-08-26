"""Paper 56 rem:kms_hodge_circle -- the thermal-time (KMS/BW) circle realises the
Hodge circle of prop:hodge_cm_point on the j=1/2 doublet.  All checks EXACT (sympy).

  K = diag(2 m_j) = diag(1,-1)   (compact BW modular generator, KMS beta = 2 pi)
  J = [[0,-1],[1,0]]             (Kramers/Hodge complex structure, prop:hodge_cm_point)
  U = Kramers frame              (columns = J eigenvectors for +i, -i)

Claims: (1) e^{2 pi i K} = I (closure); (2) e^{i pi K} = -I on the spinor doublet vs
+I on integer-l scalars (spin double cover at beta/2); (3) J = the flow's beta/4
quarter-period point; the quarter-points are mu_4 = {I, J, -I, -J}; (4) the whole flow
transports to the Hodge circle cos t I + sin t J; (5) Q[J] = Q(i), norm-1 torus.
Promoted from debug/_kms_mt_torus.py (12/12, 2026-08-21)."""
import sympy as sp

t = sp.symbols('t', real=True)
x = sp.symbols('x')
I2 = sp.eye(2)
Z2 = sp.zeros(2)
K = sp.diag(1, -1)
J = sp.Matrix([[0, -1], [1, 0]])
Q = sp.Matrix([[0, 1], [-1, 0]])
U = sp.Matrix([[1, 1], [-sp.I, sp.I]]) / sp.sqrt(2)


def test_kms_2pi_closure():
    assert sp.simplify((2 * sp.pi * sp.I * K).exp() - I2) == Z2


def test_spin_double_cover_at_half_period():
    # spinor doublet: -I at beta/2; integer-l scalar sector: +I
    assert sp.simplify((sp.pi * sp.I * K).exp() + I2) == Z2
    K_scalar = sp.diag(2, 0, -2)
    assert sp.simplify((sp.pi * sp.I * K_scalar).exp() - sp.eye(3)) == sp.zeros(3)


def test_J_is_quarter_period_and_mu4():
    assert sp.simplify(U * U.H - I2) == Z2
    assert sp.simplify(U * (sp.pi * sp.I * K / 2).exp() * U.H - J) == Z2
    quarters = [sp.simplify(U * (sp.I * sp.Rational(k, 2) * sp.pi * K).exp() * U.H)
                for k in range(4)]
    assert quarters == [I2, J, -I2, -J]


def test_flow_is_hodge_circle():
    flowK = (sp.I * t * K).exp()
    hodge = sp.cos(t) * I2 + sp.sin(t) * J
    assert sp.simplify(sp.expand_complex(U * flowK * U.H - hodge)) == Z2


def test_QJ_is_Qi_norm_one_torus():
    assert J * J == -I2
    assert sp.expand(J.charpoly(x).as_expr() - (x**2 + 1)) == 0
    assert Q * J == I2                                  # Riemann form (polarization)
    a, b, c, d = sp.symbols('a b c d', real=True)
    prod = (a * I2 + b * J) * (c * I2 + d * J)
    assert sp.simplify(prod - ((a * c - b * d) * I2 + (a * d + b * c) * J)) == Z2
    assert sp.simplify((a * I2 + b * J).det() - (a**2 + b**2)) == 0
