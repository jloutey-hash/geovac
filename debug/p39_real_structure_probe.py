"""Paper 39 real-structure probe: KO-6 real structure J_ab for the
Clifford-doubled product of two Camporesi-Higuchi S^3 (KO-3) triples.

Construction under test:
    H_ab  = H_a (x) C^2 (x) H_b
    D_ab  = D_a (x) s1 (x) I_b  +  I_a (x) s2 (x) D_b
    g_ab  = I_a (x) s3 (x) I_b                       (grading)
    A_ab  = C^inf(S^3) (x) I_2 (x) C^inf(S^3)        (a (x) I_2 (x) b)

Candidate real structure:
    J_ab  = J_a (x) (W . K_2) (x) J_b ,   W in {I2, s1, s2, s3}
antilinear;  as a matrix  J_ab(psi) = U_ab @ conj(psi)  with
    U_ab = U_a (x) W (x) U_b .

We test all four W to show which (if any) satisfies the KO-6 signs
    (eps, eps', eps'') = (+, +, -)   [van Suijlekom Table 3.1]
and that the choice is FORCED, not a 2x2 coincidence -- every sign is
checked against the ACTUAL D_a, D_b, J_a, J_b matrices (which carry the
m_j-flip permutation and phases), per the "suspiciously clean toy" warning.

No paper / module / test is modified.
"""
from __future__ import annotations

import numpy as np

from geovac.full_dirac_operator_system import (
    FullDiracTruncatedOperatorSystem,
    camporesi_higuchi_full_dirac_matrix,
)
from geovac.real_structure import build_J_full_dirac

I2 = np.eye(2, dtype=np.complex128)
S1 = np.array([[0, 1], [1, 0]], dtype=np.complex128)
S2 = np.array([[0, -1j], [1j, 0]], dtype=np.complex128)
S3 = np.array([[1, 0], [0, -1]], dtype=np.complex128)
PAULI = {"I2": I2, "s1": S1, "s2": S2, "s3": S3}

kron = np.kron


def build_factor(n_max: int):
    op = FullDiracTruncatedOperatorSystem(n_max)
    D = camporesi_higuchi_full_dirac_matrix(op.basis)
    J = build_J_full_dirac(n_max)          # J.U is the unitary part; J = U K
    return op, D, J.U


def anti_apply_op(U_ab: np.ndarray, op: np.ndarray) -> np.ndarray:
    """J op J^{-1} for antilinear J = U K:  U conj(op) U^T."""
    return U_ab @ np.conj(op) @ U_ab.T


def check_product(n_a: int, n_b: int, W_name: str, tol: float = 1e-10,
                  do_order: bool = True):
    op_a, D_a, U_a = build_factor(n_a)
    op_b, D_b, U_b = build_factor(n_b)
    dA, dB = D_a.shape[0], D_b.shape[0]
    IA, IB = np.eye(dA, dtype=np.complex128), np.eye(dB, dtype=np.complex128)
    W = PAULI[W_name]

    # product operators
    D_ab = kron(kron(D_a, S1), IB) + kron(kron(IA, S2), D_b)
    g_ab = kron(kron(IA, S3), IB)
    U_ab = kron(kron(U_a, W), U_b)                       # unitary part of J_ab

    dim = D_ab.shape[0]
    I = np.eye(dim, dtype=np.complex128)

    # --- structural sanity on D_ab, g_ab (independent of J) -------------
    D_sq = D_ab @ D_ab
    D_sq_expect = kron(kron(D_a @ D_a, I2), IB) + kron(kron(IA, I2), D_b @ D_b)
    res_Dsq = float(np.max(np.abs(D_sq - D_sq_expect)))
    # two summands anticommute
    A = kron(kron(D_a, S1), IB)
    B = kron(kron(IA, S2), D_b)
    res_anti = float(np.max(np.abs(A @ B + B @ A)))
    # grading: g^2=I, g*=g, {g,D}=0, [g,A]=0 (A commutes trivially, middle I2)
    res_g2 = float(np.max(np.abs(g_ab @ g_ab - I)))
    res_gD = float(np.max(np.abs(g_ab @ D_ab + D_ab @ g_ab)))

    # --- the three KO sign relations (antilinear J) --------------------
    # J^2 = U conj(U)
    J2 = U_ab @ np.conj(U_ab)
    eps = _classify_sign(J2, I, tol)
    res_eps = _sign_residual(J2, I, eps)

    # J D = eps' D J   <=>   U conj(D) = eps' D U
    lhs = U_ab @ np.conj(D_ab)
    epsp = _classify_sign(lhs, D_ab @ U_ab, tol)
    res_epsp = _sign_residual(lhs, D_ab @ U_ab, epsp)

    # J g = eps'' g J  <=>  U conj(g) = eps'' g U ; g real => U g = eps'' g U
    lhsg = U_ab @ np.conj(g_ab)
    epspp = _classify_sign(lhsg, g_ab @ U_ab, tol)
    res_epspp = _sign_residual(lhsg, g_ab @ U_ab, epspp)

    # unitarity of U_ab
    res_unit = float(np.max(np.abs(U_ab @ U_ab.conj().T - I)))

    # --- order-zero / order-one on A_ab -------------------------------
    # generators: a (x) I2 (x) I_b  and  I_a (x) I2 (x) b
    gens = []
    for M in op_a.multiplier_matrices:
        gens.append(kron(kron(M, I2), IB))
    for M in op_b.multiplier_matrices:
        gens.append(kron(kron(IA, I2), M))

    o0_fail, o0_max = 0, 0.0
    o1_fail, o1_max = 0, 0.0
    if not do_order:
        return dict(
            n_a=n_a, n_b=n_b, W=W_name, dim=dim, n_gens=len(gens),
            res_Dsq=res_Dsq, res_anti=res_anti, res_g2=res_g2, res_gD=res_gD,
            res_unit=res_unit,
            eps=eps, res_eps=res_eps, epsp=epsp, res_epsp=res_epsp,
            epspp=epspp, res_epspp=res_epspp,
            o0_fail=-1, o0_max=-1.0, n_pairs=len(gens) ** 2,
            o1_fail=-1, o1_max=-1.0,
        )
    Da_comm = [D_ab @ a - a @ D_ab for a in gens]
    JbJ_list = [anti_apply_op(U_ab, b) for b in gens]   # precompute once
    for a in gens:
        for JbJ in JbJ_list:
            c = a @ JbJ - JbJ @ a
            r = float(np.max(np.abs(c)))
            if r > tol:
                o0_fail += 1
                o0_max = max(o0_max, r)
    for Da in Da_comm:
        for JbJ in JbJ_list:
            c = Da @ JbJ - JbJ @ Da
            r = float(np.max(np.abs(c)))
            if r > tol:
                o1_fail += 1
                o1_max = max(o1_max, r)

    return dict(
        n_a=n_a, n_b=n_b, W=W_name, dim=dim, n_gens=len(gens),
        res_Dsq=res_Dsq, res_anti=res_anti, res_g2=res_g2, res_gD=res_gD,
        res_unit=res_unit,
        eps=eps, res_eps=res_eps,
        epsp=epsp, res_epsp=res_epsp,
        epspp=epspp, res_epspp=res_epspp,
        o0_fail=o0_fail, o0_max=o0_max, n_pairs=len(gens) ** 2,
        o1_fail=o1_fail, o1_max=o1_max,
    )


def _classify_sign(A, B, tol):
    """Return +1 if A==B, -1 if A==-B, 0 otherwise (to tol)."""
    if float(np.max(np.abs(A - B))) < tol:
        return +1
    if float(np.max(np.abs(A + B))) < tol:
        return -1
    return 0


def _sign_residual(A, B, sign):
    if sign == +1:
        return float(np.max(np.abs(A - B)))
    if sign == -1:
        return float(np.max(np.abs(A + B)))
    # neither: report distance to nearest
    return min(float(np.max(np.abs(A - B))), float(np.max(np.abs(A + B))))


if __name__ == "__main__":
    print("=" * 78)
    print("SINGLE-FACTOR order-0/order-1 finite-cutoff trend (KO-3, S^3)")
    print("=" * 78)
    from geovac.real_structure import audit_J
    for nm in (1, 2, 3):
        r = audit_J(nm, sector="full_dirac", D_mode="truthful")
        print(f"  n_max={nm}: dim_O={r['dim_O']:3d}  J^2=-I {r['J_squared_minus_I'][1]:.1e}"
              f"  JD=+DJ {r['J_commutes_D'][1]:.1e}"
              f"  order0 fails={r['order_zero'][1]:3d} max={r['order_zero'][2]:.3e}"
              f"  order1 fails={r['order_one'][1]:3d} max={r['order_one'][2]:.3e}")

    print()
    print("=" * 78)
    print("PRODUCT KO-6:  all four middle operators W, at (n_a,n_b)=(2,2)")
    print("  target KO-6 signs (eps, eps', eps'') = (+, +, -)")
    print("=" * 78)
    header = ("  W     dim  eps  res_eps   eps'  res_eps'  eps'' res_eps''"
              "  {D_a,D_b}  D_ab^2   {g,D}")
    print(header)
    for W in ("I2", "s1", "s2", "s3"):
        r = check_product(2, 2, W, do_order=False)
        print(f"  {W:4s} {r['dim']:4d}  {r['eps']:+d}   {r['res_eps']:.1e}"
              f"   {r['epsp']:+d}   {r['res_epsp']:.1e}"
              f"   {r['epspp']:+d}   {r['res_epspp']:.1e}"
              f"  {r['res_anti']:.1e}  {r['res_Dsq']:.1e}  {r['res_gD']:.1e}")

    print()
    print("=" * 78)
    print("PRODUCT with W=s1 (the forced choice):  full axiom set")
    print("=" * 78)
    for (na, nb) in [(1, 1), (1, 2), (2, 1), (2, 2)]:
        r = check_product(na, nb, "s1")
        print(f"  (n_a,n_b)=({na},{nb}) dim={r['dim']:3d} gens={r['n_gens']:2d}"
              f" | U unit {r['res_unit']:.1e}"
              f" | J^2=+I ({r['eps']:+d},{r['res_eps']:.1e})"
              f" | JD=+DJ ({r['epsp']:+d},{r['res_epsp']:.1e})"
              f" | Jg=-gJ ({r['epspp']:+d},{r['res_epspp']:.1e})")
        print(f"        order0: {r['o0_fail']}/{r['n_pairs']} fail max={r['o0_max']:.3e}"
              f"   order1: {r['o1_fail']}/{r['n_pairs']} fail max={r['o1_max']:.3e}")
