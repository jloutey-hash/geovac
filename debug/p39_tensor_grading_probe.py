"""Grading / Clifford structure checks for Paper 39's tensor Dirac.

G1  Does a chirality grading gamma_a exist on L^2(S^3, Sigma)?
    Requirement (forced by boundedness of [D_ab, a (x) b]): gamma_a commutes
    with C^inf(S^3) => it is a bundle endomorphism => at a point it must
    anticommute with the whole Clifford action c(T*_x M) = span{sigma_1,2,3}.
    Solve the linear system  {gamma, sigma_i} = 0, i = 1,2,3, over M_2(C).

G2  Clifford-doubled odd (x) odd product (the correct construction):
        H = H_a (x) C^2 (x) H_b,
        D = D_a (x) sigma_1 (x) 1  +  1 (x) sigma_2 (x) D_b,
        grading = 1 (x) sigma_3 (x) 1.
    Symbol level: A = c_a(v) (x) s1 (x) 1,  B = 1 (x) s2 (x) c_b(w).
    Check {A,B} = 0, (A+B)^2 = (|v|^2+|w|^2) I, ||A+B|| = sqrt(|v|^2+|w|^2).

G3  Paper 39's eq:anticomm with gamma_a := sign(D_a) (Paper 38's corrected
    reading of the diag(+1,-1) label): does {D_a, gamma_a} = 0?

G4  How much does the triangle bound over-count on the true Clifford symbol?
    ratio (|v| + |w|) / sqrt(|v|^2 + |w|^2)   -> sqrt(2) at |v| = |w|.

G5  Camporesi-Higuchi window bookkeeping for the exact inclusion (a).
"""
import numpy as np

I2 = np.eye(2, dtype=complex)
s1 = np.array([[0, 1], [1, 0]], dtype=complex)
s2 = np.array([[0, -1j], [1j, 0]])
s3 = np.array([[1, 0], [0, -1]], dtype=complex)
PAULI = [s1, s2, s3]
rng = np.random.default_rng(11)


def kron3(a, b, c):
    return np.kron(np.kron(a, b), c)


print("G1  solve {gamma, sigma_i} = 0 (i=1,2,3) for gamma in M_2(C)")
rows = []
for s in PAULI:
    # vec(gamma s + s gamma) = (s^T (x) I + I (x) s) vec(gamma)  [column-major]
    rows.append(np.kron(s.T, I2) + np.kron(I2, s))
Amat = np.vstack(rows)
u, sv, vh = np.linalg.svd(Amat)
print("     singular values of the 12x4 constraint matrix:",
      np.array2string(sv, precision=6))
print("     nullity (sv < 1e-10): %d  ->  only gamma = 0 solves it"
      % int(np.sum(sv < 1e-10)))
print("     => NO bundle-endomorphism grading exists on a 3-manifold spinor")
print("     algebraic reason: gamma would anticommute with the volume element")
print("     omega = -i s1 s2 s3, but omega = I is central:")
om = -1j * s1 @ s2 @ s3
print("     -i*s1*s2*s3 =", np.array2string(om, precision=6))

print("")
print("G2  Clifford-doubled product symbol")
worst_ac, worst_sq, worst_nm = 0.0, 0.0, 0.0
for _ in range(500):
    v, w = rng.normal(size=3), rng.normal(size=3)
    cv = sum(v[i] * PAULI[i] for i in range(3))
    cw = sum(w[i] * PAULI[i] for i in range(3))
    A = kron3(cv, s1, I2)
    B = kron3(I2, s2, cw)
    worst_ac = max(worst_ac, float(np.max(np.abs(A @ B + B @ A))))
    S = A + B
    tgt = (v @ v + w @ w) * np.eye(8)
    worst_sq = max(worst_sq, float(np.max(np.abs(S @ S - tgt))))
    nrm = float(np.linalg.norm(S, 2))
    worst_nm = max(worst_nm, abs(nrm - np.sqrt(v @ v + w @ w)))
print("     max |{A,B}|                          = %.3e" % worst_ac)
print("     max |(A+B)^2 - (|v|^2+|w|^2) I|      = %.3e" % worst_sq)
print("     max | ||A+B|| - sqrt(|v|^2+|w|^2) |  = %.3e" % worst_nm)
print("     => C_3^(2) = 1 EXACTLY at the continuum symbol level")

print("")
print("G3  Paper 39 eq:anticomm with gamma_a := sign(D_a)")
# toy CH-like spectrum: shells n=1,2 both chiralities
D_a = np.diag([1.5, 1.5, -1.5, -1.5, 2.5, 2.5, -2.5, -2.5])
G_a = np.diag(np.sign(np.diag(D_a)))
ac = D_a @ G_a + G_a @ D_a
print("     ||{D_a, sign(D_a)}||_op = %.6f   (eq:anticomm needs 0)"
      % float(np.linalg.norm(ac, 2)))
print("     ||[D_a, sign(D_a)]||_op = %.6f   (sign(D) COMMUTES with D)"
      % float(np.linalg.norm(D_a @ G_a - G_a @ D_a, 2)))

print("")
print("G4  triangle over-count on the Clifford symbol")
print("     %8s %8s %14s %14s %10s" % ("|v|", "|w|", "||A+B|| (true)",
                                       "|v|+|w| (P39)", "ratio"))
for (nv, nw) in [(1, 1), (1, 2), (1, 0.5), (3, 3), (1, 0.01)]:
    tru = np.hypot(nv, nw)
    print("     %8.3f %8.3f %14.6f %14.6f %10.6f"
          % (nv, nw, tru, nv + nw, (nv + nw) / tru))
print("     sqrt(2) = %.6f is exactly the |v|=|w| over-count" % np.sqrt(2))

print("")
print("G5  Camporesi-Higuchi window bookkeeping (exact inclusion, part (a))")
print("     %4s %6s %28s %28s %10s" % ("2j", "n_max", "V_j (x) V_{j+1/2}",
                                       "V_j (x) V_{j-1/2}", "in window?"))
for n_max in [2, 3, 5, 8]:
    ok = True
    for twoj in range(0, n_max):                # j <= J = (n_max-1)/2
        j = twoj / 2
        n_plus, n_minus = twoj + 1, twoj
        d_plus = int((2 * j + 1) * (2 * j + 2))
        d_minus = int((2 * j + 1) * (2 * j))
        ok_p = (d_plus == n_plus * (n_plus + 1)) and n_plus <= n_max
        ok_m = (n_minus == 0) or ((d_minus == n_minus * (n_minus + 1))
                                  and n_minus <= n_max)
        ok = ok and ok_p and ok_m
        print("     %4d %6d   shell n=%d dim %3d = n(n+1)=%3s   "
              "shell n=%d dim %3d = n(n+1)=%3s   %s"
              % (twoj, n_max, n_plus, d_plus, n_plus * (n_plus + 1),
                 n_minus, d_minus,
                 n_minus * (n_minus + 1), "yes" if (ok_p and ok_m) else "NO"))
    lhs = 2 * sum(n * n for n in range(1, n_max + 1))
    rhs = 2 * n_max * (n_max + 1) * (n_max + 2) // 3
    print("     n_max=%d: dim(window span) = 2*sum n^2 = %d  <=  "
          "dim H_nmax = %d   %s   [inclusion %s]"
          % (n_max, lhs, rhs, "OK" if lhs <= rhs else "FAIL",
             "EXACT" if ok else "BROKEN"))
