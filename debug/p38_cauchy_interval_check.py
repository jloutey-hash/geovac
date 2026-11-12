"""Independent route: eigenvalues of (-Delta)^(1/2) on (-1,1), Dirichlet exterior.
Galerkin on phi_m = sqrt(1-x^2) U_m(x), using the identity (-Delta)^(1/2) phi_m = (m+1) U_m on (-1,1).
K_ml = <phi_m, (m+1)U_l...> = (l+1) int sqrt(1-x^2) U_m U_l = (l+1) pi/2 delta; M_ml = int (1-x^2) U_m U_l."""
import numpy as np
from scipy.linalg import eigh
from numpy.polynomial import chebyshev as C

def run(M_terms, parity):
    ms = [m for m in range(M_terms) if m % 2 == parity]
    # Gauss-Chebyshev (2nd kind) quadrature: int sqrt(1-x^2) f = sum w f, exact for poly deg < 2Q
    Q = 2 * M_terms + 10
    i = np.arange(1, Q + 1); th = i * np.pi / (Q + 1); x = np.cos(th); w = np.pi / (Q + 1) * np.sin(th) ** 2
    U = np.array([np.sin((m + 1) * th) / np.sin(th) for m in ms])
    Mm = (U * w * np.sqrt(1 - x**2) ** 0) @ (U * (1 - x**2) / np.sqrt(1 - x**2) ** 0).T
    # careful: int (1-x^2) U_m U_l dx = int sqrt(1-x^2) * sqrt(1-x^2) U_m U_l dx -> weight w times sqrt(1-x^2)
    Mm = (U * w * np.sqrt(1 - x**2)) @ U.T
    K = np.diag([(m + 1) * np.pi / 2 for m in ms])
    return eigh(K, Mm, eigvals_only=True)[:2]

for M in (40, 80, 160, 320):
    print(M, "even:", run(M, 0)[0], " odd:", run(M, 1)[0])
