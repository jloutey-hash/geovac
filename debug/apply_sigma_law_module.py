"""Add the contraction-window apparatus to geovac/sturmian_sigma_law.py.

Three functions, so the Paper 60 claims they back recompute from TRACKED code
rather than from the prunable debug/ tree (the v5.10.17 coverage lesson):

  theta2_band_matrix(n)        closed-form <theta^2> in the sine band; symbol
                               (pi-chi)^2, Toeplitz minus Hankel.  No quadrature.
  weighted_blocks(s, n, W)     intra/cross blocks under a radial Fock-sphere
                               weight W; W = 1 must reproduce (I, sw_cross_block).
  generalized_sigma_max(A, B)  largest generalized singular value.

Idempotent.
"""
from __future__ import annotations

import sys

MOD = "geovac/sturmian_sigma_law.py"
ANCHOR = "def sigma_law_report("
MARKER = "def theta2_band_matrix("

NEW = '''def theta2_band_matrix(n: int) -> np.ndarray:
    """<theta^2> in the sine band, in closed form:  no quadrature, no Bessel.

    T_ab = (2/pi) int_0^pi sin(a chi) sin(b chi) (pi - chi)^2 dchi, evaluated
    exactly as a Toeplitz-minus-Hankel matrix with c_j = 2/j^2 (j != 0),
    c_0 = pi^2/3, i.e. symbol (pi - chi)^2 = theta^2:

        T_ab = 2/(a-b)^2 - 2/(a+b)^2   (a != b),
        T_aa = pi^2/3 - 1/(2 a^2).

    Its least eigenvalue is the minimal mean-square spread of a band-limited
    function about the Fock pole, and n^2 lam_min -> pi^2 -- the
    Kac-Murdock-Szego constant obtained with no Bessel function present, which
    is what shows that the pi^2 of the conditioning law is a TRUNCATION
    constant rather than Bessel content.  Backing test:
    tests/test_paper60_contraction_window.py.
    """
    a = np.arange(1, n + 1)
    A = a[:, None]
    B = a[None, :]
    d = A - B
    s = A + B
    T = 2.0 / np.where(d == 0, 1, d) ** 2 - 2.0 / s ** 2
    np.fill_diagonal(T, np.pi ** 2 / 3.0 - 1.0 / (2.0 * a ** 2))
    return T


def weighted_blocks(s: float, n: int, weight=None, M: int = 200001
                    ) -> Tuple[np.ndarray, np.ndarray]:
    """Intra/cross blocks under a radial weight W on the Fock sphere.

    A_ab = (2/pi) int sin(a chi) sin(b chi) W(chi) dchi              (intra)
    B_ab = (2/pi) int sin(a chi) sin(b chi) W(chi) j0(s cot(chi/2)) dchi (cross)

    weight=None means W = 1, for which (A, B) reproduces (I, sw_cross_block).
    The generalized symbol is the quotient W*j0/W = j0, so the conditioning law
    is predicted to be independent of any smooth positive W -- measured, with
    the vanishing-W control, in tests/test_paper60_contraction_window.py.
    """
    chi = np.linspace(1e-9, np.pi, M)
    wq = np.full(M, chi[1] - chi[0])
    wq[0] *= 0.5
    wq[-1] *= 0.5
    Wv = np.ones_like(chi) if weight is None else np.asarray(weight(chi), float)
    x = s / np.tan(chi / 2.0)
    sym = np.ones_like(x)
    nz = x != 0.0
    sym[nz] = np.sin(x[nz]) / x[nz]
    a = np.arange(1, n + 1)
    smat = np.sin(np.outer(a, chi))
    A = (2.0 / np.pi) * (smat * (wq * Wv)) @ smat.T
    B = (2.0 / np.pi) * (smat * (wq * Wv * sym)) @ smat.T
    return A, B


def generalized_sigma_max(A: np.ndarray, B: np.ndarray) -> float:
    """Largest generalized singular value: spectral radius of A^-1/2 B A^-1/2."""
    ev, Q = np.linalg.eigh(A)
    ev = np.maximum(ev, 1e-300)
    half = Q @ np.diag(ev ** -0.5) @ Q.T
    return float(np.max(np.abs(np.linalg.eigvalsh(half @ B @ half))))


'''


def main() -> int:
    with open(MOD, encoding="utf-8") as fh:
        text = fh.read()
    if MARKER in text:
        print("ALREADY APPLIED -- marker present; nothing done.")
        return 1
    if text.count(ANCHOR) != 1:
        print(f"ANCHOR not unique (count={text.count(ANCHOR)}); aborting.")
        return 2
    text = text.replace(ANCHOR, NEW + ANCHOR)
    with open(MOD, "w", encoding="utf-8") as fh:
        fh.write(text)
    print("applied: 3 functions added to geovac/sturmian_sigma_law.py")
    return 0


if __name__ == "__main__":
    sys.exit(main())
