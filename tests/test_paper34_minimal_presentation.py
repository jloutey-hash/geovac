"""Backing for Paper 34's minimal-presentation remark (rem:minimal_presentation).

Three legs:

1. IDENTITY -- the atomic one-electron Hamiltonian in the shared-k Coulomb-Sturmian
   basis carries no information beyond (integer labels, overlap metric, charge, one
   scale):  h1 = k^2 (I - S/2) - Z k diag(1/n).  Verified against the INDEPENDENT
   grid-quadrature route (geovac.transcorrelated_sturmian.build_one_body: gradient-form
   kinetic + numerical 1/r, which knows nothing of the identity).

2. TEETH -- a deliberately wrong closed form (factor I - S instead of I - S/2) must
   fail by many orders more than the true one, so the tolerance discriminates.

3. BOUNDARY -- the two-electron ERI tensor must NOT reduce to the same data: its best
   fit by Kronecker words in the metric leaves a large residual (measured 0.388 at
   ns=4), while the same fitting machinery reproduces a genuinely metric-generated
   tensor to ~1e-12 (self-validation).  The derived/independent boundary sits exactly
   at the electron-electron channel.

Provenance: debug/minimal_rep_h1_identity.py (2026-08-27) and
debug/sprint_minimal_presentation_memo.md.
"""
import numpy as np
import pytest

from geovac import transcorrelated_sturmian as TC


def _build(ns, k, Z, Ng=500):
    r, wr = TC.make_grid(k, Ng=Ng)
    S, h1, Rtab, W = TC.build_one_body(ns, r, wr, k, Z)
    return r, wr, S, h1, Rtab, W


def _closed_form(S, k, Z, ns):
    n = np.arange(1, ns + 1)
    return k * k * (np.eye(ns) - S / 2.0) - Z * k * np.diag(1.0 / n)


@pytest.mark.parametrize("ns,k,Z", [(3, 1.7, 2.0), (4, 2.4, 3.0), (5, 2.0, 3.0)])
def test_h1_is_metric_labels_scale(ns, k, Z):
    """h1 = k^2 (I - S/2) - Z k diag(1/n) at grid precision (independent route)."""
    _, _, S, h1, _, _ = _build(ns, k, Z)
    pred = _closed_form(S, k, Z, ns)
    rel = np.abs(h1 - pred).max() / np.abs(h1).max()
    assert rel < 1e-6, f"identity broken: rel dev {rel:.2e}"


def test_identity_has_teeth():
    """A wrong coefficient (I - S instead of I - S/2) must fail loudly."""
    ns, k, Z = 4, 2.4, 3.0
    _, _, S, h1, _, _ = _build(ns, k, Z)
    n = np.arange(1, ns + 1)
    wrong = k * k * (np.eye(ns) - S) - Z * k * np.diag(1.0 / n)
    rel_true = np.abs(h1 - _closed_form(S, k, Z, ns)).max() / np.abs(h1).max()
    rel_wrong = np.abs(h1 - wrong).max() / np.abs(h1).max()
    assert rel_wrong > 1e-1, "wrong form should deviate at O(1)"
    assert rel_wrong / max(rel_true, 1e-300) > 1e4, "tolerance does not discriminate"


def test_eri_is_not_metric_generated():
    """The two-electron tensor refuses the reduction; the fitter itself is validated."""
    ns, k, Z, Ng = 4, 2.4, 3.0, 500
    r, wr, S, h1, Rtab, W = _build(ns, k, Z, Ng=Ng)
    Km = TC.build_kernels(r, 0.7, nx=64)
    eri, _, _ = TC.two_body(ns, Rtab, W, Km)
    E = eri.reshape(ns * ns, ns * ns)

    P = [np.eye(ns), S, S @ S, S @ S @ S]
    M = np.array([np.kron(A, B).ravel() for A in P for B in P]).T

    def resid(target):
        coef, *_ = np.linalg.lstsq(M, target.ravel(), rcond=None)
        return float(np.linalg.norm(target.ravel() - M @ coef)
                     / np.linalg.norm(target.ravel()))

    # self-validation: a metric-generated tensor must be reproduced ~exactly
    control = 0.3 * np.kron(S, S @ S) + 1.7 * np.kron(np.eye(ns), S)
    assert resid(control) < 1e-10, "fitting machinery broken"

    # the physical ERI must NOT reduce
    r_eri = resid(E)
    assert r_eri > 0.2, f"ERI unexpectedly metric-generated (residual {r_eri:.3f})"
