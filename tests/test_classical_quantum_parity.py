"""Classical-vs-quantum parity tests for the GeoVac chemistry pipeline.

Cross-checks the continuous Level 4 multichannel + PK solver against the
balanced qubit FCI on the same molecule. The mission: any regression that
silently desyncs the two solver paths (e.g., a call-convention bug in the
R-dependence corrector that makes one bind and the other not — exactly
the Path B bug that produced Sprint W1e-Projection-Audit 2026-06-07) is
caught at the test layer instead of at the next chemistry sprint.

Marked slow because the continuous solver takes ~30s and the balanced
qubit FCI panel takes ~75s; run with ``pytest --slow``.
"""

from __future__ import annotations

import numpy as np
import pytest


@pytest.mark.slow
def test_lih_continuous_vs_balanced_qubit_both_bind() -> None:
    """LiH must bind under both the continuous Level 4 solver and the
    balanced qubit FCI, at approximately the same R_eq.

    Reference values (CHANGELOG v3.56.0 and Sprint W1e-Projection-Audit
    2026-06-07):
      - Continuous L4 + PK: E_min = -8.183 Ha at R = 3.015 bohr,
        D_e = 0.067 Ha (well depth vs R = 5.0).
      - Balanced qubit FCI: bowl-shaped on the same R panel with minimum
        at R = 3.015 bohr (panel resolution), D_e ~= 0.158 Ha (2.4x
        over-binding vs continuous).

    The Path B V_NN-double-count bug fixed in this sprint would break
    this test: pre-fix balanced gave a monotone-descending PES with
    minimum at the panel boundary, not at R = 3.015. If a future change
    re-introduces that kind of desync, this test fires.
    """
    from geovac.composed_diatomic import ComposedDiatomicSolver
    from geovac.molecular_spec import lih_spec
    from geovac.balanced_coupled import build_balanced_hamiltonian
    from geovac.coupled_composition import coupled_fci_energy

    R_panel = np.array([2.5, 3.015, 3.5, 4.0, 5.0])

    # Continuous Level 4 solver.
    solver = ComposedDiatomicSolver.LiH_ab_initio(l_max=2, verbose=False)
    solver.solve_core()
    pes_cont = solver.scan_pes(R_grid=R_panel, n_Re=300)
    E_cont = pes_cont['E_composed']

    # Balanced qubit FCI (Path A or post-fix Path B — both work after the
    # 2026-06-07 V_NN bug fix; using Path A here is robust for any caller).
    spec = lih_spec()  # default R = 3.015
    E_bal = np.zeros(len(R_panel))
    for i, R in enumerate(R_panel):
        bal = build_balanced_hamiltonian(spec, nuclei=None, R=R, verbose=False)
        E_bal[i] = coupled_fci_energy(bal, n_electrons=4, verbose=False)['E_coupled']

    # Continuous has interior minimum at R = 3.015 (index 1).
    i_min_cont = int(np.argmin(E_cont))
    assert i_min_cont == 1, (
        f"Continuous LiH PES minimum at R={R_panel[i_min_cont]}, "
        f"expected R=3.015. Panel: {dict(zip(R_panel.tolist(), E_cont.tolist()))}"
    )
    D_e_cont = E_cont[-1] - E_cont[i_min_cont]
    assert D_e_cont > 0.0, (
        f"Continuous LiH not binding: well depth = {D_e_cont} <= 0."
    )

    # Balanced qubit FCI must have interior minimum at the same R within
    # panel resolution. Pre-W1e-Projection-Audit (Path B bug), balanced
    # gave argmin at panel boundary R=5.0; the test must reject that.
    i_min_bal = int(np.argmin(E_bal))
    assert 1 <= i_min_bal <= 2, (
        f"Balanced LiH qubit FCI minimum at R={R_panel[i_min_bal]}, "
        f"expected R in [3.015, 3.5]. Panel: "
        f"{dict(zip(R_panel.tolist(), E_bal.tolist()))}.  "
        f"If at panel boundary, the Path B V_NN-double-count bug has "
        f"regressed -- check spec.R wiring in balanced_coupled.py."
    )

    D_e_bal = E_bal[-1] - E_bal[i_min_bal]
    assert D_e_bal > 0.0, (
        f"Balanced LiH qubit FCI not binding: well depth = {D_e_bal} <= 0."
    )

    # Balanced over-binds by ~2.4x relative to continuous.  Allow a wide
    # band [1.5x, 5x] so this remains a SANITY check on the relative
    # magnitudes, not a fragile pin on the precise over-binding ratio.
    ratio = D_e_bal / D_e_cont
    assert 1.5 < ratio < 5.0, (
        f"Balanced/continuous well-depth ratio = {ratio:.2f}, "
        f"expected in (1.5, 5.0).  Continuous D_e = {D_e_cont:.4f}, "
        f"balanced D_e = {D_e_bal:.4f}."
    )


@pytest.mark.slow
def test_lih_path_a_vs_path_b_bit_identical() -> None:
    """Balanced LiH must give the same FCI energy under Path A and Path B.

    Path A: ``lih_spec()`` with default R, ``build_balanced_hamiltonian(spec, R=R)``.
    Path B: ``lih_spec(R=R)`` with explicit R, ``build_balanced_hamiltonian(spec, R=R)``.

    Pre-2026-06-07 fix, Path B silently double-counted V_NN(R) - V_NN(R_default).
    Post-fix, both paths produce bit-identical energies because the corrector
    reads ``spec.R`` (set by ``hydride_spec``) instead of ``_HYDRIDE_REQ[name]``.

    This test is the regression sentinel for the V_NN double-count bug.
    """
    from geovac.molecular_spec import lih_spec
    from geovac.balanced_coupled import build_balanced_hamiltonian
    from geovac.coupled_composition import coupled_fci_energy

    R_panel = [2.5, 3.5, 5.0]  # three points; 3.015 is trivially equal

    for R in R_panel:
        spec_a = lih_spec()  # default R = 3.015
        spec_b = lih_spec(R=R)

        bal_a = build_balanced_hamiltonian(spec_a, nuclei=None, R=R,
                                           verbose=False)
        bal_b = build_balanced_hamiltonian(spec_b, nuclei=None, R=R,
                                           verbose=False)

        e_a = coupled_fci_energy(bal_a, n_electrons=4,
                                 verbose=False)['E_coupled']
        e_b = coupled_fci_energy(bal_b, n_electrons=4,
                                 verbose=False)['E_coupled']

        # Bit-identical at machine precision (no fitting tolerance:  if these
        # disagree, a regression has been introduced in the corrector).
        assert abs(e_a - e_b) < 1e-10, (
            f"Path A and Path B balanced LiH disagree at R={R}: "
            f"E_A={e_a:.10f} vs E_B={e_b:.10f}, diff={e_a - e_b:+.3e}.  "
            f"Check spec.R wiring in balanced_coupled.py's V_NN corrector."
        )
