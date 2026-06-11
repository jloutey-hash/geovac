"""
He-4 Electric Dipole Polarizability (Compact Nucleus Control)
=============================================================
Uses the Paper 23 He-4 Hamiltonian (2p + 2n, Minnesota NN + Coulomb)
to compute alpha_E. He-4 should be SMALL (~0.07 fm^3) since it's
tightly bound (B = 28.3 MeV). This is the control case for the
deuteron polarizability computation.

For He-4 (Z=2, N=2), the E1 operator is the isovector dipole:
    D_z = (e/A) * sum_i t_z(i) * z_i
where t_z is the isospin projection (+1/2 for protons, -1/2 for neutrons).
For the proton-only part:
    D_z^(p) = sum_{protons} z_i - (Z/A) * sum_{all} z_i
The second term subtracts center-of-mass motion. For a 1-body E1
operator in the (proton SD) x (neutron SD) basis, we compute the
proton z-coordinate sum minus the CoM correction.

Simpler approach: since we're in the HO basis, we can use the
Tassie-Barker prescription: the E1 response is proportional to the
isovector dipole strength function, computable from the FCI spectrum.
"""
import sys, io, json
import numpy as np

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')
sys.path.insert(0, '.')

from geovac.nuclear.nuclear_hamiltonian import (
    build_he4_hamiltonian, diagonalize_he4, enumerate_sp_states,
    DeuteronSpec, SPState, _enumerate_slater_dets, _sd_phase,
)

HBAR_C = 197.3269804  # MeV*fm
E2_MEV_FM = 1.4399764  # e^2 in MeV*fm


def ho_radial_r_me(nr1, l1, nr2, l2, b_fm):
    """Radial matrix element <n_r1 l1 | r | n_r2 l2> for 3D HO."""
    if l1 == l2 + 1:
        if nr1 == nr2:
            return b_fm * np.sqrt(nr2 + l2 + 1.5)
        elif nr1 == nr2 - 1:
            return b_fm * np.sqrt(float(nr2))
        return 0.0
    elif l1 == l2 - 1:
        return ho_radial_r_me(nr2, l2, nr1, l1, b_fm)
    return 0.0


def angular_cos_theta(l1, m1, l2, m2):
    """Angular matrix element <l1 m1 | cos(theta) | l2 m2>."""
    if m1 != m2:
        return 0.0
    m = m1
    if l1 == l2 + 1:
        return np.sqrt(((l2 + 1)**2 - m**2) / ((2*l2 + 1) * (2*l2 + 3)))
    elif l1 == l2 - 1:
        return np.sqrt((l2**2 - m**2) / ((2*l2 - 1) * (2*l2 + 1)))
    return 0.0


def sp_dipole_z_matrix(states, b_fm):
    """One-body dipole z matrix in the single-particle basis.
    <alpha | z | beta> = <alpha | r cos(theta) | beta>"""
    n = len(states)
    D = np.zeros((n, n))
    for a, sa in enumerate(states):
        for c, sc in enumerate(states):
            if sa.m_l != sc.m_l or sa.m_s != sc.m_s:
                continue
            if abs(sa.l - sc.l) != 1:
                continue
            r_me = ho_radial_r_me(sa.n_r, sa.l, sc.n_r, sc.l, b_fm)
            ang_me = angular_cos_theta(sa.l, sa.m_l, sc.l, sc.m_l)
            D[a, c] = r_me * ang_me
    return D


def build_he4_dipole_z(states_p, states_n, b_fm):
    """Build the E1 dipole operator D_z in the He-4 FCI basis.

    For He-4, the isovector E1 operator is:
        D_z = sum_{protons} z_i - (Z/A) * sum_{all} z_i
    The CoM subtraction (Z/A = 1/2) removes spurious center-of-mass excitations.

    In the SD basis: D_z acts on proton SDs via one-body z operator,
    with neutron SDs as spectators, minus 1/2 * (proton z + neutron z).
    """
    n_p = len(states_p)
    n_n = len(states_n)

    sds_p = _enumerate_slater_dets(n_p, 2)
    sds_n = _enumerate_slater_dets(n_n, 2)
    dim_p = len(sds_p)
    dim_n = len(sds_n)
    dim = dim_p * dim_n

    sd_p_idx = {sd: i for i, sd in enumerate(sds_p)}
    sd_n_idx = {sd: j for j, sd in enumerate(sds_n)}

    D_z_sp_p = sp_dipole_z_matrix(states_p, b_fm)
    D_z_sp_n = sp_dipole_z_matrix(states_n, b_fm)

    D_z = np.zeros((dim, dim))

    # Proton contribution: (1 - Z/A) * sum_protons z_i = (1/2) * sum_protons z_i
    # Neutron contribution: -(Z/A) * sum_neutrons z_i = -(1/2) * sum_neutrons z_i
    # Equivalent to isovector: D_z = (1/2) * (sum_p z_p - sum_n z_n)

    # Proton one-body part
    for ip, sd_p_bra in enumerate(sds_p):
        for a in range(n_p):
            for c in range(n_p):
                if abs(D_z_sp_p[a, c]) < 1e-15:
                    continue
                phase, new_sd = _sd_phase(sd_p_bra, a, c)
                if phase == 0:
                    continue
                jp = sd_p_idx.get(new_sd)
                if jp is None:
                    continue
                for jn in range(dim_n):
                    bra = ip * dim_n + jn
                    ket = jp * dim_n + jn
                    D_z[bra, ket] += 0.5 * phase * D_z_sp_p[a, c]

    # Neutron one-body part (negative, isovector)
    for jn, sd_n_bra in enumerate(sds_n):
        for b in range(n_n):
            for d in range(n_n):
                if abs(D_z_sp_n[b, d]) < 1e-15:
                    continue
                phase, new_sd = _sd_phase(sd_n_bra, b, d)
                if phase == 0:
                    continue
                kn = sd_n_idx.get(new_sd)
                if kn is None:
                    continue
                for ip in range(dim_p):
                    bra = ip * dim_n + jn
                    ket = ip * dim_n + kn
                    D_z[bra, ket] -= 0.5 * phase * D_z_sp_n[b, d]

    return D_z


def compute_he4_polarizability(hw_values=None):
    """Compute He-4 alpha_E at various hw values."""
    if hw_values is None:
        hw_values = [10.0, 15.0, 20.0, 25.0, 30.0, 40.0]

    print("=" * 72)
    print("He-4 ELECTRIC DIPOLE POLARIZABILITY (COMPACT NUCLEUS CONTROL)")
    print("from Paper 23 Minnesota NN + Coulomb Hamiltonian")
    print("=" * 72)

    # Literature: alpha_E(He-4) = 0.0725(15) fm^3 (Pachucki & Miskimen 2003)
    # Some refs: 0.076(8) fm^3 (Ji et al. 2003)
    # The key point: ~10x smaller than deuteron's 0.633 fm^3
    alpha_E_lit = 0.073  # fm^3 (approximate)
    alpha_E_d = 0.633   # deuteron for comparison

    print(f"\nLiterature: alpha_E(He-4) ~ {alpha_E_lit:.3f} fm^3")
    print(f"           alpha_E(d) = 0.633 fm^3 (for comparison)")
    print(f"           Ratio d/He-4 ~ {alpha_E_d/alpha_E_lit:.0f}x (compact vs extended)")

    results = []

    for hw in hw_values:
        print(f"\n  hw = {hw:.1f} MeV ...", end=" ", flush=True)
        try:
            data = diagonalize_he4(N_shells=2, hw=hw, include_coulomb=True)
        except Exception as e:
            print(f"FAILED ({e})")
            continue

        H = data['H_data']['H_matrix']
        evals, evecs = np.linalg.eigh(H)
        gs = evecs[:, 0]
        E_gs = evals[0]

        states_p = data['H_data']['states_p']
        states_n = data['H_data']['states_n']
        b_fm = data['metadata']['b_fm']
        dim = len(evals)

        D_z = build_he4_dipole_z(states_p, states_n, b_fm)
        assert np.allclose(D_z, D_z.T, atol=1e-12), "D_z not symmetric!"

        D_gs = D_z @ gs

        alpha_E = 0.0
        ewsr = 0.0
        n_contributing = 0

        for n in range(1, dim):
            d_me = np.dot(evecs[:, n], D_gs)
            delta_E = evals[n] - E_gs
            if delta_E < 1e-10:
                continue
            # alpha_E = 2 * e^2 * sum |<n|D_z|0>|^2 / (E_n - E_0)
            # D_z already includes the isovector 1/2 factor and is in fm
            # So alpha_E = 2 * e^2 * |d_me|^2 / delta_E [fm^3]
            contribution = 2 * E2_MEV_FM * d_me**2 / delta_E
            alpha_E += contribution
            ewsr += d_me**2 * delta_E
            if abs(d_me) > 1e-10:
                n_contributing += 1

        m_N = 938.918
        ewsr_trk = HBAR_C**2 * 2 / (2 * m_N * 4)  # Z=2, A=4: Z*hbar^2/(2*m*A)

        print(f"done")
        print(f"    E_gs = {E_gs:.4f} MeV, FCI dim = {dim}")
        print(f"    alpha_E = {alpha_E:.6f} fm^3")
        print(f"    EWSR fraction = {ewsr/ewsr_trk:.3f}")
        print(f"    Contributing transitions: {n_contributing}")

        # Compare to literature
        if alpha_E > 0:
            residual = (alpha_E - alpha_E_lit) / alpha_E_lit * 100
            print(f"    vs literature: {residual:+.1f}%")
        else:
            residual = float('nan')

        results.append({
            'hw': hw, 'b_fm': b_fm, 'E_gs': E_gs,
            'alpha_E_fm3': alpha_E, 'dim': dim,
            'ewsr_fraction': ewsr / ewsr_trk,
            'residual_pct': residual,
            'n_contributing': n_contributing,
        })

    # Summary comparison with deuteron
    print(f"\n{'='*72}")
    print(f"SUMMARY: He-4 vs Deuteron Polarizability")
    print(f"{'='*72}")
    print(f"\n{'hw':>6} {'alpha_E(He4)':>14} {'alpha_E(d)':>12} {'Ratio d/He4':>12}")
    print(f"-" * 50)

    # Load deuteron results
    try:
        with open('debug/data/deuteron_polarizability.json') as f:
            d_results = json.load(f)
        d_by_hw = {r['hw']: r['alpha_E_fm3'] for r in d_results if r['N_shells'] == 2}
    except Exception:
        d_by_hw = {}

    for r in results:
        d_val = d_by_hw.get(r['hw'], float('nan'))
        ratio = d_val / r['alpha_E_fm3'] if r['alpha_E_fm3'] > 0 else float('nan')
        print(f"{r['hw']:>6.1f} {r['alpha_E_fm3']:>14.6f} {d_val:>12.4f} {ratio:>12.1f}x")

    print(f"\nLiterature ratio: {alpha_E_d/alpha_E_lit:.1f}x")
    print(f"(Deuteron is {alpha_E_d/alpha_E_lit:.0f}x more polarizable than He-4)")

    # Physics interpretation
    print(f"\n--- PHYSICS ---")
    print(f"  He-4 binding: B = 28.3 MeV (tightly bound)")
    print(f"  Deuteron binding: B = 2.22 MeV (weakly bound)")
    print(f"  Ratio B(He4)/B(d) = {28.3/2.22:.0f}x")
    print(f"  Polarizability scales as ~1/B^2 for nuclear systems")
    print(f"  Predicted ratio: (28.3/2.22)^2 = {(28.3/2.22)**2:.0f}x")
    print(f"  This is why He-4 nuclear structure contributes negligibly")
    print(f"  to muonic He-4 Lamb shift (+4 meV) vs muonic deuterium (+1.69 meV)")

    with open('debug/data/he4_polarizability.json', 'w') as f:
        json.dump(results, f, indent=2, default=float)
    print(f"\nResults saved to debug/data/he4_polarizability.json")

    return results


if __name__ == '__main__':
    compute_he4_polarizability()
