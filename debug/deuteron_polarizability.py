"""
Deuteron Electric Dipole Polarizability from Paper 23 Hamiltonian
================================================================
Computes alpha_E(deuteron) using sum-over-states on the FCI spectrum
of the Minnesota NN potential deuteron Hamiltonian.

The electric dipole operator couples to PROTON position only
(neutron is uncharged). In the HO basis, the E1 operator is:
    D_z = e * z_proton = e * r_p * cos(theta_p)
which has selection rules Delta_l = +/-1, Delta_m_l = 0 (for z-component),
Delta_m_s = 0, and acts only on the proton register.

Polarizability (static, isotropic):
    alpha_E = (2/3) * sum_{n>0} |<n|D|0>|^2 / (E_n - E_0)

where the 2/3 comes from averaging over the three Cartesian components
and using |<n|D_z|0>|^2 = (1/3) * |<n|D|0>|^2 for an isotropic system.

Actually for the deuteron (J=1 ground state), the static polarizability is:
    alpha_E = (2e^2) * sum_{n>0} |<n|r_p Y_10(Omega_p)|0>|^2 / (E_n - E_0)

Units: We work in MeV (energies) and fm (distances), so alpha_E comes
out in fm^3 (the standard nuclear polarizability unit).
"""
import sys, io, json
import numpy as np

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')
sys.path.insert(0, '.')

from geovac.nuclear.nuclear_hamiltonian import (
    build_deuteron_hamiltonian, diagonalize_deuteron, enumerate_sp_states,
    DeuteronSpec, SPState,
)

# Constants
HBAR_C = 197.3269804  # MeV*fm
E2_MEV_FM = 1.4399764  # e^2 in MeV*fm (= alpha * hbar * c)


def build_dipole_z_proton(states_p, states_n, b_fm):
    """Build the z-component of the electric dipole operator D_z = e * z_p
    in the FCI basis |p_i, n_j>.

    D_z acts ONLY on the proton register. In the HO basis:
    <n'_r l' m'_l m'_s | r cos(theta) | n_r l m_l m_s>
    = delta(m_l', m_l) * delta(m_s', m_s) * R_{n'_r l', n_r l} * A_{l' m_l, l m_l}

    where R is the radial matrix element of r between HO states,
    and A is the angular matrix element of cos(theta) = C^1_0 (Racah).

    Selection rules: Delta_l = +/-1, Delta_m_l = 0, Delta_m_s = 0.
    """
    n_p = len(states_p)
    n_n = len(states_n)
    dim = n_p * n_n

    D_z = np.zeros((dim, dim))

    for i, pi_state in enumerate(states_p):
        for k, pk_state in enumerate(states_p):
            # Selection rules
            if pi_state.m_l != pk_state.m_l:
                continue
            if pi_state.m_s != pk_state.m_s:
                continue
            if abs(pi_state.l - pk_state.l) != 1:
                continue

            # Radial matrix element <n'_r l' | r | n_r l> in units of b (HO length)
            r_me = ho_radial_r_matrix_element(
                pi_state.n_r, pi_state.l,
                pk_state.n_r, pk_state.l,
                b_fm
            )

            # Angular matrix element <l' m_l | cos(theta) | l m_l>
            # = <l' m_l | C^1_0 | l m_l> (Racah spherical harmonic)
            # Using the standard formula:
            # <l+1, m | cos(theta) | l, m> = sqrt((l+1)^2 - m^2) / sqrt((2l+1)(2l+3))
            # <l-1, m | cos(theta) | l, m> = sqrt(l^2 - m^2) / sqrt((2l-1)(2l+1))
            ang_me = angular_cos_theta(pi_state.l, pi_state.m_l,
                                       pk_state.l, pk_state.m_l)

            if abs(r_me * ang_me) < 1e-15:
                continue

            # Fill FCI matrix: D_z acts on proton only, neutron is spectator
            for j in range(n_n):
                idx_bra = i * n_n + j
                idx_ket = k * n_n + j
                D_z[idx_bra, idx_ket] += r_me * ang_me

    return D_z


def ho_radial_r_matrix_element(nr1, l1, nr2, l2, b_fm):
    """Radial matrix element <n_r1 l1 | r | n_r2 l2> for 3D HO basis.

    In units of b (HO length parameter).

    Using the standard recursion for HO radial matrix elements:
    <n_r, l+1 | r | n_r, l> = sqrt((n_r + l + 3/2)) * b  [raising l by 1, same n_r]
    <n_r-1, l+1 | r | n_r, l> = sqrt(n_r) * b  [lowering n_r, raising l]

    More precisely, the matrix elements of r between HO states are:
    <n'_r l' | r | n_r l> with l' = l+1:
        = b * sqrt(2*n_r + 2*l + 3) / 2  if n'_r = n_r
        ... (general formula involves Moshinsky-like recursion)

    For the simplest cases at small N_shells, I'll use the explicit formula.
    """
    # Use the creation/annihilation operator approach
    # r = b * (a^dag + a) / sqrt(2) in 1D, but in 3D with angular momentum:
    # <n'_r l+1 | r | n_r l> = b * delta(n'_r, n_r) * sqrt(n_r + l + 3/2)
    #                         + b * delta(n'_r, n_r - 1) * sqrt(n_r)
    # (for l' = l + 1)
    # <n'_r l-1 | r | n_r l> = b * delta(n'_r, n_r) * sqrt(n_r + l + 1/2)
    #                         + b * delta(n'_r, n_r + 1) * sqrt(n_r + 1)
    # (for l' = l - 1)

    # Ref: Suhonen "From Nucleons to Nucleus" Table A.4, or
    #      de-Shalit & Talmi "Nuclear Shell Theory" Ch. 13

    if l1 == l2 + 1:
        # bra has l' = l+1, ket has l
        n_r_bra, l_bra = nr1, l1
        n_r_ket, l_ket = nr2, l2
        # <n'_r, l+1 | r | n_r, l>
        if n_r_bra == n_r_ket:
            return b_fm * np.sqrt(n_r_ket + l_ket + 1.5)
        elif n_r_bra == n_r_ket - 1:
            return b_fm * np.sqrt(float(n_r_ket))
        else:
            return 0.0
    elif l1 == l2 - 1:
        # bra has l' = l-1, ket has l → use Hermitian conjugate
        # <n_r1, l-1 | r | n_r2, l> = <n_r2, l | r | n_r1, l-1>*
        # = ho_radial_r_matrix_element(n_r2, l2, n_r1, l1, b_fm)
        return ho_radial_r_matrix_element(nr2, l2, nr1, l1, b_fm)
    else:
        return 0.0


def angular_cos_theta(l1, m1, l2, m2):
    """Angular matrix element <l1 m1 | cos(theta) | l2 m2>.

    cos(theta) = sqrt(4pi/3) * Y_10 = C^1_0 (Racah).

    Selection rules: m1 = m2, |l1 - l2| = 1.

    <l+1, m | cos(theta) | l, m> = sqrt(((l+1)^2 - m^2) / ((2l+1)*(2l+3)))
    <l-1, m | cos(theta) | l, m> = sqrt((l^2 - m^2) / ((2l-1)*(2l+1)))
    """
    if m1 != m2:
        return 0.0
    m = m1
    if l1 == l2 + 1:
        return np.sqrt(((l2 + 1)**2 - m**2) / ((2*l2 + 1) * (2*l2 + 3)))
    elif l1 == l2 - 1:
        return np.sqrt((l2**2 - m**2) / ((2*l2 - 1) * (2*l2 + 1)))
    else:
        return 0.0


def compute_polarizability(hw_values=None):
    """Compute deuteron alpha_E at various hw values."""
    if hw_values is None:
        hw_values = [5.0, 7.0, 10.0, 12.0, 15.0, 20.0, 25.0, 30.0, 40.0]

    print("=" * 72)
    print("DEUTERON ELECTRIC DIPOLE POLARIZABILITY")
    print("from Paper 23 Minnesota NN Hamiltonian")
    print("=" * 72)

    # Literature values for comparison (in fm^3)
    # Experimental: alpha_E(d) = 0.6328(17) fm^3 (Griesshammer et al. 2012)
    # Chiral EFT:   alpha_E(d) = 0.631(5) fm^3 (Beane et al. 2000)
    # The muD Lamb shift uses polarizability in meV:
    # Delta_E_pol(muD) = +1.690(20) meV (CGV 2014)
    # Delta_E_pol(muD) = +1.94 meV (CV 2011)
    # Conversion: Delta_E_pol = (2*alpha_E / (3*a_mu^3)) * const
    # More precisely, for muonic atoms:
    # Delta_E(pol) = alpha_E * <r^{-4}> * correction factors
    # The connection between fm^3 and meV is through the muonic wavefunction

    alpha_E_exp = 0.6328  # fm^3 (experimental)
    alpha_E_eft = 0.631   # fm^3 (chiral EFT)

    print(f"\nLiterature: alpha_E(d) = {alpha_E_exp:.4f}(17) fm^3 (Griesshammer 2012)")
    print(f"           alpha_E(d) = {alpha_E_eft:.3f}(5) fm^3 (chiral EFT)")

    results = []

    for N_shells in [2, 3]:
        print(f"\n{'='*60}")
        print(f"N_shells = {N_shells}")
        print(f"{'='*60}")

        for hw in hw_values:
            try:
                data = diagonalize_deuteron(N_shells=N_shells, hw=hw)
            except Exception as e:
                print(f"  hw={hw:.1f}: FAILED ({e})")
                continue

            evals = data['eigenvalues']
            evecs = data['H_data']['H_matrix']  # need eigenvectors from eigh
            H = data['H_data']['H_matrix']
            evals_full, evecs_full = np.linalg.eigh(H)

            gs = evecs_full[:, 0]
            E_gs = evals_full[0]

            states_p = data['H_data']['states_p']
            states_n = data['H_data']['states_n']
            b_fm = data['H_data']['metadata']['b_fm']

            # Build dipole operator
            D_z = build_dipole_z_proton(states_p, states_n, b_fm)

            # Check D_z is Hermitian
            assert np.allclose(D_z, D_z.T), "D_z not symmetric!"

            # Compute transition dipole moments
            # D_z|0> → project onto excited states
            D_gs = D_z @ gs  # D_z applied to ground state
            dim = len(evals_full)

            # Sum-over-states polarizability
            # alpha_E = (2*e^2/3) * sum_{n>0} |<n|D_z|0>|^2 / (E_n - E_0)
            # The factor 2/3 comes from: we compute D_z (one component),
            # and alpha_E = (2/3) * sum |<n|r|0>|^2 / Delta_E
            # Since |<n|r|0>|^2 = |<n|D_x|0>|^2 + |<n|D_y|0>|^2 + |<n|D_z|0>|^2
            # and by isotropy each is 1/3, so sum_n |<n|D_z|0>|^2 / Delta_E
            # = (1/3) * sum_n |<n|r|0>|^2 / Delta_E
            # Hence alpha_E = 2 * sum_n |<n|D_z|0>|^2 / (E_n - E_0)
            # (the factor 2 = 2e^2 where e^2 is already in the units)

            # Actually: alpha_E = 2 * sum_{n>0} |<n|D_z|0>|^2 / (E_n - E_0)
            # where D_z is in fm and E in MeV, giving alpha_E in fm^2/MeV
            # Need to convert: alpha_E [fm^3] = (e^2) * sum |<n|z_p|0>|^2 / (E_n - E_0)
            # with e^2 = 1.44 MeV*fm, z_p in fm, E in MeV
            # → alpha_E [fm^3] = 1.44 * sum [fm^2 / MeV] = fm^3

            # For isotropic system with J=1 ground state:
            # alpha_E = (2*e^2/3) * sum_{n>0} |<n||r||0>|^2 / (3*(2J+1)) / (E_n - E_0)
            # This gets complicated with reduced matrix elements.
            # Simpler: compute the response to a uniform field.

            # For the simplest sum-over-states:
            # alpha_E = 2 * sum_{n>0} |<n|D_z|0>|^2 / (E_n - E_0)
            # where D_z = z_p (proton position) in fm, E in MeV
            # This gives alpha_E in fm^2 / MeV * fm = fm^3 / MeV ← wrong units
            # Need alpha_E in fm^3. The correct formula:
            # alpha_E = 2 * e^2 * sum_{n>0} |<n|z_p|0>|^2 / (E_n - E_0)
            # = 2 * 1.44 MeV*fm * [fm^2 / MeV] = fm^3 ← correct!

            alpha_E = 0.0
            ewsr = 0.0
            n_contributing = 0
            dominant_transitions = []

            for n in range(1, dim):
                d_me = np.dot(evecs_full[:, n], D_gs)  # <n|D_z|0> in fm
                delta_E = evals_full[n] - E_gs  # MeV
                if delta_E < 1e-10:
                    continue
                contribution = 2 * E2_MEV_FM * d_me**2 / delta_E
                alpha_E += contribution
                ewsr += d_me**2 * delta_E  # energy-weighted sum rule
                if abs(d_me) > 1e-10:
                    n_contributing += 1
                    if abs(contribution) > 0.001:
                        dominant_transitions.append({
                            'n': n, 'E_n': evals_full[n],
                            'delta_E': delta_E, 'd_me': d_me,
                            'contribution': contribution,
                        })

            # EWSR check: sum_n (E_n - E_0) |<n|D_z|0>|^2 = (Z*hbar^2)/(2*m_N) * A
            # Thomas-Reiche-Kuhn: = (hbar^2 * Z * A) / (2 * m_N * A) = hbar^2*Z/(2*m_N)
            # For deuteron: Z=1, A=2, m_N = 938.918 MeV/c^2 (average nucleon)
            m_N = 938.918  # MeV
            ewsr_trk = HBAR_C**2 / (2 * m_N)  # MeV * fm^2

            print(f"\n  hw = {hw:.1f} MeV, b = {b_fm:.3f} fm")
            print(f"    E_gs = {E_gs:.4f} MeV")
            print(f"    FCI dim = {dim}")
            print(f"    alpha_E = {alpha_E:.4f} fm^3")
            print(f"    EWSR = {ewsr:.4f} MeV*fm^2 (TRK = {ewsr_trk:.4f})")
            print(f"    EWSR fraction = {ewsr/ewsr_trk:.3f}")
            print(f"    Contributing transitions: {n_contributing}")

            if dominant_transitions:
                dominant_transitions.sort(key=lambda x: -abs(x['contribution']))
                print(f"    Top transitions:")
                for t in dominant_transitions[:5]:
                    print(f"      n={t['n']:>3}, Delta_E={t['delta_E']:.3f} MeV, "
                          f"|<n|D_z|0>|={abs(t['d_me']):.4f} fm, "
                          f"contrib={t['contribution']:.4f} fm^3")

            residual = (alpha_E - alpha_E_exp) / alpha_E_exp * 100 if alpha_E > 0 else float('nan')
            print(f"    vs experiment: {residual:+.1f}%")

            results.append({
                'N_shells': N_shells,
                'hw': hw,
                'b_fm': b_fm,
                'E_gs': E_gs,
                'alpha_E_fm3': alpha_E,
                'ewsr_fraction': ewsr / ewsr_trk,
                'residual_pct': residual,
                'n_contributing': n_contributing,
            })

    # Summary
    print(f"\n{'='*72}")
    print(f"SUMMARY")
    print(f"{'='*72}")
    print(f"\n{'N':>2} {'hw':>6} {'b':>6} {'E_gs':>8} {'alpha_E':>8} {'EWSR%':>7} {'Residual':>10}")
    print(f"-" * 55)
    for r in results:
        print(f"{r['N_shells']:>2} {r['hw']:>6.1f} {r['b_fm']:>6.3f} "
              f"{r['E_gs']:>8.3f} {r['alpha_E_fm3']:>8.4f} "
              f"{r['ewsr_fraction']:>6.3f} {r['residual_pct']:>+9.1f}%")

    print(f"\nExperimental: alpha_E = {alpha_E_exp:.4f} fm^3")
    print(f"Chiral EFT:   alpha_E = {alpha_E_eft:.3f} fm^3")

    # Connection to muD Lamb shift
    # The polarizability contribution to muonic deuterium Lamb shift:
    # Delta_E_pol = (alpha_d / 3) * <|psi_mu(0)|^2> * integral_factor
    # CGV 2014: +1.690(20) meV → alpha_E ~ 0.63 fm^3
    # CV 2011: +1.94 meV → alpha_E ~ 0.72 fm^3
    # The ratio 1.94/1.69 = 1.148 would correspond to alpha_E = 0.727 fm^3
    # (15% larger than experimental)

    # Save
    with open('debug/data/deuteron_polarizability.json', 'w') as f:
        json.dump(results, f, indent=2, default=float)
    print(f"\nResults saved to debug/data/deuteron_polarizability.json")

    return results


if __name__ == '__main__':
    results = compute_polarizability()
