"""
Paper 54 consolidated recompute: DIAGONAL (eq:A_full as written) vs DOUBLE sum
=============================================================================
Resolves which gauge-field construction the paper's published numbers came
from, so the paper can be made internally consistent.

  Construction DIAGONAL (matches the WRITTEN eq:A_full, single sum):
      A = Σ_i  (M_i† ⊗ M_i†) [D_total, M_i ⊗ M_i]
  Construction DOUBLE (matches the original Prop-3 driver
  tensor_product_radial_compare.py build_connected_interaction_full_algebra):
      A = Σ_i Σ_j  (M_i† ⊗ M_j†) [D_total, M_i ⊗ M_j]

For each construction and n_max in {2,3} we compute every number that appears
in Paper 54:
  - connected fraction  ||V_conn|| / ||{D,A}||           (Prop 2)
  - Gaunt %, m-conservation %, multipole-k decomposition (Theorem 2)
  - radial vs EXACT Slater two_electron_integral:
        Pearson, CV(ratio), sign-agreement                (Prop 3)
  - heat-trace ratio at t=0.5, eps in {0,0.1,0.5,1.0}      (Table I)

Published values to match:
  Prop 2:    76.7% (n=2), 32.4% (n=3)
  Theorem 2: Gaunt 100%/100%, m-cons 100%, pure k=0 both cutoffs
  Prop 3:    Pearson 0.58 (n=2), 0.41 (n=3), CV>1600%, sign 68%/54%
  Table I:   t=0.5 ratios 1.0 / 0.9955.. / 0.9778.. / 0.9568.. (n=2)
                          1.0 / 0.8384.. / 0.4102.. / 0.1854.. (n=3)
"""
import sys, json
import numpy as np

try:
    sys.stdout.reconfigure(encoding='utf-8', errors='replace')
except Exception:
    pass
sys.path.insert(0, '.')

from geovac.operator_system import build_multiplier_matrix, allowed_multiplier_labels
from geovac.casimir_ci import two_electron_integral
from debug.tensor_product_dirac import build_single_particle_scalar


def embedded_multipliers(sp, n_max):
    """3-Y multipliers embedded (block-diagonal) in the chirality-doubled space."""
    dim, N_sc, basis = sp['dim'], sp['N_scalar'], sp['scalar_basis']
    out = []
    for (N, L, M) in allowed_multiplier_labels(n_max):
        M_sc = build_multiplier_matrix(N, L, M, basis)
        if np.linalg.norm(M_sc) < 1e-15:
            continue
        M_full = np.zeros((dim, dim), dtype=complex)
        M_full[:N_sc, :N_sc] = M_sc
        M_full[N_sc:, N_sc:] = M_sc
        out.append(M_full)
    return out


def build_A(sp1, sp2, n_max, mode):
    d1, d2 = sp1['dim'], sp2['dim']
    D1, D2, g1 = sp1['D'], sp2['D'], sp1['gamma']
    D_total = np.kron(D1, np.eye(d2)) + np.kron(g1, D2)
    mults = embedded_multipliers(sp1, n_max)
    A = np.zeros((d1 * d2, d1 * d2), dtype=complex)
    if mode == 'diagonal':
        for Mi in mults:
            ci = D1 @ Mi - Mi @ D1
            if np.linalg.norm(ci) < 1e-14:
                continue
            fluct = np.kron(ci, Mi) + np.kron(g1 @ Mi, D2 @ Mi - Mi @ D2)
            A += np.kron(Mi.conj().T, Mi.conj().T) @ fluct
    elif mode == 'double':
        for Mi in mults:
            ci = D1 @ Mi - Mi @ D1
            if np.linalg.norm(ci) < 1e-14:
                continue
            for Mj in mults:
                cj = D2 @ Mj - Mj @ D2
                if np.linalg.norm(cj) < 1e-14:
                    continue
                fluct = np.kron(ci, Mj) + np.kron(g1 @ Mi, cj)
                A += np.kron(Mi.conj().T, Mj.conj().T) @ fluct
    A = ((A + A.conj().T) / 2).real
    return D_total, A


def connected_part(DA, d1, d2):
    I1, I2 = np.eye(d1), np.eye(d2)
    t = DA.reshape(d1, d2, d1, d2)
    Tr2 = np.trace(t, axis1=1, axis2=3)
    Tr1 = np.trace(t, axis1=0, axis2=2)
    fac = np.kron(Tr2 / d2, I2) + np.kron(I1, Tr1 / d1) - (np.trace(DA) / (d1 * d2)) * np.eye(d1 * d2)
    conn = DA - fac
    frac = np.linalg.norm(conn) / np.linalg.norm(DA) if np.linalg.norm(DA) > 0 else 0.0
    return float(frac), conn


def pp_tensor(conn, N, d1):
    V = np.zeros((N, N, N, N))
    for a in range(N):
        for b in range(N):
            for c in range(N):
                for d in range(N):
                    V[a, b, c, d] = conn[a * d1 + b, c * d1 + d]
    return V


def angular_audit(V, basis):
    g_ok = g_bad = m_ok = m_bad = 0.0
    by_k = {}
    N = len(basis)
    for a in range(N):
        for b in range(N):
            for c in range(N):
                for d in range(N):
                    v = V[a, b, c, d]
                    if abs(v) < 1e-12:
                        continue
                    dl1, dl2 = abs(basis[a].l - basis[c].l), abs(basis[b].l - basis[d].l)
                    if dl1 == dl2:
                        g_ok += v * v
                        by_k[dl1] = by_k.get(dl1, 0.0) + v * v
                    else:
                        g_bad += v * v
                    if basis[a].m + basis[b].m == basis[c].m + basis[d].m:
                        m_ok += v * v
                    else:
                        m_bad += v * v
    tg, tm, tk = g_ok + g_bad, m_ok + m_bad, sum(by_k.values())
    return {
        'gaunt_pct': 100 * g_ok / tg if tg else None,
        'mcons_pct': 100 * m_ok / tm if tm else None,
        'multipole_k': {int(k): round(100 * v / tk, 2) for k, v in sorted(by_k.items())} if tk else {},
    }


def coulomb_tensor(basis):
    N = len(basis)
    C = np.zeros((N, N, N, N))
    for a in range(N):
        for b in range(N):
            for c in range(N):
                for d in range(N):
                    val = two_electron_integral(
                        basis[a].n, basis[a].l, basis[a].m,
                        basis[b].n, basis[b].l, basis[b].m,
                        basis[c].n, basis[c].l, basis[c].m,
                        basis[d].n, basis[d].l, basis[d].m)
                    if abs(val) > 1e-15:
                        C[a, b, c, d] = val
    return C


def radial_vs_slater(V, C):
    s, c = V.ravel(), C.ravel()
    both = (np.abs(s) > 1e-12) & (np.abs(c) > 1e-12)
    n = int(both.sum())
    if n < 3:
        return {'pearson': None, 'cv': None, 'sign_agree': None, 'n_overlap': n,
                'n_spec': int((np.abs(s) > 1e-12).sum()), 'n_coul': int((np.abs(c) > 1e-12).sum())}
    ss, cc = s[both], c[both]
    ratio = ss / cc
    return {
        'pearson': float(np.corrcoef(ss, cc)[0, 1]),
        'cv': float(np.std(ratio) / abs(np.mean(ratio))) if np.mean(ratio) != 0 else None,
        'sign_agree': float(np.mean(np.sign(ss) == np.sign(cc))),
        'n_overlap': n, 'n_spec': int((np.abs(s) > 1e-12).sum()), 'n_coul': int((np.abs(c) > 1e-12).sum()),
    }


def heat_trace(D_total, A, D1, t=0.5, eps_list=(0.0, 0.1, 0.5, 1.0)):
    e1 = np.linalg.eigvalsh(D1)
    Z1 = float(np.sum(np.exp(-t * e1 ** 2)))
    denom = Z1 * Z1
    out = {}
    for eps in eps_list:
        lam = np.linalg.eigvalsh(D_total + eps * A)
        num = float(np.sum(np.exp(-t * lam ** 2)))
        out[eps] = num / denom
    return out


def main():
    results = {}
    for n_max in [2, 3]:
        sp1 = build_single_particle_scalar(n_max)
        sp2 = build_single_particle_scalar(n_max)
        N, d1 = sp1['N_scalar'], sp1['dim']
        basis = sp1['scalar_basis']
        C = coulomb_tensor(basis)
        results[n_max] = {}
        for mode in ['diagonal', 'double']:
            D_total, A = build_A(sp1, sp2, n_max, mode)
            DA = D_total @ A + A @ D_total
            frac, conn = connected_part(DA, d1, sp2['dim'])
            V = pp_tensor(conn, N, d1)
            ang = angular_audit(V, basis)
            rad = radial_vs_slater(V, C)
            heat = heat_trace(D_total, A, sp1['D'])
            results[n_max][mode] = {
                'connected_fraction_pct': round(100 * frac, 1),
                **ang, 'radial': rad,
                'heat_t0.5': {str(k): round(v, 6) for k, v in heat.items()},
            }
            print(f"n_max={n_max} [{mode:8s}] conn={100*frac:5.1f}%  "
                  f"Gaunt={ang['gaunt_pct']}  k={ang['multipole_k']}  "
                  f"Pearson={rad['pearson']}  CV={rad['cv']}  sign={rad['sign_agree']}")
            print(f"            heat t=0.5: " +
                  ", ".join(f"{k}:{v:.4f}" for k, v in heat.items()))
    with open('debug/data/paper54_recompute_both.json', 'w') as f:
        json.dump(results, f, indent=2)
    print("\nWrote debug/data/paper54_recompute_both.json")
    print("\nPUBLISHED: Prop2 76.7/32.4 | Thm2 Gaunt100/100 pure-k0 | "
          "Prop3 Pearson 0.58/0.41 CV>16 sign .68/.54 | "
          "TableI n3 1/0.8384/0.4102/0.1854")


if __name__ == '__main__':
    main()
