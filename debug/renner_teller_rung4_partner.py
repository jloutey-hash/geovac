"""m-seam rung 4: identify the sigma partner that drives a2.

For each R: (a) the sigma (m = 0, z-even) states ranked by their share of the
sigma relay; (b) where the dominant partner and pi live in space -- probability
within radius RHO of either proton vs within RHO of the centre; (c) its
angular-momentum content about the centre.  The non-crossing rule holds in the
z-even sigma block (three collinear centres: no separability, so no hidden
crossings), so the k-th level keeps its index from the united atom (R -> 0,
Z = 3) to the separated atoms.

Usage: python debug/renner_teller_rung4_partner.py
"""
import importlib.util
import json
import os

import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
_s = importlib.util.spec_from_file_location('g', os.path.join(ROOT, 'tests', 'test_renner_teller_rung4.py'))
g = importlib.util.module_from_spec(_s)
_s.loader.exec_module(g)

LMAX = 19
RHO = 0.75


def bounds_for(R):
    return sorted(set([0, 0.1, 0.25, 0.5, R * 0.5, R * 0.8, R, R + 0.4, R + 1, R + 2.5, R + 5,
                       R + 9, R + 15, R + 24, R + 36, R + 52, R + 75]))


def region_weights(r, labs, vec, R):
    """P(within RHO of a proton), P(within RHO of the centre); m = 0 or 1 block."""
    Nr = len(r)
    U, WU = g.U, g.WU
    m = labs[0][1]
    C = vec.reshape(len(labs), Nr)
    # radial amplitude u_l(r_i) = c / sqrt(W_i); probability weight at node i = W_i |psi|^2 r^2 -> |c|^2 form
    Y = np.array([g.real_Y(*lab)[:, 0] for lab in labs])  # phi = 0 slice
    # azimuthal integral of |Y|^2 over phi: m = 0 -> 2 pi * (Y at phi=0)^2; m = 1 'c' -> pi * (Y/cos 0)^2
    phi_fac = 2 * np.pi if m == 0 else np.pi
    amp = C.T @ Y  # (Nr, NU): sum_l c_{l,i} Y_l(u)  (DVR coefficient units)
    dens = phi_fac * amp ** 2 * WU[None, :]  # each node already carries its radial weight
    z = r[:, None] * U[None, :]
    rho_perp = r[:, None] * np.sqrt(1 - U[None, :] ** 2)
    dprot = np.minimum(np.hypot(z - R, rho_perp), np.hypot(z + R, rho_perp))
    near_p = float(np.sum(dens[dprot < RHO]))
    near_c = float(np.sum(dens[r[:, None] * np.ones_like(U)[None, :] < RHO]))
    lw = np.sum(C ** 2, axis=1)
    return near_p, near_c, {int(lab[0]): round(float(w), 3) for lab, w in zip(labs, lw) if w > 0.01}


rows = []
for R in (0.25, 0.5, 1.0, 1.5, 2.0, 3.0, 4.0, 6.0, 8.0, 12.0, 16.0):
    g.R = R
    r, T = g.fedvr(bounds_for(R), 10)
    kmax = 2 * LMAX
    pi_c = g.harmonics(1, 'c', LMAX, +1)
    sig = g.harmonics(0, 'c', LMAX, +1)
    ev, V = np.linalg.eigh(g.h0_block(r, T, pi_c, kmax))
    Epi, pi = ev[0], V[:, 0]
    es, Vs = np.linalg.eigh(g.h0_block(r, T, sig, kmax))
    cpl = Vs.T @ (g.radial_op(r, sig, pi_c, range(1, kmax + 1, 2), g.v1_fun) @ pi)
    contrib = -cpl ** 2 / (es - Epi)
    order = np.argsort(-np.abs(contrib))
    top = int(order[0])
    pp, pc, lw = region_weights(r, sig, Vs[:, top], R)
    qp, qc, lpi = region_weights(r, pi_c, pi, R)
    row = {'R': R, 'E_pi': float(Epi), 'sigma_levels': [float(x) for x in es[:4]],
           'partner_index': top + 1, 'partner_E': float(es[top]),
           'partner_share': float(contrib[top] / np.sum(contrib)),
           'sigma_relay': float(np.sum(contrib)),
           'second_index': int(order[1]) + 1, 'second_share': float(contrib[order[1]] / np.sum(contrib)),
           'partner_near_protons': pp, 'partner_near_centre': pc, 'partner_l': lw,
           'pi_near_protons': qp, 'pi_near_centre': qc}
    rows.append(row)
    print('R=%5.2f  E_pi=% .4f  sigma=%s  partner=#%d (E=% .4f, %.0f%% of relay %+.3f; next #%d %.0f%%)'
          '  near H %.2f / centre %.2f  [pi: H %.2f / centre %.2f]  l-content %s'
          % (R, Epi, np.round(es[:4], 4), top + 1, es[top], 100 * row['partner_share'], row['sigma_relay'],
             row['second_index'], 100 * row['second_share'], pp, pc, qp, qc, lw))
with open(os.path.join(ROOT, 'debug', 'data', 'renner_teller_rung4_partner.json'), 'w') as f:
    json.dump(rows, f, indent=1)
