r"""Route C energy follow-on (v5.15.16 -> energy): the RADIAL LADDER that closes
the prolate all-electron LiH energy gap, built strictly ON TOP of the validated
C4 engine (debug/prolate_allelectron_c4.py).  Nothing in geovac/ is edited.

WHERE THIS STARTS.  C4 delivered the from-scratch prolate all-electron LiH
GEOMETRY at experiment (R_eq +0.2%, pi-converged) but the ENERGY sat ~60 mHa
above exact (-8.011 vs -8.070): every C4 orbital carries a SINGLE xi-power and a
SINGLE exponent, so the basis is radially incomplete.  C4's own note scoped the
cure as "a multi-exponent radial ladder (the prolate_recondition route)".

WHAT IS NEW HERE (no new integral machinery -- C4's build_Xtab_s already handles
arbitrary mixed exponents/weights; a radial ladder is just more basis functions):

1. ``build_channels(spec, alpha)`` -- a flexible orbital-set builder.  A ``spec``
   is a list of channels (mu, msign-list, J, L): each contributes the full grid
   of prolate primitives xi^j eta^l (weights) e^{-alpha xi} e^{i m phi} for
   j = 0..J (the RADIAL ladder) and l = 0..L (the angular ladder), over the given
   signed-m values.  mu = 0 (sigma), 1 (pi), 2 (delta), 3 (phi) are supported via
   C4's from_sigma / valence_pi_orbital / valence_delta_orbital and a local
   valence_phi_orbital (mu=3, C4's one_body_general/eri_general are general-mu).

2. ``build_C(tags, alpha)`` + ``assemble_rebased(...)`` -- the CONDITIONING FIX.
   A monomial radial ladder xi^j at fixed alpha is the classic ill-conditioned
   (Hankel) set; a plain float64 canonical solve of C4's assemble_energy_m drops
   functions and REGRESSES once the ladder is long (measured: delta J=4 gives
   cond(S) > 1e10, 66/75 kept, D_e% 99.51 -> 98.93).  This is exactly the wall
   geovac/prolate_recondition.py documents for Paper 12's H2.  The fix is the
   same: re-base the monomial primitives onto Laguerre(xi) x Legendre(eta) -- the
   SAME span, so the energy is unchanged, only the conditioning.  The change of
   basis C is block-diagonal in (mu, signed m) and factorizes radial (x) angular,
   built here from prolate_recondition.laguerre_coeffs / legendre_coeffs.  Because
   Laguerre is a WELL-conditioned target (unlike the S^{-1/2} orthogonalizer), the
   re-basing is accurate in float64: validated to reproduce the monomial energy
   bit-for-bit at d(2,2) (E = -1.173738, D_e% 99.578) while cutting cond(S)
   2.9e9 -> 1.4e5, so the long ladders no longer collapse.

Everything else (the ERIs eri_general, one-body one_body_general, the FCI
fci_energy) is imported from C4 unchanged.

Run from root:
  python debug/prolate_energy_ladder.py reduce    # G-REDUCE (ladder=1 == C4)
  python debug/prolate_energy_ladder.py h2         # G-H2RAD (the make-or-break)
  python debug/prolate_energy_ladder.py lih        # LiH energy convergence @ R=3.015
  python debug/prolate_energy_ladder.py all        # everything -> data/prolate_energy_ladder.log
"""
from __future__ import annotations

import os
import sys
import time
from collections import defaultdict
from typing import Dict, List, Sequence, Tuple

import numpy as np
import mpmath as mp

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))                    # debug/
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))   # root

from geovac import prolate_recondition as pr                       # noqa: E402
import prolate_allelectron_c4 as c4                                 # noqa: E402
from prolate_allelectron_c4 import (                                # noqa: E402
    OrbitalM, from_sigma, valence_pi_orbital, valence_delta_orbital,
    one_body_general, build_eri_tensor_m, eri_general, eri_sigma, assemble_energy_m,
    ZC_LI, _canon4,
)
from prolate_mixed_eri import valence_prolate_orbital, sto_orbital  # noqa: E402
from prolate_allelectron_analytic_fci import sto_orbital_B          # noqa: E402
from prolate_allelectron_fci import fci_energy                      # noqa: E402

mp.mp.dps = 60
DATA = os.path.join(os.path.dirname(os.path.abspath(__file__)), "data")

R_H2 = 1.40
E_EXACT_H2 = -1.174475
DE_H2 = 0.174475
R_LIH = 3.015
E_EXACT_LIH = -8.070


# ==========================================================================
# mu=3 (phi) valence orbital -- C4 supplies mu=0,1,2; the engine is general-mu.
# ==========================================================================
def valence_phi_orbital(j: int, l: int, alpha, msign: int) -> OrbitalM:
    """phi valence (mu=3), msign in {+3,-3}."""
    assert abs(msign) == 3
    return OrbitalM(j, c4._shift([mp.mpf(1)], l), alpha, mp.mpf(1), mu=3, msign=msign)


_MSIGNS = {0: [0], 1: [+1, -1], 2: [+2, -2], 3: [+3, -3]}


def _make_orb(mu: int, j: int, l: int, alpha, msign: int) -> OrbitalM:
    if mu == 0:
        return from_sigma(valence_prolate_orbital(j, l, mp.mpf(alpha)))
    if mu == 1:
        return valence_pi_orbital(j, l, mp.mpf(alpha), msign)
    if mu == 2:
        return valence_delta_orbital(j, l, mp.mpf(alpha), msign)
    if mu == 3:
        return valence_phi_orbital(j, l, mp.mpf(alpha), msign)
    raise ValueError(f"mu={mu} unsupported")


# ==========================================================================
# Orbital-set builder: a channel is (mu, J, L); every channel is a full
# radial (j=0..J) x angular (l=0..L) prolate grid, over its signed-m values.
# ==========================================================================
def sig(J: int, L: int) -> Tuple[int, int, int]:
    return (0, J, L)


def pi(J: int, L: int) -> Tuple[int, int, int]:
    return (1, J, L)


def de(J: int, L: int) -> Tuple[int, int, int]:
    return (2, J, L)


def ph(J: int, L: int) -> Tuple[int, int, int]:
    return (3, J, L)


def build_channels(spec: Sequence[Tuple[int, int, int]], alpha: float
                   ) -> Tuple[List[OrbitalM], List[Tuple]]:
    """spec = list of (mu, J, L).  Returns (orbitals, tags) with one tag
    ((mu, msign, J, L), j, l) per orbital -- tags drive the re-basing."""
    orbs: List[OrbitalM] = []
    tags: List[Tuple] = []
    for (mu, J, L) in spec:
        for ms in _MSIGNS[mu]:
            for j in range(J + 1):
                for l in range(L + 1):
                    orbs.append(_make_orb(mu, j, l, alpha, ms))
                    tags.append(((mu, ms, J, L), j, l))
    return orbs, tags


# ==========================================================================
# The conditioning fix: block-diagonal monomial -> Laguerre(xi) x Legendre(eta)
# ==========================================================================
def build_C(tags: Sequence[Tuple], alpha: float) -> np.ndarray:
    """float64 change of basis C[o, i] : primitive i=(j,l) -> orthogonal o=(nr,na)
    within each (mu, msign) channel.  C_block = Lag(nr,j) . Leg(na,l).  Square,
    block-diagonal, invertible (Laguerre/Legendre span the same polynomials)."""
    n = len(tags)
    chan: Dict[Tuple, List[Tuple[int, int, int]]] = defaultdict(list)
    for i, (key, j, l) in enumerate(tags):
        chan[key].append((i, j, l))
    C = np.zeros((n, n))
    with mp.workdps(60):
        for key, items in chan.items():
            # width from the items themselves (robust to any key structure); a
            # single-power channel (Jmax=Lmax=0) re-bases to the 1x1 identity.
            Jmax = max(j for (_, j, _) in items)
            Lmax = max(l for (_, _, l) in items)
            lag = [[float(x) for x in pr.laguerre_coeffs(nr, alpha, Jmax + 1)]
                   for nr in range(Jmax + 1)]
            leg = [[float(x) for x in pr.legendre_coeffs(na, Lmax + 1)]
                   for na in range(Lmax + 1)]
            # orthogonal fn o=(nr,na) placed on the row of the primitive with (j=nr,l=na)
            for (i, jr, lr) in items:            # this row carries orthogonal fn (nr=jr, na=lr)
                nr, na = jr, lr
                for (ip, jp, lp) in items:        # primitive column (jp, lp)
                    C[i, ip] = lag[nr][jp] * leg[na][lp]
    return C


def assemble_rebased(orbs: Sequence[OrbitalM], tags: Sequence[Tuple], R: float,
                     Z_A: float, Z_B: float, nelec: int, Vnn: float,
                     alpha: float, tol: float = 1e-12, verbose: bool = False
                     ) -> Tuple[float, float, int, int, float]:
    """Robust re-based float64 FCI.  Build primitive S/h1/ERI (exact mpf -> f64),
    re-base to Laguerre x Legendre (well-conditioned), canonical-orthogonalize,
    fci_energy.  Reduces to C4's assemble_energy_m for a single-power ladder."""
    t0 = time.time()
    S, h1 = one_body_general(orbs, R, Z_A, Z_B)
    eri = build_eri_tensor_m(orbs, R, verbose=verbose)
    C = build_C(tags, alpha)
    So = C @ S @ C.T
    h1o = C @ h1 @ C.T
    erio = np.einsum('pa,qb,rc,sd,abcd->pqrs', C, C, C, C, eri, optimize=True)
    w, U = np.linalg.eigh(So)
    wmax = w[-1]
    keep = w > tol * wmax
    X = U[:, keep] / np.sqrt(w[keep])
    Mk = int(keep.sum())
    h1_f = X.T @ h1o @ X
    eri_f = np.einsum('ap,bq,cr,ds,abcd->pqrs', X, X, X, X, erio, optimize=True)
    E_elec, ndet = fci_energy(h1_f, eri_f, Mk, nelec)
    E_tot = E_elec + Vnn
    cond = wmax / max(w[keep].min(), 1e-300)
    if verbose:
        print(f"    M={len(orbs)} kept={Mk} ndet={ndet} cond(S_o)={cond:.1e} "
              f"E_tot={E_tot:.6f}  [{time.time()-t0:.0f}s]", flush=True)
    return E_tot, E_elec, Mk, ndet, cond


def _fmt_spec(spec: Sequence[Tuple[int, int, int]]) -> str:
    names = {0: "s", 1: "pi", 2: "d", 3: "phi"}
    return "+".join(f"{names[mu]}({J},{L})" for (mu, J, L) in spec)


# ==========================================================================
# G-REDUCE:  ladder size 1 (single power per channel) == C4 numbers
# ==========================================================================
def gate_reduce() -> bool:
    """The re-based assemble on a single-power ladder must reproduce C4's
    assemble_energy_m bit-for-bit (change of basis on a 1x1 radial block is a
    scalar).  Uses the C4 H2 sigma+pi set."""
    print("=" * 72)
    print("G-REDUCE  re-based assemble on a size-1 ladder == C4 assemble_energy_m")
    print("=" * 72)
    # a small sigma+pi set at J=L=1 per channel, single exponent
    spec = [sig(1, 1), pi(1, 1)]
    alpha = 1.0
    orbs, tags = build_channels(spec, alpha)
    Et_r, _, Mk_r, _, _ = assemble_rebased(orbs, tags, R_H2, 1.0, 1.0, 2, 1.0 / R_H2, alpha)
    Et_c, _, Mk_c, _, _ = assemble_energy_m(orbs, R_H2, 1.0, 1.0, 2, 1.0 / R_H2)
    d = abs(Et_r - Et_c)
    ok = d < 1e-9
    print(f"  {_fmt_spec(spec)}:  rebased E={Et_r:.8f} (M={Mk_r})   "
          f"C4 E={Et_c:.8f} (M={Mk_c})   |dE|={d:.1e}")
    # and a size-1 ladder is literally the C4 basis -> also check the C4 H2 builder
    orbs2 = c4.h2_orbitals_c4(M_sigma=6, n_pi=1, alpha=1.0, alpha_pi=1.0)
    Et2, _, _, _, _ = assemble_energy_m(orbs2, R_H2, 1.0, 1.0, 2, 1.0 / R_H2)
    print(f"  [C4 h2 sigma6+1pi baseline E={Et2:.6f}  D_e%={100*(-1-Et2)/DE_H2:.2f}]")
    print(f"  G-REDUCE: {'PASS' if ok else 'FAIL'}")
    return ok


# ==========================================================================
# G-H2RAD (make-or-break): H2 radial+azimuthal ladder convergence
# ==========================================================================
def _de_pct(Et: float) -> float:
    return 100.0 * (-1.0 - Et) / DE_H2


def run_h2(spec, alpha, tol=1e-12, tag=None):
    orbs, tags = build_channels(spec, alpha)
    Et, Ee, Mk, nd, cond = assemble_rebased(orbs, tags, R_H2, 1.0, 1.0, 2,
                                             1.0 / R_H2, alpha)
    dep = _de_pct(Et)
    err = (E_EXACT_H2 - Et) * 1000.0
    name = tag or _fmt_spec(spec)
    print(f"  {name:28s} a={alpha:.2f} M={len(orbs):3d} kept={Mk:3d} nd={nd:5d} "
          f"cond={cond:.1e} E={Et:.6f} D_e%={dep:6.3f} err={err:+.3f}mHa", flush=True)
    return dep, Et, len(orbs)


def gate_h2(configs=None):
    """The convergence table.  Radial ladder saturates fast; the azimuthal (pi,
    delta, phi) channels are the lever; the re-basing keeps every function."""
    print("=" * 72)
    print(f"G-H2RAD  H2 radial+azimuthal ladder  (R={R_H2}, exact {E_EXACT_H2}, "
          f"D_e {DE_H2})")
    print("=" * 72)
    if configs is None:
        configs = [
            ("sigma only (radial saturates)", [sig(4, 2)], 1.2),
            ("+pi", [sig(4, 2), pi(4, 2)], 1.2),
            ("+delta", [sig(4, 2), pi(4, 2), de(3, 2)], 1.2),
            ("+delta radial (rebased)", [sig(4, 2), pi(4, 2), de(4, 2)], 1.2),
            ("+phi", [sig(4, 2), pi(4, 2), de(3, 2), ph(1, 1)], 1.2),
            ("deep", [sig(4, 2), pi(5, 2), de(4, 2), ph(2, 1)], 1.2),
        ]
    best = 0.0
    for tag, spec, alpha in configs:
        dep, Et, M = run_h2(spec, alpha, tag=tag)
        best = max(best, dep)
    print(f"  best D_e% = {best:.3f}  "
          f"({'PASS >=99.8' if best >= 99.8 else 'chemical acc' if best >= 99.0 else 'LOW'})")
    return best


# ==========================================================================
# LiH energy convergence at fixed geometry (R = 3.015)
# ==========================================================================
def lih_core_orbitals(R) -> Tuple[List[OrbitalM], List[Tuple]]:
    """The Route-A analytic core + minimal separated-atom sigma anchors.
    Tagged mu=0 single-power (no radial re-basing within the core anchors --
    they are distinct physical exponents, kept as their own 1x1 blocks)."""
    orbs, tags = [], []

    def add(o, key):
        orbs.append(o)
        tags.append((key, 0, 0))
    add(from_sigma(sto_orbital(ZC_LI, R, is_core=True)), ('coreLi', 0, 0))
    add(from_sigma(sto_orbital(mp.mpf('0.65'), R)), ('Li2s', 0, 0))
    add(from_sigma(sto_orbital_B(mp.mpf('1.0'), R)), ('H1s', 0, 0))
    add(from_sigma(sto_orbital_B(mp.mpf('0.70'), R)), ('Hm', 0, 0))
    return orbs, tags


def lih_orbitals_ladder(R, Jbond=2, Lbond=1, npi=1, Jpi=1, Lpi=0, alpha=1.0):
    """LiH set: analytic core + separated-atom STO anchors + a bond-centred
    sigma radial ladder (j=0..Jbond, l=0..Lbond) + pi shells.  The bond ladder
    is one re-basable channel (single alpha); the atom-centred anchors are their
    own 1x1 blocks."""
    orbs, tags = lih_core_orbitals(R)
    # bond-centred sigma radial+angular ladder (re-basable channel, key includes exps)
    for j in range(Jbond + 1):
        for l in range(Lbond + 1):
            orbs.append(from_sigma(valence_prolate_orbital(j, l, mp.mpf(alpha))))
            tags.append((('bond', 0, Jbond, Lbond), j, l))
    # pi shells (bond-centred), re-basable channel
    if npi > 0:
        for ms in (+1, -1):
            for j in range(Jpi + 1):
                for l in range(Lpi + 1):
                    orbs.append(valence_pi_orbital(j, l, mp.mpf(alpha), ms))
                    tags.append((('bondpi', ms, Jpi, Lpi), j, l))
    return orbs, tags


def _build_C_lih(tags, alpha):
    """LiH re-basing: the atom-centred anchors ('coreLi', etc.) are 1x1 identity
    blocks; the bond / bondpi channels re-base to Laguerre x Legendre."""
    return build_C(tags, alpha)


def run_lih_point(R, Jbond, Lbond, npi, Jpi, Lpi, alpha, tol=1e-11, verbose=True):
    orbs, tags = lih_orbitals_ladder(R, Jbond, Lbond, npi, Jpi, Lpi, alpha)
    Et, Ee, Mk, nd, cond = assemble_rebased(orbs, tags, R, 3.0, 1.0, 4, 3.0 / R,
                                            alpha, tol=tol)
    err = (E_EXACT_LIH - Et) * 1000.0
    if verbose:
        print(f"  bond({Jbond},{Lbond})+{npi}pi({Jpi},{Lpi}) a={alpha:.2f}  "
              f"M={len(orbs):2d} kept={Mk:2d} nd={nd:5d} cond={cond:.1e}  "
              f"E={Et:.5f}  err={err:+.1f}mHa", flush=True)
    return Et, len(orbs), nd


def scan_lih_energy(R=R_LIH, det_cap=17000):
    """LiH energy convergence at fixed R, growing the bond radial ladder + pi.

    HARD CONSTRAINT: LiH is a 4-electron FCI, so ndet = C(M,2)^2 (na=nb=2).  A
    dense ndet x ndet Hamiltonian caps the basis at M ~ 16 (ndet ~ 14400); this
    is the binding limit for LiH, not conditioning or ERI cost.  Configs past the
    cap are skipped and reported.  Each config = (Jbond,Lbond, npi,Jpi,Lpi, alpha)."""
    print("=" * 72)
    print(f"LiH energy convergence at R={R}  (exact {E_EXACT_LIH}; C4 sigma+pi -8.011)")
    print(f"  4e determinant wall: ndet=C(M,2)^2; cap M~16 (ndet~14400)")
    print("=" * 72)
    configs = [
        (1, 0, 0, 0, 0, 1.0),     # core+anchors + minimal bond, no pi
        (2, 1, 0, 0, 0, 1.0),     # bond sigma radial+angular ladder
        (2, 1, 1, 1, 0, 1.0),     # + 1 pi shell (~ C4 sigma+pi region)
        (3, 1, 1, 1, 0, 1.0),     # deeper bond radial ladder
        (2, 1, 1, 2, 0, 1.0),     # deeper pi radial ladder
        (3, 1, 1, 2, 0, 1.0),     # both ladders (at/near the det wall)
    ]
    rows = []
    for (Jb, Lb, npi, Jpi, Lpi, a) in configs:
        orbs, _ = lih_orbitals_ladder(R, Jb, Lb, npi, Jpi, Lpi, a)
        M = len(orbs)
        nd_est = (M * (M - 1) // 2) ** 2
        if nd_est > det_cap:
            print(f"  bond({Jb},{Lb})+{npi}pi({Jpi},{Lpi})  M={M} ndet~{nd_est} "
                  f"> cap {det_cap}  SKIPPED (4e determinant wall)", flush=True)
            continue
        Et, M, nd = run_lih_point(R, Jb, Lb, npi, Jpi, Lpi, a)
        rows.append((Jb, Lb, npi, Jpi, Lpi, Et, M, nd))
    if rows:
        best = min(r[5] for r in rows)
        print(f"  best LiH E_tot = {best:.5f}  (err {(E_EXACT_LIH-best)*1000:+.1f} mHa "
              f"vs exact {E_EXACT_LIH})")
    return rows


# ==========================================================================
# drivers
# ==========================================================================
def run_all():
    os.makedirs(DATA, exist_ok=True)
    log_path = os.path.join(DATA, "prolate_energy_ladder.log")
    import io
    buf = io.StringIO()

    class Tee:
        def write(self, s):
            sys.__stdout__.write(s); buf.write(s)

        def flush(self):
            sys.__stdout__.flush()
    sys.stdout = Tee()
    try:
        print(f"Route C energy follow-on: radial ladder (mp.dps={mp.mp.dps})")
        print(f"date 2026-09-22\n")
        gate_reduce()
        print()
        gate_h2()
        print()
        scan_lih_energy()
    finally:
        with open(log_path, "w") as f:
            f.write(buf.getvalue())
        sys.stdout = sys.__stdout__
        print(f"\n[log written to {log_path}]")


if __name__ == "__main__":
    arg = sys.argv[1] if len(sys.argv) > 1 else "reduce"
    if arg == "reduce":
        gate_reduce()
    elif arg == "h2":
        gate_h2()
    elif arg == "lih":
        scan_lih_energy()
    elif arg == "all":
        run_all()
    else:
        gate_reduce()
