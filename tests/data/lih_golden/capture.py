r"""Golden-value capture for the LiH / prolate productionization (step 1 of
debug/lih_productionization_plan.md).

Calls the CURRENT debug/ entry points (and the already-promoted geovac/lih_r12ci) and stores
full-precision outputs as the oracle the promoted geovac packages must reproduce.

Every item runs in its OWN fresh Python process (module-level state -- mpmath dps set at import,
geovac.lih_r12ci.kernels.USE_EXACT_NEUMANN, the triangle / randomized-SVD RNG streams seeded at
import, ARPACK's start-vector state, lru caches -- is part of what a value depends on, so a fresh
process per item is what makes a value reproducible and attributable to one call sequence).

BLAS / OpenMP threads are pinned to 1 BEFORE numpy is imported (the production drivers do the
same; multi-threaded BLAS changes summation order and therefore the low bits).

Usage (from the repo root):
    python tests/data/lih_golden/capture.py --list
    python tests/data/lih_golden/capture.py --run-item NAME [--tag rep2] [--stage DIR]
    python tests/data/lih_golden/capture.py --compare NAME [--stage DIR]      # NAME.json vs NAME.rep2.json
    python tests/data/lih_golden/capture.py --merge [--stage DIR]             # -> golden.json + golden.npz

Serialization: floats as repr() strings (round-trip exact); arrays in golden.npz under
"<item>__<key>" with the sha256 of their C-contiguous bytes; very large arrays are recorded by
sha256 + summary statistics only ("stored": false).

Nothing under debug/ or geovac/ is modified: drivers that write artifacts are redirected (module
globals patched IN-PROCESS) to a temporary directory.
"""
from __future__ import annotations

import os

for _v in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS", "NUMEXPR_NUM_THREADS",
           "VECLIB_MAXIMUM_THREADS"):
    os.environ[_v] = "1"

import argparse
import datetime
import glob
import hashlib
import importlib.util
import json
import platform
import subprocess
import sys
import tempfile
import time
from typing import Any, Callable, Dict, List, Optional

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(os.path.dirname(HERE)))
DEBUG = os.path.join(ROOT, "debug")
DATA = os.path.join(DEBUG, "data")
TESTS = os.path.join(ROOT, "tests")
DEFAULT_STAGE = os.path.join(HERE, "_stage")
GOLDEN_JSON = os.path.join(HERE, "golden.json")
GOLDEN_NPZ = os.path.join(HERE, "golden.npz")
STORE_LIMIT_BYTES = 8_000_000          # arrays larger than this are recorded by sha256 + stats only


# =========================================================================== #
# registry
# =========================================================================== #
ITEMS: Dict[str, dict] = {}


def item(name: str, tier: str, call: str, seeds: Optional[dict] = None, notes: str = ""):
    """Register an item.  tier = 'fast' (<60 s expected) or 'slow'.  `call` is the exact
    debug/ (or geovac/) entry-point sequence, for the promotion agents."""
    def deco(fn: Callable[[], dict]):
        ITEMS[name] = dict(fn=fn, tier=tier, call=call, seeds=seeds or {}, notes=notes)
        return fn
    return deco


class ShaOnly:
    """Wrap an array to be recorded by sha256 + stats only (not stored in golden.npz)."""
    def __init__(self, a):
        self.a = np.ascontiguousarray(np.asarray(a))


def _paths() -> None:
    sys.argv = [sys.argv[0]]                 # debug drivers parse argv at import
    for p in (ROOT, DEBUG):
        if p not in sys.path:
            sys.path.insert(0, p)


def _load_module(path: str, name: str):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod
    spec.loader.exec_module(mod)
    return mod


def _sha(a: np.ndarray) -> str:
    return hashlib.sha256(np.ascontiguousarray(a).tobytes()).hexdigest()


def _stats(a: np.ndarray) -> dict:
    if a.size == 0 or not np.issubdtype(a.dtype, np.number):
        return {}
    b = a.astype(complex) if np.iscomplexobj(a) else a.astype(float)
    if np.iscomplexobj(b):
        return dict(sum=repr(complex(b.sum())), absmax=repr(float(np.abs(b).max())))
    return dict(sum=repr(float(b.sum())), min=repr(float(b.min())), max=repr(float(b.max())),
                l2=repr(float(np.sqrt((b * b).sum()))))


def _ser(key: str, v: Any, arrays: Dict[str, np.ndarray]) -> dict:
    if isinstance(v, ShaOnly):
        a = v.a
        return dict(type="array", shape=list(a.shape), dtype=str(a.dtype), sha256=_sha(a),
                    stored=False, stats=_stats(a))
    if isinstance(v, (bool, np.bool_)):
        return dict(type="bool", value=bool(v))
    if isinstance(v, (int, np.integer)):
        return dict(type="int", value=int(v))
    if isinstance(v, (float, np.floating)):
        return dict(type="float", repr=repr(float(v)))
    if isinstance(v, (complex, np.complexfloating)):
        return dict(type="complex", repr=repr(complex(v)))
    if isinstance(v, str):
        return dict(type="str", value=v)
    if v is None:
        return dict(type="none")
    if isinstance(v, dict):
        return dict(type="dict", value={str(k): _ser(f"{key}.{k}", x, arrays) for k, x in v.items()})
    if isinstance(v, (list, tuple)):
        try:
            a = np.asarray(v)
            if a.dtype != object and np.issubdtype(a.dtype, np.number):
                v = a
            else:
                return dict(type="list", value=[_ser(f"{key}[{i}]", x, arrays) for i, x in enumerate(v)])
        except Exception:
            return dict(type="list", value=[_ser(f"{key}[{i}]", x, arrays) for i, x in enumerate(v)])
    if isinstance(v, np.ndarray):
        a = np.ascontiguousarray(v)
        if a.dtype == object:
            return dict(type="list", value=[_ser(f"{key}[{i}]", x, arrays) for i, x in enumerate(a.tolist())])
        rec = dict(type="array", shape=list(a.shape), dtype=str(a.dtype), sha256=_sha(a), stats=_stats(a))
        if a.nbytes <= STORE_LIMIT_BYTES:
            arrays[key] = a
            rec["stored"] = True
        else:
            rec["stored"] = False
        return rec
    return dict(type="repr", value=repr(v))


def _peak_rss_gb() -> float:
    try:
        import psutil
        mi = psutil.Process().memory_info()
        return float(getattr(mi, "peak_wset", mi.rss)) / 1e9
    except Exception:
        return float("nan")


def _rel(a: float, b: float) -> float:
    return abs(a - b) / max(abs(b), 1e-300)


# =========================================================================== #
# A. prolate_allelectron_fci  (tests/test_prolate_allelectron_fci.py)
# =========================================================================== #
def _h2_orb(na: int, R: float = 1.4):
    from geovac.prolate_scf import get_orbital_on_grid
    o = get_orbital_on_grid(R=R, Z_A=1.0, Z_B=1.0, n_angular=na, m=0,
                            N_xi_solve=4000, N_xi_grid=36, N_eta_grid=36, xi_max_grid=12.0)
    o['R'] = R
    return o


@item("pafci_azimuthal_kernels", "fast",
      "geovac.prolate_scf.get_orbital_on_grid(R=1.4,Z_A=Z_B=1,n_angular=0|1,m=0,N_xi_solve=4000,N_xi_grid=36,"
      "N_eta_grid=36,xi_max_grid=12); prolate_allelectron_fci._azimuthal_kernels(o0,mu_max=2); "
      "prolate_allelectron_fci.vee_m(a,b,c,d,0,0,0,0,K) vs compute_vee_integral(a,c,b,d)",
      notes="test_prolate_allelectron_fci fast tests (mu=0 regression + toroidal kernel vs dphi quadrature)")
def _pafci_kernels() -> dict:
    _paths()
    from geovac.prolate_allelectron import fci_grid as A          # promoted (step 2a)
    from geovac.prolate_scf import compute_vee_integral
    R = 1.4
    o0, o1 = _h2_orb(0), _h2_orb(1)
    K = A._azimuthal_kernels(o0, mu_max=2)
    out: Dict[str, Any] = {}
    for tag, (a, b, c, d) in {"0000": (o0, o0, o0, o0), "0101": (o0, o1, o0, o1),
                              "0011": (o0, o0, o1, o1)}.items():
        out[f"vee_ref_{tag}"] = float(compute_vee_integral(a, c, b, d))
        out[f"vee_m_{tag}"] = float(A.vee_m(a, b, c, d, 0, 0, 0, 0, K))
    XI, ETA = np.meshgrid(o0['xi'], o0['eta'], indexing='ij')
    rho = (R / 2) * np.sqrt(np.maximum((XI ** 2 - 1) * (1 - ETA ** 2), 0.0)).ravel()
    z = (R / 2) * (XI * ETA).ravel()
    ng = rho.size
    order = np.argsort(rho)
    ia, ib = int(order[-1]), int(order[-25])
    aa = rho[ia] ** 2 + rho[ib] ** 2 + (z[ia] - z[ib]) ** 2
    bb = 2 * rho[ia] * rho[ib]
    dphi = np.linspace(0, 2 * np.pi, 16000, endpoint=False)
    dd = dphi[1] - dphi[0]
    for mu in (0, 1, 2):
        out[f"K{mu}_elem"] = float(K[mu].reshape(ng, ng)[ia, ib])
        out[f"K{mu}_numerical"] = float(2 * np.pi * np.sum(np.cos(mu * dphi) / np.sqrt(aa - bb * np.cos(dphi))) * dd)
        out[f"K{mu}"] = ShaOnly(np.asarray(K[mu]))
    out["o0_psi"] = np.asarray(o0['psi'])
    out["o1_psi"] = np.asarray(o1['psi'])
    out["ia_ib"] = np.array([ia, ib])
    return out


@item("pafci_hf_vs_eckart", "slow",
      "geovac.prolate_scf.eckart_scf_energy(1.4,1.0,N_xi_solve=5000,N_grid=44,xi_max_grid=13.0); "
      "prolate_allelectron_fci.build_mo_integrals(1.4,2,1.0,1.0,N_xi_solve=5000,N_grid=44,xi_max=13.0)",
      notes="test_prolate_allelectron_fci::test_hf_matches_eckart (slow)")
def _pafci_hf() -> dict:
    _paths()
    from geovac.prolate_allelectron import fci_grid as A          # promoted (step 2a)
    from geovac.prolate_scf import eckart_scf_energy
    R = 1.4
    ek = eckart_scf_energy(R, 1.0, N_xi_solve=5000, N_grid=44, xi_max_grid=13.0)
    h1, eri, M, cond = A.build_mo_integrals(R, 2, 1.0, 1.0, N_xi_solve=5000, N_grid=44, xi_max=13.0)
    e_det = 2 * h1[0, 0] + eri[0, 0, 0, 0]
    return dict(E_HF_eckart=float(ek['E_HF']), h1=np.asarray(h1), eri=np.asarray(eri), M=M,
                cond=float(cond), e_det_plus_Vnn=float(e_det + 1.0 / R))


# =========================================================================== #
# B. Route C energy engine  (tests/test_prolate_route_c_lih.py + Paper 19 banked numbers)
# =========================================================================== #
@item("routec_fci_fast_vs_dense_random", "fast",
      "rng=np.random.default_rng(1); random symmetric h1 (8x8), chemist eri*0.1 (8^4); "
      "prolate_allelectron_fci.fci_energy(h1,e,8,4) vs fci_fast.fci_energy_fast(h1,e,8,4)",
      seeds={"numpy_default_rng": 1},
      notes="test_prolate_route_c_lih::test_sparse_fci_bitexact_vs_dense. eigsh (ARPACK) start vector "
            "is internal -> energies may differ in the last bits run-to-run")
def _routec_fci() -> dict:
    _paths()
    from geovac.prolate_allelectron import fci_fast               # promoted (step 2a)
    from geovac.prolate_allelectron.fci_grid import fci_energy
    rng = np.random.default_rng(1)
    M = 8
    h1 = rng.standard_normal((M, M)); h1 = h1 + h1.T
    e = rng.standard_normal((M, M, M, M)) * 0.1
    e = e + e.transpose(1, 0, 3, 2)
    e = e + e.transpose(2, 3, 0, 1)
    Ed, nd = fci_energy(h1, e, M, 4)
    Ef, ndf = fci_fast.fci_energy_fast(h1, e, M, 4)
    return dict(E_dense=float(Ed), nd_dense=int(nd), E_fast=float(Ef), nd_fast=int(ndf))


@item("routec_helike_cusp", "fast",
      "lih_r12_ceiling_probe.energy(Z,zeta,kinds) for (2.0,27/16,RICH),(3.0,2.6875,RADIAL),(3.0,2.6875,RICH); "
      "RADIAL=['t2','s'], RICH=['u','t2','s','u2','ut2']",
      notes="test_prolate_route_c_lih He control + Li core cusp (the 30.1 mHa additive-F12 input)")
def _routec_helike() -> dict:
    _paths()
    from geovac.prolate_allelectron.helike_cusp import energy as helike_energy   # promoted (step 2a)
    RADIAL = ['t2', 's']
    RICH = ['u', 't2', 's', 'u2', 'ut2']
    he = helike_energy(2.0, 27 / 16, RICH)
    li_rad = helike_energy(3.0, 2.6875, RADIAL)
    li_rich = helike_energy(3.0, 2.6875, RICH)
    return dict(He_rich=[float(x) for x in he], Li_radial=[float(x) for x in li_rad],
                Li_rich=[float(x) for x in li_rich],
                Li_core_cusp_Ha=float(li_rad[1] - li_rich[1]))


def _fast_engine():
    """Patch the ladder to the productionized fast engine, exactly as the slow route-C test does."""
    from geovac.prolate_allelectron import energy_ladder as L     # promoted (step 2a)
    from geovac.prolate_allelectron import float_eri as F
    from geovac.prolate_allelectron import fci_fast
    L.build_eri_tensor_m = F.build_eri_tensor_m_f
    L.fci_energy = fci_fast.fci_energy_fast
    from geovac.prolate_allelectron import routec as P            # promoted (step 2a)
    return P


def _routec_run(Jb, Lb, npi, Jpi, Lpi, core2) -> dict:
    _paths()
    P = _fast_engine()
    import mpmath as mp
    Et, M, Mk, nd, cond, dt = P.run(Jb=Jb, Lb=Lb, npi=npi, Jpi=Jpi, Lpi=Lpi, alpha=1.0, core2=core2)
    return dict(E_tot=float(Et), M=int(M), Mk=int(Mk), ndet=int(nd), cond=float(cond),
                mp_dps=int(mp.mp.dps))


_FAST_CALL = ("prolate_energy_ladder.build_eri_tensor_m=prolate_float_eri.build_eri_tensor_m_f; "
              "prolate_energy_ladder.fci_energy=fci_fast.fci_energy_fast; lih_core2exp_probe.run(")


@item("routec_float64_M6", "slow", _FAST_CALL + "Jb=1,Lb=0,npi=0,Jpi=0,Lpi=0,alpha=1.0,core2=None)",
      notes="test_prolate_route_c_lih::test_float64_engine_reproduces_mpf_lih (banked mpf -7.99468). "
            "eigsh start vector internal")
def _routec_m6() -> dict:
    return _routec_run(1, 0, 0, 0, 0, None)


@item("routec_M12_sigma_ref", "slow", _FAST_CALL + "Jb=2,Lb=1,npi=0,Jpi=0,Lpi=0,alpha=1.0,core2=[4.5,1.6])",
      notes="the sigma core-enriched reference of the sigma marriage, Paper 19 E=-8.02335")
def _routec_m12() -> dict:
    return _routec_run(2, 1, 0, 0, 0, [4.5, 1.6])


@item("routec_M16_ceiling", "slow", _FAST_CALL + "Jb=2,Lb=1,npi=1,Jpi=1,Lpi=0,alpha=1.0,core2=[4.5,1.6])",
      notes="M=16 orbital ceiling -8.02905 (Paper 19; the pi marriage reference basis)")
def _routec_m16() -> dict:
    return _routec_run(2, 1, 1, 1, 0, [4.5, 1.6])


@item("routec_M20_push", "slow", _FAST_CALL + "Jb=3,Lb=1,npi=1,Jpi=2,Lpi=0,alpha=1.0,core2=[4.5,1.6])",
      notes="pure-orbital ceiling -8.03170 (prolate_float_eri.push '+both'; Paper 19 '-8.032')")
def _routec_m20() -> dict:
    return _routec_run(3, 1, 1, 2, 0, [4.5, 1.6])


C4_RS = ["2.70", "2.85", "3.015", "3.20", "3.45"]


def _make_c4_item(Rs: str):
    @item(f"routec_c4_req_R{Rs}", "slow",
          f"prolate_allelectron_c4.assemble_energy_m(lih_orbitals_c4({Rs},extended=True,n_pi=1,alpha_pi=1.1),"
          f"{Rs},3.0,1.0,4,3.0/{Rs})  [mpf dps=60 analytic engine]",
          notes="Route C R_eq curve (ext+1pi(a=1.1)); R_eq = parabola through the 3 points around the min "
                "(prolate_allelectron_c4._fit_req) -> 3.020 bohr, +0.2% (debug/data/lih_analytic_c4.log)")
    def _fn() -> dict:
        _paths()
        import mpmath as mp
        from geovac.prolate_allelectron import c4                 # promoted (step 2a)
        R = float(Rs)
        Et, Ee, Mk, nd, cond = c4.assemble_energy_m(
            c4.lih_orbitals_c4(R, extended=True, n_pi=1, alpha_pi=1.1), R, 3.0, 1.0, 4, 3.0 / R)
        return dict(R=R, E_tot=float(Et), E_elec=float(Ee), Mk=int(Mk), ndet=int(nd), cond=float(cond),
                    mp_dps=int(mp.mp.dps))
    return _fn


for _Rs in C4_RS:
    _make_c4_item(_Rs)


# =========================================================================== #
# C. VMC  (tests/test_lih_vmc.py + short fixed-seed runs)
# =========================================================================== #
def _vmc_orbs():
    import mpmath as mp
    from geovac.prolate_allelectron.mixed_eri import sto_orbital, valence_prolate_orbital, ZC_LI
    from geovac.prolate_allelectron.analytic_fci import sto_orbital_B
    from geovac.prolate_allelectron.c4 import from_sigma, valence_pi_orbital
    R = 3.015
    return R, {
        "Li1s_core": from_sigma(sto_orbital(ZC_LI, R, is_core=True)),
        "core2_4.5": from_sigma(sto_orbital(mp.mpf('4.5'), R, is_core=True)),
        "H1s": from_sigma(sto_orbital_B(mp.mpf('1.0'), R)),
        "bond21": from_sigma(valence_prolate_orbital(2, 1, mp.mpf('1.0'))),
        "pi_p00": valence_pi_orbital(0, 0, mp.mpf('1.0'), +1),
        "pi_m10": valence_pi_orbital(1, 0, mp.mpf('1.0'), -1),
    }


@item("vmc_orbitals", "fast",
      "lih_vmc.extract_prim(orb); lih_vmc.orbital_value_grad(prim,pts,R) [pts: default_rng(1)]; "
      "lih_vmc.orbital_vgl(prim,pts,R) [pts: default_rng(3)]; prolate_allelectron_c4._orb_on_grid(orb,"
      "_grid_template(R,N_grid=24,xi_max=10.0)) at the middle node",
      seeds={"numpy_default_rng_grad": 1, "numpy_default_rng_lap": 3},
      notes="test_lih_vmc orbital gates; orbitals = Li1s_core, core2_4.5, H1s, bond21, pi_p00, pi_m10")
def _vmc_orbitals() -> dict:
    _paths()
    import geovac.lih_vmc as V                                    # promoted (step 2b)
    from geovac.prolate_allelectron.c4 import _orb_on_grid, _grid_template
    R, orbs = _vmc_orbs()
    rng = np.random.default_rng(1)
    pts1 = rng.uniform(-2, 2, size=(6, 3)); pts1[:, 2] += rng.uniform(-1, 1, 6)
    rng = np.random.default_rng(3)
    pts3 = rng.uniform(-2, 2, size=(6, 3)); pts3[:, 2] += rng.uniform(-1, 1, 6)
    tmpl = _grid_template(R, N_grid=24, xi_max=10.0)
    xi, eta = tmpl['xi'], tmpl['eta']
    ii, jj = len(xi) // 2, len(eta) // 2
    xi0, eta0 = xi[ii], eta[jj]
    rc = np.array([[(R / 2) * np.sqrt(max((xi0 ** 2 - 1) * (1 - eta0 ** 2), 0.0)), 0.0, (R / 2) * xi0 * eta0]])
    out: Dict[str, Any] = dict(pts_grad=pts1, pts_lap=pts3)
    for name, orb in orbs.items():
        prim = V.extract_prim(orb)
        out[f"{name}.prim"] = {k: (np.asarray(v, dtype=float) if isinstance(v, (list, np.ndarray)) else v)
                               for k, v in prim.items()}
        v1, g1 = V.orbital_value_grad(prim, pts1, R)
        v3, g3, l3 = V.orbital_vgl(prim, pts3, R)
        out[f"{name}.val"] = np.asarray(v1); out[f"{name}.grad"] = np.asarray(g1)
        out[f"{name}.vgl_val"] = np.asarray(v3); out[f"{name}.vgl_grad"] = np.asarray(g3)
        out[f"{name}.vgl_lap"] = np.asarray(l3)
        val_rc, _ = V.orbital_value_grad(prim, rc, R)
        out[f"{name}.realspace_at_node"] = complex(val_rc[0]) if np.iscomplexobj(val_rc) else float(val_rc[0])
        g = _orb_on_grid(orb, tmpl)['psi']
        out[f"{name}.engine_grid_at_node"] = complex(g[ii, jj]) if np.iscomplexobj(g) else float(g[ii, jj])
        out[f"{name}.engine_grid_psi"] = np.asarray(g)
    return out


@item("vmc_jastrow_cusp", "fast",
      "synthetic 1-det LiHWavefunction (Li1s_core, H1s; T=I, CS=[[1]]); PadeJastrow(1.5,2,2); "
      "_init_walkers(wf,1,default_rng(0)); local_energy_analytic at r_02 in (0.5,0.1,0.02,0.005)",
      seeds={"numpy_default_rng": 0}, notes="test_lih_vmc::test_jastrow_cusp_finite_local_energy")
def _vmc_cusp() -> dict:
    _paths()
    import mpmath as mp
    import geovac.lih_vmc as V                                    # promoted (step 2b)
    from geovac.prolate_allelectron.mixed_eri import sto_orbital, ZC_LI
    from geovac.prolate_allelectron.analytic_fci import sto_orbital_B
    from geovac.prolate_allelectron.c4 import from_sigma
    R = 3.015
    orbs = [from_sigma(sto_orbital(ZC_LI, R, is_core=True)), from_sigma(sto_orbital_B(mp.mpf('1.0'), R))]
    prims = [V.extract_prim(o) for o in orbs]
    wf = V.LiHWavefunction(prims, np.eye(2), np.array([1.0]), [(0, 1)], [(0, 1)], np.array([[1.0]]),
                           R, 3.0, 1.0, 3.0 / R, -8.0, 2, 2, 2)
    jas = V.PadeJastrow(1.5, 2, 2)
    rng = np.random.default_rng(0)
    r = V._init_walkers(wf, 1, rng)
    out: Dict[str, Any] = {"r_init": r.copy()}
    for dr in (0.5, 0.1, 0.02, 0.005):
        r[:, 2] = r[:, 0] + np.array([dr, 0.0, 0.0])
        EL, extra = V.local_energy_analytic(wf, r, jas)
        out[f"EL_dr{dr}"] = float(np.asarray(EL)[0])
    return out


@item("vmc_cached_small_short", "fast",
      "lih_vmc.cached_wavefunction('small') [reads debug/data/lih_vmc_wf_small.pkl]; "
      "lih_vmc.vmc(wf,nwalk=200,nsweep=300,nburn=100,jas=None,seed=0); "
      "lih_vmc.vmc(wf,...,jas=PadeJastrow(1.5,2,2),seed=0); "
      "lih_vmc.vmc_linear(wf,LinearGeminalJastrow('linexp',0.5,0.2,3.5,2,2),nwalk=200,nsweep=300,nburn=100,seed=0)",
      seeds={"vmc_seed": 0},
      notes="gate-6 machinery (VMC(J=1) ~ FCI) at short, exactly reproducible fixed-seed statistics")
def _vmc_small_short() -> dict:
    _paths()
    import geovac.lih_vmc as V                                    # promoted (step 2b)
    wf = V.cached_wavefunction("small")
    out: Dict[str, Any] = dict(E_fci_tot=float(wf.E_fci + wf.Vnn), Mk=int(wf.Mk), nd=int(len(wf.c)))
    E, err, acc = V.vmc(wf, nwalk=200, nsweep=300, nburn=100, jas=None, seed=0, verbose=False)
    out.update(gate6_E=float(E), gate6_err=float(err), gate6_acc=float(acc))
    E, err, acc = V.vmc(wf, nwalk=200, nsweep=300, nburn=100, jas=V.PadeJastrow(1.5, 2, 2), seed=0,
                        verbose=False)
    out.update(pade_E=float(E), pade_err=float(err), pade_acc=float(acc))
    lin = V.LinearGeminalJastrow('linexp', 0.5, 0.2, 3.5, 2, 2)
    E, err, acc = V.vmc_linear(wf, lin, nwalk=200, nsweep=300, nburn=100, seed=0, verbose=False)
    out.update(linear_E=float(E), linear_err=float(err), linear_acc=float(acc))
    return out


@item("vmc_build_small_wf", "slow",
      "lih_vmc.build_lih_wavefunction(Jb=1,Lb=0,npi=0,core2=None,verbose=False) [from scratch, float64 ERI + "
      "sparse FCI eigenvector] vs the cached debug/data/lih_vmc_wf_small.pkl",
      notes="eigsh eigenvector: sign + last bits are not fixed by the call (no v0); compare |c| / up to sign")
def _vmc_build_small() -> dict:
    _paths()
    import geovac.lih_vmc as V                                    # promoted (step 2b)
    wf = V.build_lih_wavefunction(Jb=1, Lb=0, npi=0, core2=None, verbose=False)
    wc = V.cached_wavefunction("small")
    s = np.sign(np.dot(wf.c, wc.c))
    return dict(E_fci_elec=float(wf.E_fci), E_fci_tot=float(wf.E_fci + wf.Vnn), Mk=int(wf.Mk),
                c=np.asarray(wf.c), T=np.asarray(wf.T), CS=np.asarray(wf.CS),
                cached_E_fci_elec=float(wc.E_fci),
                max_abs_dc_vs_cache_up_to_sign=float(np.max(np.abs(s * wf.c - wc.c))),
                max_abs_dT_vs_cache=float(np.max(np.abs(np.asarray(wf.T) - np.asarray(wc.T)))))


@item("vmc_build_m16_wf", "slow",
      "lih_vmc.build_lih_wavefunction(Jb=2,Lb=1,npi=1,Jpi=1,Lpi=0,core2=(4.5,1.6),verbose=False) vs "
      "debug/data/lih_vmc_wf_m16.pkl",
      notes="the M=16 VMC-over-FCI wavefunction (-8.029 basis); eigsh sign/last bits not fixed")
def _vmc_build_m16() -> dict:
    _paths()
    import geovac.lih_vmc as V                                    # promoted (step 2b)
    wf = V.build_lih_wavefunction(Jb=2, Lb=1, npi=1, Jpi=1, Lpi=0, core2=(4.5, 1.6), verbose=False)
    wc = V.cached_wavefunction("m16")
    s = np.sign(np.real(np.vdot(wc.c, wf.c)))
    return dict(E_fci_elec=float(wf.E_fci), E_fci_tot=float(wf.E_fci + wf.Vnn), Mk=int(wf.Mk),
                nd=int(len(wf.c)), c=np.asarray(wf.c), T=np.asarray(wf.T),
                cached_E_fci_elec=float(wc.E_fci),
                max_abs_dc_vs_cache_up_to_sign=float(np.max(np.abs(s * wf.c - wc.c))),
                max_abs_dT_vs_cache=float(np.max(np.abs(np.asarray(wf.T) - np.asarray(wc.T)))))


@item("vmc_cached_m16_short", "slow",
      "lih_vmc.cached_wavefunction('m16') [debug/data/lih_vmc_wf_m16.pkl]; "
      "lih_vmc.vmc(wf,nwalk=100,nsweep=120,nburn=20,jas=None,seed=10); "
      "lih_vmc.vmc(wf,nwalk=100,nsweep=120,nburn=20,jas=PadeJastrow(1.5,2,2),seed=10)",
      seeds={"vmc_seed": 10},
      notes="the -8.047 VMC-over-FCI machinery (Pade gamma=1.5) at short, exactly reproducible statistics; "
            "the banked -8.04731 +/- 0.00073 is final_run(nwalk=4000,nsweep=8000,nburn=1500,seeds 10,11,12) "
            "-- NOT captured (multi-hour, stochastic)")
def _vmc_m16_short() -> dict:
    _paths()
    import geovac.lih_vmc as V                                    # promoted (step 2b)
    wf = V.cached_wavefunction("m16")
    out: Dict[str, Any] = dict(E_fci_tot=float(wf.E_fci + wf.Vnn), Mk=int(wf.Mk), nd=int(len(wf.c)))
    E, err, acc = V.vmc(wf, nwalk=100, nsweep=120, nburn=20, jas=None, seed=10, verbose=False)
    out.update(gate6_E=float(E), gate6_err=float(err), gate6_acc=float(acc))
    E, err, acc = V.vmc(wf, nwalk=100, nsweep=120, nburn=20, jas=V.PadeJastrow(1.5, 2, 2), seed=10,
                        verbose=False)
    out.update(pade_E=float(E), pade_err=float(err), pade_acc=float(acc))
    return out


# =========================================================================== #
# D. Bridge (Phase 1, gate G-leaf)  (tests/test_lih_r12ci_bridge.py)
# =========================================================================== #
@item("bridge_reductions", "fast",
      "tests/test_lih_r12ci_bridge.py helpers on debug/data/lih_marriage_phase1.npz: "
      "_reduce(D)=grid_int(D*neumann_potential(D,34,exact=True)); isotropic-shortcut + L<=2 fire-test values",
      notes="the bridge guard's recomputed quantities (exact-Neumann bridge on the frozen dressings)")
def _bridge_red() -> dict:
    _paths()
    T = _load_module(os.path.join(TESTS, "test_lih_r12ci_bridge.py"), "_golden_bridge_test")
    import geovac.lih_r12ci.kernels as KG
    art = dict(np.load(T.OUT, allow_pickle=True))
    alpha = float(art["ALPHA"])
    rho_bridge = (alpha ** 3 / np.pi) * np.exp(-2.0 * alpha * KG.rB.ravel())
    out: Dict[str, Any] = {}
    out["C0_I_prod_recomputed"] = T._reduce(rho_bridge * art["C0_Psi_prod"])
    out["C1_I_prod_recomputed"] = T._reduce(art["rho_11"] * art["C1_Psi_prod"])
    out["C1_I_conv_recomputed"] = T._reduce(art["rho_11"] * art["C1_Psi_hi"])
    for k in ("C0_I_prod", "C1_I_prod", "C1_I_conv", "C1_I_mc", "C1_s_mc"):
        out[f"artifact_{k}"] = float(art[k])
    gam = float(art["GAM"]); R = float(art["R"])
    mod = T._load_driver()
    rb = mod.graded_bounds(0.0, 30.0, [0.0, R], 0.004, 1.7)
    rs, wr = mod.gl_on_panels(rb, 10)
    zeta = 2.6875
    rho_1s = (zeta ** 3 / np.pi) * np.exp(-2.0 * zeta * rs)
    s_probe = np.array([0.3, 1.0, 2.5, 6.0])
    out["iso_probe"] = T._isotropic_dressing(s_probe, rs, wr, rho_1s, gam)
    out["iso_closed"] = mod.psi_f_iso(s_probe, zeta, gam)
    Psi_iso = T._isotropic_dressing(KG.rA.ravel(), rs, wr, art["C1_leaf_sph_rho"], gam)
    out["C1_I_isotropic_shortcut"] = T._reduce(art["rho_11"] * Psi_iso)
    cl = np.asarray(art["C1_cl_hi"], dtype=float)
    out["C1_I_L2"] = float(cl[:3].sum()); out["C1_I_L34"] = float(cl[:35].sum())
    return out


@item("bridge_phase1_dressings", "slow",
      "lih_marriage_phase1.main(['--configs','C0,C1','--skip-mc','--no-scratch','--no-cache','--workers','12']) "
      "with lih_marriage_phase1.OUT / .CACHE redirected to a temp dir",
      notes="regenerates the Phase-1 dressings (direct 3-D quadrature R1/RC + production gVee.dress) from "
            "scratch and compares with the frozen debug/data/lih_marriage_phase1.npz (MC skipped)")
def _bridge_dress() -> dict:
    _paths()
    import geovac.lih_r12ci.bridge_gate as P1                     # promoted (step 2c)
    tmp = tempfile.mkdtemp(prefix="golden_p1_")
    P1.OUT = os.path.join(tmp, "p1.npz")
    P1.CACHE = os.path.join(tmp, "p1_cache.npz")
    rc = P1.main(['--configs', 'C0,C1', '--skip-mc', '--no-scratch', '--no-cache', '--workers', '12'])
    new = dict(np.load(P1.OUT, allow_pickle=True))
    old = dict(np.load(os.path.join(DATA, "lih_marriage_phase1.npz"), allow_pickle=True))
    out: Dict[str, Any] = dict(main_rc=rc if rc is not None else -1)
    for k in sorted(new):
        v = new[k]
        if v.dtype == object or v.dtype.kind in "US":
            continue
        if v.ndim == 0:
            out[f"new.{k}"] = v.item()
        else:
            out[f"new.{k}"] = v
        if k in old and old[k].shape == v.shape and old[k].dtype.kind in "fc" and v.dtype.kind in "fc":
            d = np.max(np.abs(np.asarray(old[k]) - v)) if v.size else 0.0
            out[f"vs_frozen.maxabs.{k}"] = float(d)
    return out


# =========================================================================== #
# E. exact Neumann operator (already promoted; regression baseline)
# =========================================================================== #
@item("neumann_exact", "fast",
      "tests/test_lih_r12ci_neumann_exact.py helpers: _self_coulomb(zeta,RA,exact) for zeta in (4.5,2.6875); "
      "(aa|bb) dressings; m=1 coul_mode_potential exact/legacy vs solid-harmonic closed form",
      notes="geovac/lih_r12ci regression baseline (no debug/ dependency)")
def _neumann_exact() -> dict:
    _paths()
    T = _load_module(os.path.join(TESTS, "test_lih_r12ci_neumann_exact.py"), "_golden_ne_test")
    import geovac.lih_r12ci.kernels as KG
    from geovac.lih_r12ci.hVee import neumann_potential, d_aa, d_bb
    from geovac.lih_r12ci.triangle import coul_mode_potential
    out: Dict[str, Any] = {}
    for zeta in (4.5, 2.6875):
        out[f"selfcoul_exact_z{zeta}"] = T._self_coulomb(zeta, T.RA, exact=True)
        out[f"selfcoul_legacy_z{zeta}"] = T._self_coulomb(zeta, T.RA, exact=False)
    out["aabb_anchor"] = float(KG.grid_int(d_bb * T._hartree_1s(T.RA, T.ZA)))
    V_aa = neumann_potential(d_aa.reshape(T.SHAPE), exact=True).reshape(-1)
    V_bb = neumann_potential(d_bb.reshape(T.SHAPE), exact=True).reshape(-1)
    out["aabb_dress_aa"] = float(KG.grid_int(d_bb * V_aa))
    out["aabb_dress_bb"] = float(KG.grid_int(d_aa * V_bb))
    zeta, m = 2.6875, 1
    Bm = T.RHO_CYL ** m * np.exp(-2.0 * zeta * T.RA)
    Ve = coul_mode_potential(Bm.reshape(T.SHAPE), m, exact=True).reshape(-1)
    Vl = coul_mode_potential(Bm.reshape(T.SHAPE), m, exact=False).reshape(-1)
    out["mode1_exact"] = Ve; out["mode1_legacy"] = Vl
    out["mode1_closed"] = T._solid_harmonic_potential(T.RA, zeta, m)
    return out


@item("fourbody_sigma_reduced", "fast",
      "lih_r12_4body_integral.reduced(); control_norm_rhoB(); control_neumann_selfcoulomb() "
      "[module builds its 260x100 grid at import]",
      notes="the two-centre sigma 4-body bridge gate (closed-form 5 alpha/8 anchor); "
            "geovac/lih_r12ci/fourbody.py is a verbatim copy of this module")
def _fourbody() -> dict:
    _paths()
    import geovac.lih_r12ci.fourbody as FB                        # single home (step 2c)
    I, per_l = FB.reduced()
    sc, sc_per_l = FB.control_neumann_selfcoulomb()[0], None
    return dict(reduced=float(I), reduced_per_l=np.asarray(per_l, dtype=float),
                norm_rhoB=float(FB.control_norm_rhoB()), selfcoul=float(sc),
                selfcoul_closed=float(5 * FB.ALPHA / 8.0))


# =========================================================================== #
# F. geovac/lih_r12ci 2-MO engine (already promoted; regression baseline)
# =========================================================================== #
def _r12_result(r) -> dict:
    out: Dict[str, Any] = {}
    for k, v in vars(r).items():
        if k == "pieces":
            out["pieces"] = {kk: float(vv) for kk, vv in v.items() if isinstance(vv, (int, float, np.floating))}
        elif isinstance(v, (bool, int, float, np.floating, np.bool_)):
            out[k] = v
    return out


@item("r12ci_2mo_legacy", "slow",
      "geovac.lih_r12ci.energy('exp'); energy('linexp')  [default USE_EXACT_NEUMANN=False]",
      notes="tests/test_lih_r12ci.py (banked -7.9420 / -7.9168)")
def _r12_legacy() -> dict:
    _paths()
    from geovac.lih_r12ci.assembly import energy
    import geovac.lih_r12ci.kernels as KG
    assert KG.USE_EXACT_NEUMANN is False
    return {"exp": _r12_result(energy("exp")), "linexp": _r12_result(energy("linexp"))}


@item("r12ci_2mo_exact", "slow",
      "geovac.lih_r12ci.kernels.USE_EXACT_NEUMANN=True; energy('linexp'); energy('exp')",
      notes="Phase-0b exact-kernel 2-MO energies (-7.917199 linexp / -7.943605 exp)")
def _r12_exact() -> dict:
    _paths()
    import geovac.lih_r12ci.kernels as KG
    KG.USE_EXACT_NEUMANN = True
    from geovac.lih_r12ci.assembly import energy
    return {"linexp": _r12_result(energy("linexp")), "exp": _r12_result(energy("exp"))}


# =========================================================================== #
# G. Paper 19 pure-product core  (tests/test_paper19_pure_product_core.py; geovac-only)
# =========================================================================== #
@item("pure_product_core", "fast",
      "tests/test_paper19_pure_product_core.cc_error_mha(N,'3.015',40+8N) for N in (15,10)",
      notes="(cc|cc)-15/8 in mHa and the 1s energy error; geovac-only")
def _ppc() -> dict:
    T = _load_module(os.path.join(TESTS, "test_paper19_pure_product_core.py"), "_golden_ppc")
    e15 = T.cc_error_mha(15, "3.015", T._dps(15))
    e10 = T.cc_error_mha(10, "3.015", T._dps(10))
    return dict(N15=[float(x) for x in e15], N10=[float(x) for x in e10])


@item("pure_product_core_rspread", "slow",
      "tests/test_paper19_pure_product_core.cc_error_mha(15,R,160) for R in ('2.5','3.5')",
      notes="test_r_spread_at_15x15 (slow); the R=3.015 point is in item pure_product_core")
def _ppc_r() -> dict:
    T = _load_module(os.path.join(TESTS, "test_paper19_pure_product_core.py"), "_golden_ppc")
    return {f"R{R}": [float(x) for x in T.cc_error_mha(15, R, T._dps(15))] for R in ("2.5", "3.5")}


# =========================================================================== #
# H. sigma marriage (Phase 0 / 2 / 3)
# =========================================================================== #
_P0_KEYS_SKIP = {"prims", "labels", "cfg", "grid_params", "verdict", "no_char", "g0_rows"}


@item("marriage_phase0_sigma_ref", "fast",
      "lih_marriage_phase0.main(exact=True) with OUT/OUT_EXACT redirected to a temp dir "
      "[uses the committed prim cache debug/data/lih_marriage_phase0_prim.npz]",
      notes="74-det sigma reference: E_full -8.023354, E_trunc -8.022918, pair tensor T, NO grid densities. "
            "Compared array-by-array with the committed debug/data/lih_marriage_phase0_ref_exact.npz. "
            "SIDE EFFECT: grid_setup() sets geovac.lih_r12ci.kernels.USE_EXACT_NEUMANN")
def _p0() -> dict:
    _paths()
    import geovac.lih_r12ci.ref_tensor as P0                      # promoted (step 2c)
    tmp = tempfile.mkdtemp(prefix="golden_p0_")
    P0.OUT = os.path.join(tmp, "ref.npz"); P0.OUT_EXACT = os.path.join(tmp, "ref_exact.npz")
    verdict = P0.main(exact=True)
    new = dict(np.load(P0.OUT_EXACT, allow_pickle=True))
    old = dict(np.load(os.path.join(DATA, "lih_marriage_phase0_ref_exact.npz"), allow_pickle=True))
    out: Dict[str, Any] = dict(verdict=str(verdict))
    for k in sorted(new):
        if k in _P0_KEYS_SKIP:
            continue
        v = new[k]
        if v.dtype == object or v.dtype.kind in "US":
            continue
        out[f"new.{k}"] = v.item() if v.ndim == 0 else (ShaOnly(v) if k in ("eri_prim",) else v)
        if k in old and old[k].shape == v.shape and old[k].dtype.kind in "fcib" and v.dtype.kind in "fcib":
            ov = np.asarray(old[k]); nv = np.asarray(v)
            out[f"vs_committed.bitequal.{k}"] = bool(np.array_equal(ov, nv))
            if v.dtype.kind in "fc" and v.size:
                out[f"vs_committed.maxabs.{k}"] = float(np.max(np.abs(ov - nv)))
    return out


@item("marriage_phase0_prims_scratch", "slow",
      "lih_marriage_phase0.PRIM_CACHE -> temp (absent); lih_marriage_phase0.build_primitives(verbose=False) "
      "[one_body_general + prolate_float_eri.build_eri_tensor_m_f from scratch]",
      notes="from-scratch sigma M=12 primitive S/h1/T1/ERI vs the committed cache")
def _p0_prims() -> dict:
    _paths()
    import geovac.lih_r12ci.ref_tensor as P0                      # promoted (step 2c)
    tmp = tempfile.mkdtemp(prefix="golden_p0p_")
    P0.PRIM_CACHE = os.path.join(tmp, "prim.npz")
    ref = P0.build_primitives(verbose=False)
    z = np.load(os.path.join(DATA, "lih_marriage_phase0_prim.npz"), allow_pickle=True)
    out: Dict[str, Any] = dict(M=int(ref['M']), Mk=int(ref['Mk']), cond=float(ref['cond']))
    for k in ("S", "h1", "T1", "eri"):
        a = np.asarray(ref[k])
        out[k] = a
        out[f"vs_cache.maxabs.{k}"] = float(np.max(np.abs(a - z[k])))
        out[f"vs_cache.bitequal.{k}"] = bool(np.array_equal(a, z[k]))
    for k in ("h1_f", "eri_f", "Tmap"):
        out[k] = np.asarray(ref[k])
    return out


@item("marriage_2mo_generalized", "slow",
      "KG.USE_EXACT_NEUMANN=True; lih_marriage_phase2.build_ref_2mo(); for g in (linexp,exp): "
      "P2.all_pieces(P2.Geminal(g),ref2,T0,verbose=False); P2.assemble(p,-7.887822) and P2.assemble(p,p['E0_grid'])",
      notes="Phase-2 gate G1 / Phase-3 M1a: the 2-MO reference through the generalized enumerators "
            "(frozen exact-path -7.917198 linexp / -7.943536 exp). Import of lih_marriage_phase2 ~2 min "
            "(gVee builds its phi kernels at import)")
def _m2mo() -> dict:
    _paths()
    import geovac.lih_r12ci.kernels as KG
    import geovac.lih_r12ci.enumerators as P2                     # promoted (step 2c)
    from geovac.lih_r12ci import gVee as GV
    from geovac.lih_r12ci.triangle import _MODE_CACHE
    KG.USE_EXACT_NEUMANN = True
    GV._NEU_CACHE.clear(); _MODE_CACHE.clear(); P2._TRI_MAT_CACHE.clear()
    ref2 = P2.build_ref_2mo()
    out: Dict[str, Any] = {}
    for g in ("linexp", "exp"):
        p = P2.all_pieces(P2.Geminal(g), ref2, time.time(), verbose=False)
        Ef, h = P2.assemble(p, -7.887822)
        Eg, _ = P2.assemble(p, p['E0_grid'])
        out[g] = dict(pieces={k: float(v) for k, v in p.items() if isinstance(v, (int, float, np.floating))},
                      E_R12_fixedE0=float(Ef), E_R12_gridE0=float(Eg), h=float(h))
    return out


def _sigma74(gname: str) -> dict:
    _paths()
    import geovac.lih_r12ci.kernels as KG
    import geovac.lih_r12ci.enumerators as P2                     # promoted (step 2c)
    import geovac.lih_r12ci.sigma_marriage as P3                  # promoted (step 2c)
    from geovac.lih_r12ci import gVee as GV
    from geovac.lih_r12ci.triangle import _MODE_CACHE
    KG.USE_EXACT_NEUMANN = True
    GV._NEU_CACHE.clear(); _MODE_CACHE.clear(); P2._TRI_MAT_CACHE.clear()
    refA = P2.build_ref_artifact()
    p = P2.all_pieces(P2.Geminal(gname), refA, time.time(), verbose=False)
    E_R12, c_mix, info = P3.mixing_c(p, p['E0_grid'], P2.V_NN)
    return dict(E_R12=float(E_R12), c_mix=float(c_mix), nact=int(refA.nact),
                pieces={k: float(v) for k, v in p.items() if isinstance(v, (int, float, np.floating))},
                info={k: float(v) for k, v in info.items()})


for _g, _banked in (("linexp", "-8.037236"), ("exp", "-8.037088")):
    item(f"marriage_sigma74_{_g}", "slow",
         f"KG.USE_EXACT_NEUMANN=True; refA=lih_marriage_phase2.build_ref_artifact() [reads "
         f"debug/data/lih_marriage_phase0_ref_exact.npz]; p=P2.all_pieces(P2.Geminal('{_g}'),refA,T0); "
         f"lih_marriage_phase3.mixing_c(p,p['E0_grid'],P2.V_NN)  [fresh process, this geminal only]",
         notes=f"the 74-det sigma marriage E_R12 (banked {_banked}, Paper 19). The banked value came from a "
               f"process that ran M1a (2-MO) first; the triangle's randomized eigh draws from "
               f"geovac.lih_r12ci.triangle.rng, so a different call order moves the result at the ~1e-7 level")(
        (lambda gg: (lambda: _sigma74(gg)))(_g))


@item("marriage_sigma74_vmc_short", "fast",
      "lih_marriage_phase3.build_truncated_wf(); lih_vmc.vmc_linear(wf,LinearGeminalJastrow(g,0.5,c,Fbar,2,2),"
      "nwalk=200,nsweep=200,nburn=50,seed=10) with c,Fbar from debug/data/lih_marriage_phase3.npz (m1_*); "
      "plus c=0 gate-6",
      seeds={"vmc_seed": 10},
      notes="G2 machinery on the 74-det reference at short, exactly reproducible statistics")
def _s74_vmc() -> dict:
    _paths()
    import geovac.lih_vmc as V                                    # promoted (step 2b)
    import geovac.lih_r12ci.sigma_marriage as P3                  # promoted (step 2c)
    wf, E_trunc, meta = P3.build_truncated_wf()
    z = np.load(os.path.join(DATA, "lih_marriage_phase3.npz"), allow_pickle=True)
    out: Dict[str, Any] = dict(E_trunc=float(E_trunc), meta={k: int(v) for k, v in meta.items()})
    lin0 = V.LinearGeminalJastrow('linexp', 0.5, 0.0, 0.0, wf.na, wf.nb)
    E, err, acc = V.vmc_linear(wf, lin0, nwalk=200, nsweep=200, nburn=50, seed=10, verbose=False)
    out.update(gate6_E=float(E), gate6_err=float(err))
    for g in ("linexp", "exp"):
        c = float(z[f"m1_{g}_c_mix"]); Fb = float(z[f"m1_{g}_Fbar"])
        lin = V.LinearGeminalJastrow(g, 0.5, c, Fb, wf.na, wf.nb)
        E, err, acc = V.vmc_linear(wf, lin, nwalk=200, nsweep=200, nburn=50, seed=10, verbose=False)
        out[g] = dict(c=c, Fbar=Fb, E=float(E), err=float(err), acc=float(acc))
    return out


# =========================================================================== #
# I. pi / M=16 marriage
# =========================================================================== #
def _pi_ref(trunc: float):
    import geovac.lih_r12ci.pi_m16 as M                           # promoted (step 2c)
    ref = M.build_primitives(M.CFG, M.PRIM_CACHE_PI, verbose=False)
    refm, info = M.build_ref_from_prims(ref, "M=16 pi natural-orbital", cfg_tag='pi', trunc=trunc,
                                        verbose=False)
    return M, ref, refm, info


@item("pi_ref_truncations", "slow",
      "lih_pi_m16.build_primitives(CFG,PRIM_CACHE_PI) [cached debug/data/lih_pi_m16_prim.npz]; "
      "build_ref_from_prims(ref,tag,cfg_tag='pi',trunc=t) for t in (5e-2,1e-2,5e-3,1e-3) "
      "[cached FCI vectors debug/data/lih_pi_m16_fci_pi_{canon,no}.npz]",
      notes="the pi reference at each truncation: E_full -8.02905, E_trunc (-8.025184 @1e-2, -8.028633 @1e-3), "
            "nact, the m-resolved tensor T. SIDE EFFECT: sets KG.USE_EXACT_NEUMANN=True (via P0.grid_setup)")
def _pi_refs() -> dict:
    _paths()
    import geovac.lih_r12ci.pi_m16 as M                           # promoted (step 2c)
    import geovac.lih_r12ci.kernels as KG
    ref = M.build_primitives(M.CFG, M.PRIM_CACHE_PI, verbose=False)
    out: Dict[str, Any] = dict(Mk=int(ref['Mk']), M=int(ref['M']), cond=float(ref['cond']),
                               m_prim=ref['m_prim'], m_mo=ref['m_mo'], h1_f=ref['h1_f'],
                               eri_f=ShaOnly(ref['eri_f']), Tmap=ref['Tmap'])
    for t in (5e-2, 1e-2, 5e-3, 1e-3):
        refm, info = M.build_ref_from_prims(ref, "golden", cfg_tag='pi', trunc=t, verbose=False)
        liA, liB, liC, liD, val = refm.T_iA, refm.T_iB, refm.T_iC, refm.T_iD, refm.T_val
        d = dict(E_full=float(info['E_full']), E_trunc=float(info['E_trunc']), nkept=int(len(info['kept'])),
                 nact=int(refm.nact), npi_pairs=int(np.sum(refm.mA != 0)), T_nnz=int(len(val)),
                 occ=info['occ'], m_no=info['m_no'], T_no=info['T_no'], act=np.asarray(refm.act),
                 mA=np.asarray(refm.mA), T_val=np.asarray(val),
                 T_idx=np.stack([np.asarray(liA), np.asarray(liB), np.asarray(liC), np.asarray(liD)]),
                 rho_ord=ShaOnly(refm.rho), grad_ord=ShaOnly(refm.grad))
        out[f"trunc_{t:.0e}"] = d
    out["USE_EXACT_NEUMANN_after"] = bool(KG.USE_EXACT_NEUMANN)
    return out


@item("pi_prims_scratch", "slow",
      "lih_pi_m16.build_primitives(lih_pi_m16.CFG, <temp path>) [from scratch: one_body_general + "
      "prolate_float_eri.build_eri_tensor_m_f] vs debug/data/lih_pi_m16_prim.npz; then "
      "lih_pi_m16.fci_solve(h1_f,eri_f,Mk) vs debug/data/lih_pi_m16_fci_pi_canon.npz",
      notes="from-scratch M=16 pi primitives + 14400-det FCI (eigsh tol=0, sign fixed by max-|c| > 0)")
def _pi_prims() -> dict:
    _paths()
    import geovac.lih_r12ci.pi_m16 as M                           # promoted (step 2c)
    tmp = tempfile.mkdtemp(prefix="golden_pip_")
    ref = M.build_primitives(M.CFG, os.path.join(tmp, "prim.npz"), verbose=False)
    z = np.load(M.PRIM_CACHE_PI, allow_pickle=True)
    out: Dict[str, Any] = dict(Mk=int(ref['Mk']), cond=float(ref['cond']))
    for k in ("S", "h1", "T1", "eri"):
        a_new = np.load(os.path.join(tmp, "prim.npz"), allow_pickle=True)[k]
        out[k] = ShaOnly(a_new) if k == "eri" else a_new
        out[f"vs_cache.maxabs.{k}"] = float(np.max(np.abs(a_new - z[k])))
        out[f"vs_cache.bitequal.{k}"] = bool(np.array_equal(a_new, z[k]))
    E0, c, dets, H = M.fci_solve(ref['h1_f'], ref['eri_f'], ref['Mk'])
    zc = np.load(os.path.join(DATA, "lih_pi_m16_fci_pi_canon.npz"))
    out.update(fci_E0=float(E0), fci_E_tot=float(E0 + ref['Vnn']), fci_c=np.asarray(c),
               cached_fci_E0=float(zc['E0']), fci_maxabs_dc_vs_cache=float(np.max(np.abs(c - zc['c']))))
    return out


def _gvee_5e2_setup():
    _paths()
    import geovac.lih_r12ci.kernels as KG
    import geovac.lih_r12ci.triangle as TRI
    import geovac.lih_r12ci.pi_m16 as M                           # promoted (step 2c)
    import geovac.lih_r12ci.gvee_sparse as S                      # promoted (step 2c)
    KG.USE_EXACT_NEUMANN = True
    ref = M.build_primitives(M.CFG, M.PRIM_CACHE_PI, verbose=False)
    refm, _ = M.build_ref_from_prims(ref, "guard pi tiny", cfg_tag='pi', trunc=5e-2, verbose=False)
    gem = M.GeminalM('linexp')
    Fbar = M.fbar(gem, refm)

    def reset_rng():
        TRI.rng = np.random.default_rng(20260921)
        M._RSVD_RNG = np.random.default_rng(20260923)
        TRI._MODE_CACHE.clear()
        M._TRIM_CACHE.clear()
    return M, S, gem, refm, Fbar, reset_rng


@item("gvee_sparse_vs_dense_5e-2", "slow",
      "tests/test_lih_gvee_sparse.py fixture: KG.USE_EXACT_NEUMANN=True; pi ref trunc 5e-2; GeminalM('linexp'); "
      "reset_rng(); lih_pi_m16.hg_Vee_m(gem,refm,Fbar,force_general=True); reset_rng(); "
      "lih_gvee_sparse.hg_Vee_sparse(gem,refm,Fbar,mem_cap_gb=2.0|0.02); _DROP_NESTED fire test",
      seeds={"triangle.rng": 20260921, "lih_pi_m16._RSVD_RNG": 20260923},
      notes="(h_Vee,g_Vee,Vee,FVee,F2Vee) dense vs sparse; equal to roundoff only with identical RNG streams")
def _gvee_sd() -> dict:
    M, S, gem, refm, Fbar, reset_rng = _gvee_5e2_setup()
    names = ['h_Vee', 'g_Vee', 'Vee', 'FVee', 'F2Vee']
    out: Dict[str, Any] = dict(Fbar=float(Fbar), nact=int(refm.nact))
    reset_rng(); dense = M.hg_Vee_m(gem, refm, Fbar, force_general=True)
    reset_rng(); sp = S.hg_Vee_sparse(gem, refm, Fbar, mem_cap_gb=2.0)
    reset_rng(); sp_small = S.hg_Vee_sparse(gem, refm, Fbar, mem_cap_gb=0.02)
    S._DROP_NESTED = True
    try:
        reset_rng(); sp_drop = S.hg_Vee_sparse(gem, refm, Fbar, mem_cap_gb=2.0)
    finally:
        S._DROP_NESTED = False
    out["dense"] = dict(zip(names, [float(x) for x in dense]))
    out["sparse"] = dict(zip(names, [float(x) for x in sp]))
    out["sparse_tinychunks"] = dict(zip(names, [float(x) for x in sp_small]))
    out["sparse_drop_tree4"] = dict(zip(names, [float(x) for x in sp_drop]))
    return out


@item("tri_symmetry_5e-2", "slow",
      "tests/test_lih_gvee_triangle_symmetry.py: pi ref trunc 5e-2, GeminalM('linexp'); for each ORBIT triple: "
      "S.Ctx(gem,refm,1.0); S.TriCollector(); S._triple_args; S.vee_reduce_sparse; S.eval_triangles; "
      "with M._WCOMB_LEGACY False and True; M.tri_wcomb table",
      notes="no RNG reset in the test (the orbit values use the RNG state as it falls)")
def _trisym() -> dict:
    _paths()
    import geovac.lih_r12ci.kernels as KG
    import geovac.lih_r12ci.pi_m16 as M                           # promoted (step 2c)
    import geovac.lih_r12ci.gvee_sparse as S                      # promoted (step 2c)
    KG.USE_EXACT_NEUMANN = True
    ref = M.build_primitives(M.CFG, M.PRIM_CACHE_PI, verbose=False)
    refm, _ = M.build_ref_from_prims(ref, "guard pi tiny tri-sym", cfg_tag='pi', trunc=5e-2, verbose=False)
    gem = M.GeminalM('linexp')
    ORBIT = [((1, 2), (1, 3), (2, 3)), ((1, 2), (1, 4), (2, 4)), ((1, 2), (2, 3), (1, 3)),
             ((1, 2), (2, 4), (1, 4)), ((3, 4), (1, 3), (1, 4)), ((3, 4), (1, 4), (1, 3)),
             ((3, 4), (2, 3), (2, 4)), ((3, 4), (2, 4), (2, 3))]

    def orbit():
        vals = []
        for p, q, r in ORBIT:
            ctx = S.Ctx(gem, refm, 1.0)
            col = S.TriCollector()
            rem_f, kk, gk = S._triple_args(p, q, r, gem)
            assert S.vee_reduce_sparse(rem_f, r, kk, gk, ctx, col) == 0.0 and col.recs
            vals.append(S.eval_triangles(ctx, col))
        return np.array(vals, dtype=float)
    wtab = {f"{A},{sa},{B},{sb}": float(M.tri_wcomb(A, sa, B, sb))
            for A in range(0, 5) for sa in (1, -1) for B in range(0, 5) for sb in (1, -1)}
    v = orbit()
    M._WCOMB_LEGACY = True
    try:
        v_leg = orbit()
    finally:
        M._WCOMB_LEGACY = False
    return dict(orbit=v, orbit_legacy=v_leg, tri_wcomb=wtab)


def _pi_marriage(gname: str, trunc: float, gam: Optional[float], mem_cap: float = 6.0) -> dict:
    """Mirror of debug/lih_pi_m16_energy_sparse.py (gam None) / debug/lih_exp_geminal_run.py (gam set),
    one geminal, fresh process, nothing written to debug/data."""
    _paths()
    if gam is not None:
        from geovac.lih_r12ci import configure_geminal   # real API (replaced the step-2c source override)
        configure_geminal(gam)                     # before the consumer imports: the golden build order
    import geovac.lih_r12ci.kernels as KG
    import geovac.lih_r12ci.pi_m16 as M                           # promoted (step 2c)
    import geovac.lih_r12ci.gvee_sparse as S                      # promoted (step 2c)
    if gam is not None:
        from geovac.lih_r12ci import gVee as GV
        import geovac.lih_r12ci.enumerators as P2                     # promoted (step 2c)
        assert P2.GAM == M.GAM == GV.GAM == KG.GAM == gam
    T0 = time.time()
    KG.USE_EXACT_NEUMANN = True
    ref = M.build_primitives(M.CFG, M.PRIM_CACHE_PI, verbose=False)
    refm, info = M.build_ref_from_prims(ref, "M=16 pi natural-orbital", cfg_tag='pi', trunc=trunc,
                                        verbose=False)
    gem = M.GeminalM(gname)
    p = S.all_pieces_sparse(gem, refm, T0, mem_cap_gb=mem_cap)
    E_R12, c, h = M.assemble_c(p, p['E0_grid'])
    return dict(E_R12=float(E_R12), c=float(c), h=float(h), E_trunc=float(info['E_trunc']),
                E_full=float(info['E_full']), nact=int(refm.nact), GAM=float(M.GAM), trunc=trunc,
                mem_cap_gb=mem_cap,
                pieces={k: float(v) for k, v in p.items() if isinstance(v, (int, float, np.floating))})


_PI_CALL = ("KG.USE_EXACT_NEUMANN=True; ref=lih_pi_m16.build_primitives(CFG,PRIM_CACHE_PI); "
            "refm,info=build_ref_from_prims(ref,tag,cfg_tag='pi',trunc={t}); "
            "p=lih_gvee_sparse.all_pieces_sparse(lih_pi_m16.GeminalM('{g}'),refm,T0,mem_cap_gb=6.0); "
            "lih_pi_m16.assemble_c(p,p['E0_grid'])")
_XR = ("geovac.lih_r12ci.configure_geminal({gam}) [before the consumer imports; golden captured with "
       "the equivalent kernels.py GAM source override lih_exp_geminal_run.install_kernels({gam},None)]; ")

item("pi_marriage_linexp_1e-2", "slow", _PI_CALL.format(t="1e-2", g="linexp"),
     notes="linexp r e^{-0.5 r} on the 9-det pi ref. Banked pre-wcomb-fix -8.042181 (Sep 23, dense); "
           "post-fix Paper 19 table -8.042304")(lambda: _pi_marriage("linexp", 1e-2, None))
item("pi_marriage_exp1.5_1e-2", "slow", _XR.format(gam=1.5) + _PI_CALL.format(t="1e-2", g="exp"),
     notes="Slater e^{-1.5 r} on the 9-det pi ref: Paper 19 Table tab:geminal_gamma -8.048141")(
    lambda: _pi_marriage("exp", 1e-2, 1.5))
item("pi_marriage_exp1.5_1e-3", "slow", _XR.format(gam=1.5) + _PI_CALL.format(t="1e-3", g="exp"),
     notes="THE corpus-best LiH energy -8.04905 (Paper 19; -8.049045089 post-fix). ~15 GB peak, ~3 h")(
    lambda: _pi_marriage("exp", 1e-3, 1.5))
item("pi_marriage_linexp_1e-3", "slow", _PI_CALL.format(t="1e-3", g="linexp"),
     notes="linexp on the full pi ref: banked -8.043494 (Sep 25, PRE triangle-weight fix). ~9 GB peak, ~2.5 h")(
    lambda: _pi_marriage("linexp", 1e-3, None))


@item("pi_vmc_short_1e-2", "slow",
      "pi ref trunc 1e-2 (cached); wf=lih_pi_m16_vmc._build_wf_from_no(info['ct'],info['T_no'],Mk,E_trunc); "
      "lih_vmc.vmc(wf,nwalk=100,nsweep=120,nburn=20,jas=None,seed=0); "
      "lih_vmc.vmc_linear(wf,LinearGeminalJastrow('exp',1.5,-0.2269,1.78,2,2),nwalk=100,nsweep=120,nburn=20,seed=0)",
      seeds={"vmc_seed": 0},
      notes="G2-pi machinery (complex pi orbitals) at short, exactly reproducible statistics; c, Fbar are fixed "
            "illustrative values, not the analytic optimum")
def _pi_vmc() -> dict:
    _paths()
    import geovac.lih_vmc as V                                    # promoted (step 2b)
    M, ref, refm, info = _pi_ref(1e-2)
    wf = M.build_wf_from_no(info['ct'], info['T_no'], ref['Mk'], info['E_trunc'])
    out: Dict[str, Any] = dict(E_trunc=float(info['E_trunc']))
    E, err, acc = V.vmc(wf, nwalk=100, nsweep=120, nburn=20, jas=None, seed=0, verbose=False)
    out.update(gate6_E=float(E), gate6_err=float(err))
    lin = V.LinearGeminalJastrow('exp', 1.5, -0.2269, 1.78, wf.na, wf.nb)
    E, err, acc = V.vmc_linear(wf, lin, nwalk=100, nsweep=120, nburn=20, seed=0, verbose=False)
    out.update(linear_E=float(E), linear_err=float(err), linear_acc=float(acc))
    return out


# =========================================================================== #
# runner / compare / merge
# =========================================================================== #
def run_item(name: str, stage: str, tag: str = "") -> int:
    os.makedirs(stage, exist_ok=True)
    spec = ITEMS[name]
    t0 = time.perf_counter()
    vals = spec["fn"]()
    rt = time.perf_counter() - t0
    arrays: Dict[str, np.ndarray] = {}
    rec = dict(name=name, tier=spec["tier"], call=spec["call"], seeds=spec["seeds"], notes=spec["notes"],
               runtime_s=round(rt, 2), peak_rss_gb=round(_peak_rss_gb(), 3),
               finished=datetime.datetime.now().isoformat(timespec="seconds"),
               values={k: _ser(k, v, arrays) for k, v in vals.items()})
    base = os.path.join(stage, name + (f".{tag}" if tag else ""))
    np.savez(base + ".npz", **arrays)
    with open(base + ".json", "w", encoding="utf-8") as f:
        json.dump(rec, f, indent=1)
    print(f"[capture] {name}{'.' + tag if tag else ''}: {rt:.1f} s, peak {rec['peak_rss_gb']} GB -> {base}.json")
    return 0


def _flatten(prefix: str, rec: dict, out: dict) -> None:
    t = rec["type"]
    if t == "dict":
        for k, v in rec["value"].items():
            _flatten(f"{prefix}.{k}", v, out)
    elif t == "list":
        for i, v in enumerate(rec["value"]):
            _flatten(f"{prefix}[{i}]", v, out)
    else:
        out[prefix] = rec


def compare(name: str, stage: str, tag: str = "rep2") -> dict:
    with open(os.path.join(stage, name + ".json"), encoding="utf-8") as f:
        a = json.load(f)
    with open(os.path.join(stage, f"{name}.{tag}.json"), encoding="utf-8") as f:
        b = json.load(f)
    za = np.load(os.path.join(stage, name + ".npz"))
    zb = np.load(os.path.join(stage, f"{name}.{tag}.npz"))
    fa: dict = {}; fb: dict = {}
    for k, v in a["values"].items():
        _flatten(k, v, fa)
    for k, v in b["values"].items():
        _flatten(k, v, fb)
    diffs = {}
    for k in sorted(set(fa) | set(fb)):
        ra, rb = fa.get(k), fb.get(k)
        if ra is None or rb is None:
            diffs[k] = "missing in one run"
            continue
        if ra["type"] == "float":
            if ra["repr"] != rb["repr"]:
                x, y = float(ra["repr"]), float(rb["repr"])
                diffs[k] = dict(abs=abs(x - y), rel=_rel(x, y))
        elif ra["type"] == "array":
            if ra["sha256"] != rb["sha256"]:
                d: Dict[str, Any] = dict(sha_differs=True)
                if ra.get("stored") and rb.get("stored") and k in za.files and k in zb.files:
                    xa, xb = za[k], zb[k]
                    if xa.shape == xb.shape and xa.dtype.kind in "fc":
                        d["maxabs"] = float(np.max(np.abs(xa - xb)))
                        d["maxabs_up_to_sign"] = float(min(np.max(np.abs(xa - xb)), np.max(np.abs(xa + xb))))
                        d["scale"] = float(np.max(np.abs(xa)))
                diffs[k] = d
        else:
            if ra != rb:
                diffs[k] = dict(a=str(ra)[:200], b=str(rb)[:200])
    res = dict(bit_identical=(len(diffs) == 0), n_compared=len(fa), diffs=diffs)
    with open(os.path.join(stage, f"{name}.determinism.json"), "w", encoding="utf-8") as f:
        json.dump(res, f, indent=1)
    print(f"[compare] {name}: bit_identical={res['bit_identical']} ({len(diffs)} differing of {len(fa)})")
    return res


def compare_golden(name: str, stage: str) -> dict:
    """Stage ``NAME.json`` (+ .npz) vs the merged ``golden.json`` item ``NAME`` (+ golden.npz
    arrays ``NAME__<key>``).  Same bit-level rules as ``compare``: floats by repr, arrays by
    sha256 (max-abs reported when both are stored), everything else by record equality.
    Writes ``NAME.vsgolden.json`` in the stage dir."""
    with open(os.path.join(stage, name + ".json"), encoding="utf-8") as f:
        a = json.load(f)
    with open(GOLDEN_JSON, encoding="utf-8") as f:
        b = json.load(f)["items"][name]
    za = np.load(os.path.join(stage, name + ".npz"))
    zg = np.load(GOLDEN_NPZ)
    fa: dict = {}; fb: dict = {}
    for k, v in a["values"].items():
        _flatten(k, v, fa)
    for k, v in b["values"].items():
        _flatten(k, v, fb)
    diffs = {}
    for k in sorted(set(fa) | set(fb)):
        ra, rb = fa.get(k), fb.get(k)
        if ra is None or rb is None:
            diffs[k] = "missing in one run"
            continue
        if ra["type"] == "float":
            if ra["repr"] != rb["repr"]:
                x, y = float(ra["repr"]), float(rb["repr"])
                diffs[k] = dict(new=ra["repr"], golden=rb["repr"], abs=abs(x - y), rel=_rel(x, y))
        elif ra["type"] == "array":
            if ra["sha256"] != rb["sha256"]:
                d: Dict[str, Any] = dict(sha_differs=True)
                gk = f"{name}__{k}"
                if ra.get("stored") and rb.get("stored") and k in za.files and gk in zg.files:
                    xa, xb = za[k], zg[gk]
                    if xa.shape == xb.shape and xa.dtype.kind in "fc":
                        d["maxabs"] = float(np.max(np.abs(xa - xb)))
                        d["scale"] = float(np.max(np.abs(xb)))
                diffs[k] = d
        else:
            if ra != rb:
                diffs[k] = dict(new=str(ra)[:200], golden=str(rb)[:200])
    res = dict(bit_identical=(len(diffs) == 0), n_compared=len(fa), diffs=diffs,
               runtime_s=a.get("runtime_s"), golden_runtime_s=b.get("runtime_s"))
    with open(os.path.join(stage, f"{name}.vsgolden.json"), "w", encoding="utf-8") as f:
        json.dump(res, f, indent=1)
    print(f"[compare-golden] {name}: bit_identical={res['bit_identical']} "
          f"({len(diffs)} differing of {len(fa)}); runtime {a.get('runtime_s')} s vs golden "
          f"{b.get('runtime_s')} s")
    return res


def _fit_req(rows: List[tuple]) -> Optional[float]:
    """prolate_allelectron_c4._fit_req verbatim (parabola through the 3 points around the minimum)."""
    band = [(R, E) for (R, E) in rows if 2.4 <= R <= 3.6]
    if len(band) < 3:
        return None
    Rr = np.array([b[0] for b in band]); Ee = np.array([b[1] for b in band])
    imin = int(np.argmin(Ee))
    lo, hi = max(0, imin - 1), min(len(band), imin + 2)
    if hi - lo < 3:
        lo, hi = 0, min(3, len(band))
    c = np.polyfit(Rr[lo:hi], Ee[lo:hi], 2)
    if c[0] <= 0:
        return None
    return float(-c[1] / (2 * c[0]))


def merge(stage: str) -> None:
    items: Dict[str, dict] = {}
    arrays: Dict[str, np.ndarray] = {}
    for path in sorted(glob.glob(os.path.join(stage, "*.json"))):
        base = os.path.basename(path)[:-5]
        with open(path, encoding="utf-8") as f:
            rec = json.load(f)
        if rec.get("name") != base:
            continue                                 # rep / determinism side files
        name = rec["name"]
        det_path = os.path.join(stage, f"{name}.determinism.json")
        if os.path.exists(det_path):
            with open(det_path, encoding="utf-8") as f:
                rec["determinism"] = json.load(f)
        else:
            rec["determinism"] = dict(checked=False)
        z = np.load(os.path.join(stage, name + ".npz"))
        for k in z.files:
            arrays[f"{name}__{k}"] = z[k]
        items[name] = rec
    # derived: Route C R_eq from the c4 points
    rows = []
    for Rs in C4_RS:
        r = items.get(f"routec_c4_req_R{Rs}")
        if r:
            rows.append((float(Rs), float(r["values"]["E_tot"]["repr"])))
    derived = {}
    if len(rows) >= 3:
        Req = _fit_req(rows)
        derived["routec_c4_Req"] = dict(points=[[R, repr(E)] for R, E in rows], R_eq=repr(Req),
                                        drift_pct_vs_3p015=repr(100.0 * (Req - 3.015) / 3.015) if Req else None,
                                        method="prolate_allelectron_c4._fit_req (parabola, 3 points around min)")
    try:
        commit = subprocess.check_output(["git", "-C", ROOT, "rev-parse", "HEAD"], text=True).strip()
        describe = subprocess.check_output(["git", "-C", ROOT, "describe", "--tags", "--always"], text=True).strip()
    except Exception:
        commit = describe = "unknown"
    import scipy
    import mpmath
    meta = dict(created=datetime.datetime.now().isoformat(timespec="seconds"), git_commit=commit,
                git_describe=describe, python=sys.version.split()[0], numpy=np.__version__,
                scipy=scipy.__version__, mpmath=mpmath.__version__, platform=platform.platform(),
                threads="OMP/MKL/OPENBLAS/NUMEXPR/VECLIB_NUM_THREADS=1 (set before numpy import)",
                float_format="repr() round-trip strings", arrays="golden.npz key '<item>__<value key>'",
                process_model="one fresh python process per item")
    with open(GOLDEN_JSON, "w", encoding="utf-8") as f:
        json.dump(dict(meta=meta, items=items, derived=derived), f, indent=1)
    np.savez_compressed(GOLDEN_NPZ, **arrays)
    print(f"[merge] {len(items)} items -> {GOLDEN_JSON} ; {len(arrays)} arrays -> {GOLDEN_NPZ}")


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--list", action="store_true")
    ap.add_argument("--run-item")
    ap.add_argument("--tag", default="")
    ap.add_argument("--compare")
    ap.add_argument("--compare-golden", help="stage NAME.json vs golden.json item NAME")
    ap.add_argument("--merge", action="store_true")
    ap.add_argument("--stage", default=DEFAULT_STAGE)
    args = ap.parse_args()
    if args.list:
        for n, s in ITEMS.items():
            print(f"{s['tier']:5s} {n}")
        return 0
    if args.run_item:
        return run_item(args.run_item, args.stage, args.tag)
    if args.compare:
        compare(args.compare, args.stage)
        return 0
    if args.compare_golden:
        res = compare_golden(args.compare_golden, args.stage)
        return 0 if res["bit_identical"] else 2
    if args.merge:
        merge(args.stage)
        return 0
    ap.print_help()
    return 1


if __name__ == "__main__":
    sys.exit(main())
