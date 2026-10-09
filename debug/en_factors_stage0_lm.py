r"""Stage 0 of the analytic e-n extension: the step-3 VMC linear-method harness rerun with the
PARAMETER-FREE core switch (PI ruling 2026-10-05; geovac/lih_r12ci/core_switch.py).

Same sampler, guide density, estimators and statistics as debug/pair_geminal_step3_lm.py (imported;
its `switch` is monkeypatched).  Switch keys (the `r0` slot of the step-3 class columns):
    0.6, 1.0   step-3 radial sigmoid 1/(1+(rho/r0)^6)          (correlated comparison)
    -1.0       DENS: w = |phi_c|^2 / |phi_c(R_Li)|^2            (the ruled switch; core NO density)
    -2.0       FRAC: w = n_c |phi_c|^2 / sum_p n_p |phi_p|^2    (variant: core share of the density)
    -3.0       AMP:  w = phi_c / phi_c(R_Li)                    (option-b amplitude switch, PI 2026-10-09)
Class columns G^{key}_{cls,gamma} Psi_CI, cls in {cc,cv,vv}, gamma in {0, 1.5}.  gamma = 0 spans
{W1 Psi, W2 Psi} (with Psi); gamma = 1.5 gives F_cc, F_cv, F_vv (sum = F).  Blocks: Phi_I, F1.5 Phi_I.

Usage: python debug/en_factors_stage0_lm.py --seed 1 [--ref 9det] [--nsweep 9000 --nburn 4000]
       [--keys 0.6,-1,-3 --tag enab_lm_s1]   (option b: density + amplitude switches)
       python debug/en_factors_stage0_lm.py --mode vmc --vec debug/data/pairgem_en0_vecs.npz --vkey X --seed 301
Raw block matrices -> debug/data/pairgem_en0_lm_s<seed>.npz (gitignored pattern pairgem_*_lm_s*.npz).
"""
from __future__ import annotations

import os
import sys

for _v in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS", "NUMEXPR_NUM_THREADS",
           "VECLIB_MAXIMUM_THREADS"):
    os.environ[_v] = "1"

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, ROOT); sys.path.insert(0, HERE)
DATA = os.path.join(HERE, "data")

import pair_geminal_step3_lm as S3  # noqa: E402
from geovac.lih_r12ci.core_switch import CoreSwitch, CoreFractionSwitch, AmpSwitch  # noqa: E402

SWREF = os.path.join(DATA, "en_switch_ref.npz")
KEYS = (0.6, 1.0, -1.0, -2.0)
GAMMAS = (0.0, 1.5)


def switch_reference():
    """T_no, occ, m_no of the M=16 pi reference (full-CI NOs, occupation-ordered); cached."""
    if not os.path.exists(SWREF):
        from geovac.lih_r12ci import pi_m16 as PM
        ref = PM.build_primitives(PM.CFG, PM.PRIM_CACHE_PI, verbose=False)
        _refm, info = PM.build_ref_from_prims(ref, "M=16 pi natural-orbital", cfg_tag='pi', trunc=1e-3,
                                              verbose=False)
        np.savez(SWREF, T_no=info['T_no'], occ=info['occ'], m_no=info['m_no'], R=info['R'])
    z = np.load(SWREF)
    return z['T_no'], z['occ'], z['m_no'], float(z['R'])


_SW = {}
_orig_switch = S3.switch


def install(prims):
    T_no, occ, m_no, R = switch_reference()
    zr = np.load(S3.S2.REFS["90det"], allow_pickle=True)
    assert np.allclose(zr['T_no'], T_no, atol=1e-12), "switch NOs != reference NOs"
    dens = CoreSwitch.from_reference(prims, T_no, occ, m_no, R)
    _SW[-1.0] = dens
    _SW[-2.0] = CoreFractionSwitch(prims, T_no, occ, dens.no_index, R)
    _SW[-3.0] = AmpSwitch.from_reference(prims, T_no, occ, m_no, R)
    print(f"  core NO = {dens.no_index} (occ {occ[dens.no_index]:.6f}), phi_c(Li) = {dens.phi0:+.6f}",
          flush=True)

    def switch(r, R_, r0, n=S3.NSW):
        if r0 in _SW:
            return _SW[r0].vgl(r)
        return _orig_switch(r, R_, r0, n)
    S3.switch = switch


def main():
    keys = KEYS
    if "--keys" in sys.argv:                       # e.g. --keys 0.6,-1,-3 (option-b run)
        i = sys.argv.index("--keys"); keys = tuple(float(x) for x in sys.argv[i + 1].split(","))
        del sys.argv[i:i + 2]
    sys.argv = [sys.argv[0]] + ["--nodress", "--radii", ",".join(str(k) for k in keys),
                                "--gammas", ",".join(str(g) for g in GAMMAS)] + sys.argv[1:]
    if "--tag" not in sys.argv and "--mode" not in sys.argv:
        seed = sys.argv[sys.argv.index("--seed") + 1] if "--seed" in sys.argv else "0"
        ref = sys.argv[sys.argv.index("--ref") + 1] if "--ref" in sys.argv else "90det"
        sys.argv += ["--tag", ("en0" if ref == "90det" else f"en0{ref}") + f"_lm_s{seed}"]
    _orig_load = S3.load_ref

    def load_ref(name, trunc=1e-3):
        if name == "9det":
            z = np.load(S3.S2.REFS["9det"], allow_pickle=True)
            wf9 = S3.PM.build_wf_from_no(z['ct'], z['T_no'], int(z['Mk']), float(z['E_trunc']))
            return wf9, float(z['E_trunc'])
        return _orig_load(name, trunc)
    S3.load_ref = load_ref
    ref = sys.argv[sys.argv.index("--ref") + 1] if "--ref" in sys.argv else "90det"
    wf, _ = S3.load_ref(ref)
    install(wf.prims)
    S3.main()


if __name__ == "__main__":
    main()
