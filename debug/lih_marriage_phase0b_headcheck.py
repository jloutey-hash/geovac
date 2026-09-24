r"""Phase 0b legacy BIT-IDENTITY check against the COMMITTED code (2026-09-23).

The exact prolate-Neumann operator (geovac/lih_r12ci/neumann_exact.py) was wired into
kernels.py / hVee.py / triangle.py / gVee.py as an OPT-IN path.  The load-bearing promise is that
the DEFAULT path is bit-identical to the code before the edit, so tests/test_lih_r12ci.py's
-7.9168 / -7.9420 anchors are untouched.  This script proves it against `git HEAD`, not against
a capture whose provenance has to be trusted:

  1. `git archive HEAD geovac/lih_r12ci` is extracted into a temp dir as the package `headpkg.lih_r12ci`
     (the package uses only relative imports, so the copy runs unmodified);
  2. two subprocesses dump the same 27 legacy-path arrays -- neumann_potential on d_aa/d_ab/d_bb,
     four isotropic 1s STO densities and the Phase-0 rho_(NO0,NO1); coul_mode_potential m=0..4 on
     rho01*rho_cyl^m and m=0..2 on P00/rho; psi_coul/psi_yuk; make_kernel_coul's W and Psi_ab; and
     the Stage-1 names V_aaaa..V_bbab / Psi_aa / Psi_bb (eager in HEAD, PEP-562 lazy now) --
     one from the HEAD copy, one from the working tree;
  3. max |Delta| over every array must be 0.0 exactly.  The previous attempt's capture
     debug/data/lih_marriage_phase0b_legacy_pristine.npz is compared to the HEAD dump as well, which
     establishes its provenance (it is what debug/lih_marriage_phase0b.py section 0 reads).

Result 2026-09-23: 27/27 arrays 0.0 exactly (both comparisons); HEAD `import kernels` 60.5 s vs
working tree 0.5 s (the Stage-1 build is deferred to first attribute access, 60.6 s there).

Run from root:  python debug/lih_marriage_phase0b_headcheck.py > debug/data/lih_marriage_phase0b_headcheck.log 2>&1
"""
from __future__ import annotations

import importlib
import io
import os
import shutil
import subprocess
import sys
import tarfile
import tempfile
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
DATA = os.path.join(HERE, "data")
REF = os.path.join(DATA, "lih_marriage_phase0_ref.npz")
PRISTINE = os.path.join(DATA, "lih_marriage_phase0b_legacy_pristine.npz")


def dump(pkg: str, root: str, out: str) -> None:
    """Legacy-path outputs of <pkg>.lih_r12ci (default switches) -> out.npz; prints import times."""
    sys.path.insert(0, root)
    t0 = time.time()
    KG = importlib.import_module(f"{pkg}.lih_r12ci.kernels")
    t_k = time.time() - t0
    hVee = importlib.import_module(f"{pkg}.lih_r12ci.hVee")
    tri = importlib.import_module(f"{pkg}.lih_r12ci.triangle")
    t0 = time.time()
    gVee = importlib.import_module(f"{pkg}.lih_r12ci.gVee")
    t_g = time.time() - t0
    EN = importlib.import_module(f"{pkg}.lih_r12ci.energy")
    print(f"[{pkg}] import kernels {t_k:.1f} s, gVee {t_g:.1f} s; USE_EXACT_NEUMANN = "
          f"{getattr(KG, 'USE_EXACT_NEUMANN', '(absent)')}", flush=True)
    shape = (hVee.NXI, hVee.NETA)
    rA = KG.rA.ravel(); rB = KG.rB.ravel()
    z0 = np.load(REF, allow_pickle=True)
    pairs = [tuple(p) for p in z0["pairs"]]
    rho01 = z0["rho_no"][pairs.index((0, 1))]
    now = {}
    now["NP_daa"] = hVee.neumann_potential(hVee.d_aa.reshape(shape))
    now["NP_dab"] = hVee.neumann_potential(hVee.d_ab.reshape(shape))
    now["NP_dbb"] = hVee.neumann_potential(hVee.d_bb.reshape(shape))
    for zeta, cen in ((1.0, "B"), (1.6, "A"), (2.6875, "A"), (4.5, "A")):
        r = rA if cen == "A" else rB
        now[f"NP_rho1s_{zeta}"] = hVee.neumann_potential(((zeta ** 3 / np.pi) * np.exp(-2.0 * zeta * r)).reshape(shape))
    now["NP_rho01"] = hVee.neumann_potential(rho01.reshape(shape))
    rc2 = np.maximum((KG.Xg.ravel() ** 2 - 1.0) * (1.0 - KG.Eg.ravel() ** 2), 0.0)
    for m in range(tri.MMAX + 1):
        now[f"CMP_rho01_m{m}"] = tri.coul_mode_potential((rho01 * rc2 ** (0.5 * m)).reshape(shape), m)
    for m in range(3):
        now[f"CMP_P00_m{m}"] = tri.coul_mode_potential(hVee.P00.reshape(shape), m)
        now[f"CMP_rho_m{m}"] = tri.coul_mode_potential(hVee.rho_g.reshape(shape), m)
    now["PSI_coul_dab"] = gVee.psi_coul(hVee.d_ab)
    now["PSI_yuk_dab"] = gVee.psi_yuk(hVee.d_ab, KG.GAM)
    now["PSI_yuk_daa"] = gVee.psi_yuk(hVee.d_aa, KG.GAM)
    Kc = hVee.make_kernel_coul()
    now["W_coul"] = Kc["W"]
    now["Psi_coul_ab"] = Kc["Psi"]["ab"]
    t0 = time.time()
    now["stage1"] = np.array([EN.V_aaaa, EN.V_bbbb, EN.V_aabb, EN.V_aaab, EN.V_bbab])
    now["Psi_aa"] = np.asarray(EN.Psi_aa)
    now["Psi_bb"] = np.asarray(EN.Psi_bb)
    ns: dict = {}
    exec(f"from {pkg}.lih_r12ci.energy import V_aaaa as _v", ns)        # PEP 562 from-import path
    assert ns["_v"] == EN.V_aaaa
    print(f"[{pkg}] Stage-1 attribute access {time.time() - t0:.1f} s (deferred build if > 0); "
          f"from-import resolves", flush=True)
    np.savez(out, **now)


def main() -> int:
    tmp = tempfile.mkdtemp(prefix="lih_headcheck_")
    try:
        raw = subprocess.run(["git", "archive", "--format=tar", "HEAD", "geovac/lih_r12ci"],
                             cwd=ROOT, check=True, capture_output=True).stdout
        with tarfile.open(fileobj=io.BytesIO(raw)) as tf:
            tf.extractall(tmp)
        os.makedirs(os.path.join(tmp, "headpkg"), exist_ok=True)
        shutil.move(os.path.join(tmp, "geovac", "lih_r12ci"), os.path.join(tmp, "headpkg", "lih_r12ci"))
        open(os.path.join(tmp, "headpkg", "__init__.py"), "w").close()
        head_sha = subprocess.run(["git", "rev-parse", "--short", "HEAD"], cwd=ROOT, check=True,
                                  capture_output=True, text=True).stdout.strip()
        print(f"HEAD = {head_sha}; committed geovac/lih_r12ci extracted to {tmp}\\headpkg", flush=True)
        outs = {}
        for pkg, root in (("headpkg", tmp), ("geovac", ROOT)):
            out = os.path.join(tmp, f"dump_{pkg}.npz")
            t0 = time.time()
            rc = subprocess.call([sys.executable, os.path.abspath(__file__), "--dump", pkg, root, out], cwd=ROOT)
            print(f"dump {pkg}: exit {rc} [{time.time() - t0:.0f} s]", flush=True)
            if rc != 0:
                return 2
            outs[pkg] = np.load(out)
        head, cur = outs["headpkg"], outs["geovac"]
        prist = np.load(PRISTINE, allow_pickle=True) if os.path.exists(PRISTINE) else None
        worst_cur = 0.0; worst_pr = 0.0
        print(f"\n{'array':18s} {'tree vs HEAD':>13s} {'pristine vs HEAD':>17s}")
        for k in head.files:
            d = float(np.max(np.abs(cur[k] - head[k]))); worst_cur = max(worst_cur, d)
            if prist is not None and k in prist.files:
                dp = float(np.max(np.abs(prist[k] - head[k]))); worst_pr = max(worst_pr, dp); ps = f"{dp:17.1e}"
            else:
                ps = f"{'(absent)':>17s}"
            print(f"{k:18s} {d:13.1e} {ps}")
        print(f"\nRESULT working-tree legacy path vs HEAD {head_sha}: max|Delta| over {len(head.files)} arrays = "
              f"{worst_cur:.1e} -> {'BIT-IDENTICAL' if worst_cur == 0.0 else 'DIFFERS'}")
        if prist is not None:
            print(f"RESULT pristine capture vs HEAD: max|Delta| = {worst_pr:.1e} -> "
                  f"{'produced from HEAD code' if worst_pr == 0.0 else 'NOT from HEAD code'}")
        return 0 if worst_cur == 0.0 else 1
    finally:
        shutil.rmtree(tmp, ignore_errors=True)


if __name__ == "__main__":
    if len(sys.argv) > 1 and sys.argv[1] == "--dump":
        dump(sys.argv[2], sys.argv[3], sys.argv[4])
        sys.exit(0)
    sys.exit(main())
