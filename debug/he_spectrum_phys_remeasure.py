"""Re-measure the Paper 34 He excited-state catalogue entries on the PHYSICAL
Hamiltonian (no (1/16)*A one-body hop).

Retraction context (2026-09-25): geovac.casimir_ci._build_graph_h1 adds a
(1/16)*A off-diagonal hop that is NOT part of the atom's Hamiltonian (in the
k=Z hydrogenic basis T+V_ne is exactly diagonal).  The Paper 34 He
2^1S-2^3S splitting (+11.86%) and 2^1P->1^1S oscillator strength (+61%) were
computed on that hop-augmented model.  This driver:

  1. reproduces the published hop-model numbers (control), then
  2. re-runs the identical drivers with include_hop=False (physical).

Drivers reused unchanged:
  debug/precision_catalogue_he_2s_singlet_triplet.py  (ss-only sub-blocks)
  debug/calc_track_he_oscillator_v1.py  (resurrected from git f8f6884)
The hop is switched off by patching the module-level _build_graph_h1 name.

Single-threaded.  Log: debug/data/residual_sweep.log (appended).
"""
from __future__ import annotations

import os
for _v in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS",
           "NUMEXPR_NUM_THREADS"):
    os.environ[_v] = "1"

import json
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "debug"))

import geovac.casimir_ci as cc  # noqa: E402
import precision_catalogue_he_2s_singlet_triplet as st  # noqa: E402
import calc_track_he_oscillator_v1 as osc  # noqa: E402

_ORIG = cc._build_graph_h1


def _phys(Z, n_max, include_hop=True):  # noqa: ANN001
    return _ORIG(Z, n_max, include_hop=False)


LOG = ROOT / "debug" / "data" / "residual_sweep.log"


def log(msg: str) -> None:
    print(msg, flush=True)
    with open(LOG, "a", encoding="utf-8") as fp:
        fp.write(msg + "\n")


def set_mode(mode: str) -> None:
    fn = _ORIG if mode == "hop" else _phys
    cc._build_graph_h1 = fn          # used by the splitting driver (local import)
    osc._build_graph_h1 = fn         # bound name in the oscillator driver


def main() -> None:
    out = {"splitting": {}, "oscillator": {}}
    log("\n=== he_spectrum_phys_remeasure.py  (%s) ===" % time.strftime("%Y-%m-%d %H:%M"))
    # Budget (single thread, cores shared): the hop-model splitting control
    # was run n=2..7 in a first pass and reproduced the published
    # debug/data/precision_catalogue_he_2s_singlet_triplet.json rows to all
    # printed digits (46651.0/9737.6/8125.7/7308.8/7250.3/7192.2 cm-1);
    # here hop = oscillator control only, splitting n<=8 phys.
    # Oscillator: the hop-model control at n=5 reproduced the published
    # f=0.445 (0.4458, +61.4%) in a second pass; n>=6 costs >10 min per
    # point single-threaded (sympy radial dipoles), so phys stops at n=5.
    split_n = {"hop": (), "phys": (2, 3, 4, 5, 6, 7, 8)}
    osc_n = {"hop": (), "phys": (2, 3, 4, 5)}
    for mode in ("hop", "phys"):
        set_mode(mode)
        log(f"--- mode={mode}: 2^1S-2^3S splitting (ss-only sub-block) ---")
        rows = []
        for n in split_n[mode]:
            t0 = time.time()
            r = st.compute_st_at_nmax(n)
            rows.append({k: r[k] for k in ("n_max", "E_1_1S_Ha", "E_2_3S_Ha",
                                           "E_2_1S_Ha", "dE_st_cm_inverse",
                                           "err_vs_NIST_pct")})
            log(f"  n={n:2d} E11S={r['E_1_1S_Ha']:.6f} E23S={r['E_2_3S_Ha']:.6f} "
                f"(Drake {st.DRAKE_E_2_3S_HA}) E21S={r['E_2_1S_Ha']:.6f} "
                f"(Drake {st.DRAKE_E_2_1S_HA}) dE={r['dE_st_cm_inverse']:.1f} cm-1 "
                f"err_NIST={r['err_vs_NIST_pct']:+.2f}%  t={time.time()-t0:.1f}s")
        out["splitting"][mode] = rows
        log(f"--- mode={mode}: 2^1P->1^1S oscillator strength (k=Z=2) ---")
        orows = []
        for n in osc_n[mode]:
            t0 = time.time()
            r = osc.compute_oscillator_strength(n_max=n, Z=2, k_orb=2.0)
            err = (r["f_length"] - osc.F_DRAKE_2P_1S) / osc.F_DRAKE_2P_1S * 100.0
            r["error_pct"] = err
            orows.append({k: float(r[k]) for k in ("n_max", "n_configs_2P",
                                                   "E_1S_Ha", "E_2P_Ha",
                                                   "omega_Ha", "dipole_squared_au",
                                                   "f_length", "error_pct")})
            log(f"  n={n} cfg2P={r['n_configs_2P']} E1S={r['E_1S_Ha']:.6f} "
                f"E2P={r['E_2P_Ha']:.6f} (Drake -2.123843) omega={r['omega_Ha']:.4f} "
                f"|z|^2={r['dipole_squared_au']:.4f} f={r['f_length']:.4f} "
                f"err={err:+.1f}%  t={time.time()-t0:.1f}s")
        out["oscillator"][mode] = orows
    set_mode("hop")
    p = ROOT / "debug" / "data" / "he_spectrum_phys_remeasure.json"
    p.write_text(json.dumps(out, indent=2))
    log(f"saved {p}")


if __name__ == "__main__":
    main()
