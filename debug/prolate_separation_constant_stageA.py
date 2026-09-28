"""Stage A driver: the prolate separation constant Lambda on one-electron two-centre systems.

Systems: H2+ (Z_A = Z_B = 1, eta-parity blocks g/u) and (Z_A, Z_B) = (3, 1).
For m = 0, 1, 2 and R on a log grid [0.01, 80]:
  1. ||[H, Lambda]|| (W-orthonormal frame, relative) + the sign-flipped-Z_- control;
  2. joint eigenpairs, labelled by the exact nodal numbers (n_xi, n_eta) (Klein oscillation
     theorem on the separated 1D problems) -- the labels are R-invariant, so they TRACK states;
  3. united-atom limit: lam -> l(l+1), E -> -Z_UA^2 / 2 n^2 with l = n_eta + m,
     n = n_xi + n_eta + m + 1;
  4. separated-atom limit: fit lam + p^2 = s R + c0 + c1/R (p^2 = -R^2 E / 2) on the three
     largest R; s - Z_other = (Z_self/n) k must land on the spectrum of the SO(4) Runge-Lenz
     A_z restricted to shell n, fixed m (built from tests/test_graph_s3_nonconvergence.py);
     cross-checked against the first-order Stark shift;
  5. same-symmetry crossings (true crossings of tracks with different Lambda).
Output: debug/data/prolate_separation_constant_stageA.json + stdout log.
"""
from __future__ import annotations

import json
import os
import sys
import time
from typing import Dict, List, Tuple

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.dirname(HERE))

from prolate_separation_constant import OneE, joint_labeled, restrict  # noqa: E402
sys.path.insert(0, os.path.join(os.path.dirname(HERE), "tests"))
from test_graph_s3_nonconvergence import so4_generators  # noqa: E402

NX, NE, K = 30, 40, 10
RGRID = np.round(np.logspace(np.log10(0.01), np.log10(80.0), 90), 6)


def alpha_of(R: float) -> float:
    return max(0.005, 0.5 * R)


def az_spectrum(n: int, m: int) -> List[int]:
    """Integer spectrum of the energy-scaled Runge-Lenz A_z on shell n, fixed m."""
    _, As, states = so4_generators(n)
    Az = As[2].real
    idx = [i for i, (nn, l, mm) in enumerate(states) if nn == n and mm == m]
    ev = np.linalg.eigvalsh(Az[np.ix_(idx, idx)])
    return sorted(int(round(e)) for e in ev), float(np.abs(ev - np.round(ev)).max())


def blocks(ZA: float, ZB: float) -> List[Tuple[str, int]]:
    return [("g", 0), ("u", 1)] if ZA == ZB else [("het", -1)]


def run_system(ZA: float, ZB: float, log) -> Dict:
    tracks: Dict[str, Dict] = {}
    comm = []
    for m in (0, 1, 2):
        for R in RGRID:
            al = alpha_of(R)
            op = OneE(R, ZA, ZB, m, NX, NE, al)
            if R in (0.01, RGRID[30], RGRID[60], RGRID[-1]) or abs(R - 1.0) < 0.06:
                c = op.commutator()
                cf = OneE(R, ZA, ZB, m, NX, NE, al, flip_minus=True).commutator()
                comm.append(dict(m=m, R=float(R), comm=c, comm_flip=cf))
                log(f"  Z=({ZA},{ZB}) m={m} R={R:8.4f}  ||[H,L]||rel={c:.2e}  flipZ-={cf:.2e}")
            for name, par in blocks(ZA, ZB):
                o = restrict(op, par) if par >= 0 else op
                for st in joint_labeled(o, K, al):
                    key = f"m{m}_{name}_x{st['nxi']}_e{st['neta']}"
                    tr = tracks.setdefault(key, dict(m=m, blk=name, nxi=st["nxi"],
                                                     neta=st["neta"], R=[], E=[], lam=[],
                                                     eta=[], rank1=[], sep=[]))
                    if tr["R"] and tr["R"][-1] == float(R):
                        tr["dup"] = True          # label collision (should never happen)
                        continue
                    tr["R"].append(float(R)); tr["E"].append(st["E"]); tr["lam"].append(st["lam"])
                    tr["eta"].append(st["eta"]); tr["rank1"].append(st["rank1"])
                    tr["sep"].append(st["sep_dist"])
    return dict(tracks=tracks, comm=comm)


def analyse(ZA: float, ZB: float, res: Dict, log) -> Dict:
    ZUA = ZA + ZB
    rows = []
    for key, tr in sorted(res["tracks"].items(),
                          key=lambda kv: (kv[1]["m"], kv[1]["E"][0])):
        R = np.array(tr["R"]); E = np.array(tr["E"]); lam = np.array(tr["lam"])
        if R[0] > 0.011 or R[-1] < 79.0 or len(R) < len(RGRID):
            continue                       # only states present over the whole range
        m, nx, ne = tr["m"], tr["nxi"], tr["neta"]
        l = ne + m
        n_ua = nx + ne + m + 1
        ua_lam_err = [float(lam[i] - l * (l + 1)) for i in range(3)]
        ua_E_err = [float(E[i] + ZUA ** 2 / (2 * n_ua ** 2)) for i in range(3)]
        # separated atom
        Rf, Ef, lf = R[-3:], E[-3:], lam[-3:]
        y = lf - Rf ** 2 * Ef / 2.0
        A = np.vstack([Rf, np.ones(3), 1.0 / Rf]).T
        s, c0, c1 = np.linalg.solve(A, y)
        eta_end = tr["eta"][-1]
        if ZA == ZB:
            atom, Zs, Zo = "A|B", ZA, ZB
        else:
            atom, Zs, Zo = ("A", ZA, ZB) if eta_end < 0 else ("B", ZB, ZA)
        Eat = E[-1] + Zo / R[-1]
        n_sa = int(round(Zs / np.sqrt(-2 * Eat)))
        mu = s - Zo
        k_lam = mu * n_sa / Zs
        spec, spec_err = az_spectrum(n_sa, m) if n_sa > m else ([], 0.0)
        k_int = int(round(k_lam))
        # first-order Stark: E + Zo/R + Zs^2/2n^2 ~ -(3 n k_E / 2 Zs) Zo / R^2 (fit 2 pts)
        dE = E[-2:] + Zo / R[-2:] + Zs ** 2 / (2 * n_sa ** 2)
        k_E = float(-np.mean(dE * R[-2:] ** 2) / (1.5 * n_sa / Zs * Zo))
        row = dict(key=key, m=m, blk=tr["blk"], nxi=nx, neta=ne,
                   UA=f"{n_ua}{'spdfghiklmnoq'[l]}", UA_lam_err=ua_lam_err, UA_E_err=ua_E_err,
                   SA_atom=atom, SA_n=n_sa, s=float(s), c0=float(c0), mu=float(mu),
                   k_lambda=float(k_lam), k_int=k_int, k_in_Az_spectrum=k_int in spec,
                   Az_spectrum=spec, k_lambda_dev=float(abs(k_lam - k_int)), k_stark=k_E,
                   max_rank1=float(max(tr["rank1"])), max_sep=float(max(tr["sep"])),
                   dup=bool(tr.get("dup", False)))
        rows.append(row)
        log(f"  {key:22s} UA {row['UA']:4s} dlam(R=.01)={ua_lam_err[0]:+.2e} "
            f"dE={ua_E_err[0]:+.1e} | SA atom {atom} n={n_sa} k_Lam={k_lam:+.4f} "
            f"(A_z spec {spec}) k_Stark={k_E:+.3f} | rank1<={row['max_rank1']:.0e} "
            f"sep<={row['max_sep']:.0e}{' DUP' if row['dup'] else ''}")
    return dict(rows=rows)


def crossings(ZA: float, ZB: float, res: Dict, log) -> List[Dict]:
    out = []
    tr = res["tracks"]
    by_blk: Dict[Tuple[int, str], List[str]] = {}
    for k, t in tr.items():
        if len(t["R"]) == len(RGRID):
            by_blk.setdefault((t["m"], t["blk"]), []).append(k)
    for (m, blk), keys in sorted(by_blk.items()):
        keys = sorted(keys, key=lambda k: tr[k]["E"][0])[:8]
        for i in range(len(keys)):
            for j in range(i + 1, len(keys)):
                a, b = tr[keys[i]], tr[keys[j]]
                d = np.array(a["E"]) - np.array(b["E"])
                flips = np.where(np.sign(d[1:]) != np.sign(d[:-1]))[0]
                for f in flips:
                    R1, R2 = RGRID[f], RGRID[f + 1]
                    # linear interpolation of the crossing point + Lambda there
                    t_ = d[f] / (d[f] - d[f + 1])
                    Rx = R1 + t_ * (R2 - R1)
                    la = a["lam"][f] + t_ * (a["lam"][f + 1] - a["lam"][f])
                    lb = b["lam"][f] + t_ * (b["lam"][f + 1] - b["lam"][f])
                    Ex = a["E"][f] + t_ * (a["E"][f + 1] - a["E"][f])
                    out.append(dict(m=m, blk=blk, a=keys[i], b=keys[j], R=float(Rx),
                                    E=float(Ex), lam_a=float(la), lam_b=float(lb)))
                    log(f"  crossing m={m} {blk}: {keys[i]} x {keys[j]} at R~{Rx:.3f} "
                        f"E~{Ex:.5f}  Lambda {la:.3f} vs {lb:.3f}")
    return out


def main() -> None:
    os.makedirs(os.path.join(HERE, "data"), exist_ok=True)
    lines: List[str] = []

    def log(s: str) -> None:
        print(s, flush=True)
        lines.append(s)

    out = {}
    t0 = time.time()
    for (ZA, ZB) in [(1.0, 1.0), (3.0, 1.0)]:
        log(f"=== Z_A={ZA} Z_B={ZB} ===")
        res = run_system(ZA, ZB, log)
        an = analyse(ZA, ZB, res, log)
        cr = crossings(ZA, ZB, res, log)
        # basis check at the largest R (states tracked above)
        R = float(RGRID[-1])
        e1 = OneE(R, ZA, ZB, 0, NX, NE, alpha_of(R)).joint(8)["E"]
        e2 = OneE(R, ZA, ZB, 0, 40, 70, alpha_of(R)).joint(8)["E"]
        log(f"  basis check R={R}: max|dE| (30x40 vs 40x70) = {np.abs(e1 - e2).max():.1e}")
        out[f"{int(ZA)}_{int(ZB)}"] = dict(comm=res["comm"], rows=an["rows"], crossings=cr,
                                           basis_check=float(np.abs(e1 - e2).max()),
                                           tracks=res["tracks"])
        log(f"  [{time.time() - t0:.0f}s]")
    with open(os.path.join(HERE, "data", "prolate_separation_constant_stageA.json"), "w") as f:
        json.dump(out, f, indent=1)
    with open(os.path.join(HERE, "data", "prolate_separation_constant_stageA.log"), "w") as f:
        f.write("\n".join(lines))


if __name__ == "__main__":
    main()
