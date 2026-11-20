"""NOCI sandbox probe, step 1: H2 validation of non-orthogonal CI machinery.

Sandbox exploration (branch sandbox/noci, repo frozen at v4.76.0 -- no release path).
Scoping step for the named open follow-on in debug/sprint_commutator_and_explorer_memo.md:
"NOCI/VB-on-GeoVac feasibility ... scope before any build."

Question 1 here is pure machinery validation: does a tiny non-orthogonal config basis
(covalent + ionic, atom-centered 1s orbitals) bind H2 correctly, and how does the same
config count behave after orthogonalization?  Question 2 is the measurement ledger:
what would a device have to measure (S_IJ, H_IJ counts) in an NOQE-style scheme.

Basis note: orbitals are STO-nG Gaussian fits of 1s exponentials so every integral is
a rock-solid closed form (s-type Gaussians, Boys F0).  This is deliberately NOT a
GeoVac-native integral set -- step 1 validates the NOCI machinery only; the GeoVac
sparsity/native-integral confrontation is step 2 (LiH).

Anchors: Heitler-London with true 1s STOs at zeta=1 gives D_e ~ 3.14 eV, R_eq ~ 1.64 a0.
H-atom self-check: E(1 center, zeta=1) must sit just above -0.5 Ha.
"""

from __future__ import annotations

import json
import math
from typing import Dict, List, Tuple

import numpy as np
from scipy.linalg import eigh
from scipy.special import erf

HARTREE_TO_EV = 27.211386

# STO-nG expansions of a zeta=1 1s Slater function (alpha scales as zeta^2).
STO6G_1S = (
    np.array([23.31030, 4.235916, 1.185057, 0.4070989, 0.1580884, 0.06510954]),
    np.array([0.00916360, 0.04936150, 0.16853830, 0.37056280, 0.41649150, 0.13033400]),
)
STO3G_1S = (
    np.array([3.42525091, 0.62391373, 0.16885540]),
    np.array([0.15432897, 0.53532814, 0.44463454]),
)


def boys_f0(x: np.ndarray) -> np.ndarray:
    x = np.asarray(x, dtype=float)
    out = np.empty_like(x)
    small = x < 1e-12
    out[small] = 1.0 - x[small] / 3.0
    xs = x[~small]
    out[~small] = 0.5 * np.sqrt(np.pi / xs) * erf(np.sqrt(xs))
    return out


class Basis1s:
    """Contracted s-type Gaussian fit of a 1s Slater orbital at a given center."""

    def __init__(self, center: np.ndarray, zeta: float, alphas: np.ndarray, ds: np.ndarray):
        self.center = np.asarray(center, dtype=float)
        self.alphas = alphas * zeta**2
        norms = (2.0 * self.alphas / np.pi) ** 0.75
        self.coeffs = ds * norms
        al, bl = np.meshgrid(self.alphas, self.alphas, indexing="ij")
        self_ov = float(np.einsum("i,j,ij->", self.coeffs, self.coeffs, (np.pi / (al + bl)) ** 1.5))
        self.coeffs = self.coeffs / math.sqrt(self_ov)


def overlap(a: Basis1s, b: Basis1s) -> float:
    ra, rb = a.center, b.center
    r2 = float(np.dot(ra - rb, ra - rb))
    al, bl = np.meshgrid(a.alphas, b.alphas, indexing="ij")
    p = al + bl
    pref = (np.pi / p) ** 1.5 * np.exp(-al * bl / p * r2)
    return float(np.einsum("i,j,ij->", a.coeffs, b.coeffs, pref))


def kinetic(a: Basis1s, b: Basis1s) -> float:
    ra, rb = a.center, b.center
    r2 = float(np.dot(ra - rb, ra - rb))
    al, bl = np.meshgrid(a.alphas, b.alphas, indexing="ij")
    p = al + bl
    mu = al * bl / p
    s = (np.pi / p) ** 1.5 * np.exp(-mu * r2)
    t = mu * (3.0 - 2.0 * mu * r2) * s
    return float(np.einsum("i,j,ij->", a.coeffs, b.coeffs, t))


def nuclear(a: Basis1s, b: Basis1s, nuc: np.ndarray, z: float) -> float:
    ra, rb = a.center, b.center
    r2 = float(np.dot(ra - rb, ra - rb))
    val = 0.0
    for ai, ci in zip(a.alphas, a.coeffs):
        for bj, cj in zip(b.alphas, b.coeffs):
            p = ai + bj
            pc = (ai * ra + bj * rb) / p
            d2 = float(np.dot(pc - nuc, pc - nuc))
            v = -z * (2.0 * np.pi / p) * math.exp(-ai * bj / p * r2) * float(boys_f0(np.array([p * d2]))[0])
            val += ci * cj * v
    return val


def eri(a: Basis1s, b: Basis1s, c: Basis1s, d: Basis1s) -> float:
    """Chemist notation (ab|cd) = int a(1)b(1) 1/r12 c(2)d(2)."""
    rab2 = float(np.dot(a.center - b.center, a.center - b.center))
    rcd2 = float(np.dot(c.center - d.center, c.center - d.center))
    val = 0.0
    for ai, ci in zip(a.alphas, a.coeffs):
        for bj, cj in zip(b.alphas, b.coeffs):
            p = ai + bj
            rp = (ai * a.center + bj * b.center) / p
            ke_ab = math.exp(-ai * bj / p * rab2)
            for ck, cck in zip(c.alphas, c.coeffs):
                for dl, cdl in zip(d.alphas, d.coeffs):
                    q = ck + dl
                    rq = (ck * c.center + dl * d.center) / q
                    ke_cd = math.exp(-ck * dl / q * rcd2)
                    d2 = float(np.dot(rp - rq, rp - rq))
                    pref = 2.0 * np.pi**2.5 / (p * q * math.sqrt(p + q))
                    f0 = float(boys_f0(np.array([p * q / (p + q) * d2]))[0])
                    val += ci * cj * cck * cdl * pref * ke_ab * ke_cd * f0
    return val


def integral_set(orbs: List[Basis1s], nuclei: List[Tuple[np.ndarray, float]]):
    n = len(orbs)
    s = np.zeros((n, n))
    h = np.zeros((n, n))
    for i in range(n):
        for j in range(n):
            s[i, j] = overlap(orbs[i], orbs[j])
            h[i, j] = kinetic(orbs[i], orbs[j]) + sum(
                nuclear(orbs[i], orbs[j], pos, z) for pos, z in nuclei
            )
    g = np.zeros((n, n, n, n))
    for i in range(n):
        for j in range(n):
            for k in range(n):
                for l in range(n):
                    g[i, j, k, l] = eri(orbs[i], orbs[j], orbs[k], orbs[l])
    return s, h, g


def config_matrices(configs: List[Tuple[int, int]], s: np.ndarray, h: np.ndarray, g: np.ndarray):
    """Singlet spatial pair functions Psi_I ~ i1(1) i2(2) + i2(1) i1(2) (unnormalized).

    Returns generalized (Hmat, Smat) over the config list.  Chemist eri g[p,r,q,s]
    = (pr|qs) pairs electron 1 with (p,r), electron 2 with (q,s).
    """
    nc = len(configs)
    smat = np.zeros((nc, nc))
    hmat = np.zeros((nc, nc))
    for a_idx, (i1, i2) in enumerate(configs):
        for b_idx, (j1, j2) in enumerate(configs):
            s_ab = 0.0
            h_ab = 0.0
            for p, q in ((i1, i2), (i2, i1)):
                for r, t in ((j1, j2), (j2, j1)):
                    s_ab += s[p, r] * s[q, t]
                    h_ab += h[p, r] * s[q, t] + s[p, r] * h[q, t]
                    h_ab += g[p, r, q, t]
            smat[a_idx, b_idx] = s_ab
            hmat[a_idx, b_idx] = h_ab
    return hmat, smat


def gen_eig_ground(hmat: np.ndarray, smat: np.ndarray, tol: float = 1e-10) -> Tuple[float, float]:
    """Ground eigenvalue via canonical orthogonalization; returns (E, cond(S))."""
    w, v = np.linalg.eigh(smat)
    keep = w > tol * w.max()
    x = v[:, keep] / np.sqrt(w[keep])
    hp = x.T @ hmat @ x
    e = np.linalg.eigvalsh(hp)[0]
    cond = w.max() / w.min() if w.min() > 0 else np.inf
    return float(e), float(cond)


def h2_energies(r: float, zeta: float, alphas: np.ndarray, ds: np.ndarray) -> Dict[str, float]:
    pos_a = np.array([0.0, 0.0, 0.0])
    pos_b = np.array([0.0, 0.0, r])
    nuclei = [(pos_a, 1.0), (pos_b, 1.0)]
    chi_a = Basis1s(pos_a, zeta, alphas, ds)
    chi_b = Basis1s(pos_b, zeta, alphas, ds)
    s, h, g = integral_set([chi_a, chi_b], nuclei)
    vnn = 1.0 / r

    cov = (0, 1)
    ion_a = (0, 0)
    ion_b = (1, 1)

    out: Dict[str, float] = {"R": r, "S_ab": s[0, 1]}

    # (1) covalent-only NOCI (Heitler-London)
    hm, sm = config_matrices([cov], s, h, g)
    out["E_cov_only"] = hm[0, 0] / sm[0, 0] + vnn

    # (2) full non-orthogonal 3x3 (covalent + both ionic)
    hm, sm = config_matrices([cov, ion_a, ion_b], s, h, g)
    e3, cond3 = gen_eig_ground(hm, sm)
    out["E_noci3"] = e3 + vnn
    out["cond_S_config"] = cond3

    # (3) covalent + symmetric-ionic (2 effective configs; ionic combo forced by g-symmetry)
    hm2, sm2 = config_matrices([cov, ion_a], s, h, g)  # asymmetric truncation, diagnostic only
    e2, _ = gen_eig_ground(hm2, sm2)
    out["E_noci2_asym"] = e2 + vnn

    # (4) orthogonal-side comparators built from Loewdin-orthogonalized orbitals
    w, v = np.linalg.eigh(s)
    s_invhalf = v @ np.diag(w**-0.5) @ v.T
    ht = s_invhalf @ h @ s_invhalf
    gt = np.einsum("pi,qj,rk,sl,ijkl->pqrs", s_invhalf, s_invhalf, s_invhalf, s_invhalf, g)
    st = np.eye(2)
    hm, sm = config_matrices([cov], st, ht, gt)
    out["E_lowdin_cov_only"] = hm[0, 0] / sm[0, 0] + vnn
    hm, sm = config_matrices([cov, ion_a, ion_b], st, ht, gt)
    e3t, _ = gen_eig_ground(hm, sm)
    out["E_lowdin3"] = e3t + vnn  # span identity: must equal E_noci3

    # (5) MO single determinant (RHF-like sigma_g^2)
    ng = 1.0 / math.sqrt(2.0 + 2.0 * s[0, 1])
    cg = np.array([ng, ng])
    hgg = cg @ h @ cg
    gggg = np.einsum("i,j,k,l,ijkl->", cg, cg, cg, cg, g)
    out["E_mo_det"] = 2.0 * hgg + gggg + vnn
    return out


def main() -> None:
    # H-atom self-check selects the basis fit quality.
    chosen = None
    for name, (alphas, ds) in (("STO-6G", STO6G_1S), ("STO-3G", STO3G_1S)):
        pos = np.array([0.0, 0.0, 0.0])
        chi = Basis1s(pos, 1.0, alphas, ds)
        s, h, _ = integral_set([chi], [(pos, 1.0)])
        norm = s[0, 0]
        e_atom = h[0, 0] / norm
        print(f"[self-check] {name}: <chi|chi> = {norm:.8f}, E(H atom, zeta=1) = {e_atom:.6f} Ha")
        ok = abs(norm - 1.0) < 1e-6 and -0.5 <= e_atom < -0.497
        if ok and chosen is None:
            chosen = (name, alphas, ds, e_atom)
    if chosen is None:
        name, alphas, ds = "STO-3G", *STO3G_1S
        pos = np.array([0.0, 0.0, 0.0])
        chi = Basis1s(pos, 1.0, alphas, ds)
        s, h, _ = integral_set([chi], [(pos, 1.0)])
        e_atom = h[0, 0] / s[0, 0]
        chosen = (name, alphas, ds, e_atom)
        print(f"[self-check] falling back to {name} (E_atom = {e_atom:.6f})")
    name, alphas, ds, e_atom = chosen
    print(f"[basis] using {name}; dissociation reference 2*E_atom = {2 * e_atom:.6f} Ha\n")

    zeta = 1.0
    grid = [round(r, 2) for r in np.concatenate([np.arange(0.8, 3.01, 0.1), [4.0, 5.0, 6.0, 8.0]])]
    rows = [h2_energies(r, zeta, alphas, ds) for r in grid]

    hdr = ("R", "S_ab", "E_cov_only", "E_noci3", "E_lowdin_cov", "E_mo_det", "span_diff")
    print("{:>5} {:>8} {:>12} {:>12} {:>13} {:>12} {:>10}".format(*hdr))
    for row in rows:
        print(
            "{:>5.2f} {:>8.4f} {:>12.6f} {:>12.6f} {:>13.6f} {:>12.6f} {:>10.2e}".format(
                row["R"], row["S_ab"], row["E_cov_only"], row["E_noci3"],
                row["E_lowdin_cov_only"], row["E_mo_det"],
                abs(row["E_noci3"] - row["E_lowdin3"]),
            )
        )

    e_inf = 2.0 * e_atom

    def well(key: str) -> Tuple[float, float, float]:
        es = np.array([row[key] for row in rows])
        idx = int(np.argmin(es))
        r_eq, e_min = rows[idx]["R"], float(es[idx])
        if 0 < idx < len(rows) - 1:
            x = np.array([rows[idx - 1]["R"], rows[idx]["R"], rows[idx + 1]["R"]])
            y = np.array([es[idx - 1], es[idx], es[idx + 1]])
            c = np.polyfit(x, y, 2)
            r_eq = float(-c[1] / (2.0 * c[0]))
            e_min = float(np.polyval(c, r_eq))
        return r_eq, e_min, (e_inf - e_min)

    print("\n[wells]  (D_e vs 2*E_atom of this basis; anchors: HL-STO zeta=1 D_e~3.14 eV, R_eq~1.64 a0)")
    summary = {}
    for key, label in (
        ("E_cov_only", "covalent-only NOCI (Heitler-London)"),
        ("E_noci3", "3-config NOCI (cov + 2 ionic)"),
        ("E_lowdin_cov_only", "Loewdin-orbital covalent-only"),
        ("E_mo_det", "MO single determinant (RHF-like)"),
    ):
        r_eq, e_min, d_e = well(key)
        bound = d_e > 0
        print(
            f"  {label:38s} R_eq = {r_eq:6.3f} a0   E_min = {e_min:10.6f} Ha   "
            f"D_e = {d_e:8.5f} Ha = {d_e * HARTREE_TO_EV:6.3f} eV   {'BINDS' if bound else 'UNBOUND'}"
        )
        summary[key] = {"R_eq": r_eq, "E_min": e_min, "D_e_Ha": d_e, "binds": bool(bound)}

    span_max = max(abs(r["E_noci3"] - r["E_lowdin3"]) for r in rows)
    cond_max = max(r["cond_S_config"] for r in rows)
    print(f"\n[checks] span identity max|E_noci3 - E_lowdin3| = {span_max:.2e} (must be ~0)")
    print(f"[checks] worst config-overlap condition number   = {cond_max:.2e}")

    n_conf = 3
    print("\n[measurement ledger, NOQE-style, per PES point]")
    print(f"  configs: {n_conf}  ->  S_IJ to measure: {n_conf * (n_conf + 1) // 2}"
          f"   H_IJ to measure: {n_conf * (n_conf + 1) // 2}")
    print("  each H_IJ = sum over the NATIVE (untransformed) Pauli decomposition of H;")
    print("  no S^-1/2 congruence anywhere -> integral sparsity untouched; cost relocated")
    print("  to state-overlap circuits (Hadamard/SWAP-type) + generalized-eig conditioning.")

    out = {
        "basis": name,
        "zeta": zeta,
        "E_atom": e_atom,
        "rows": rows,
        "wells": summary,
        "span_identity_max": span_max,
        "cond_S_config_max": cond_max,
    }
    with open("debug/data/noci_h2_probe_results.json", "w") as fh:
        json.dump(out, fh, indent=1)
    print("\n[saved] debug/data/noci_h2_probe_results.json")


if __name__ == "__main__":
    main()
