"""I/O ladder accounting — Rung 1 (the spine).

DIAGNOSTIC ONLY. Measures LOAD-vs-GENERATE: how many INDEPENDENT SCALARS must be
shipped to a device to specify a molecular electronic Hamiltonian, in two
encodings:

  (i)  GAUSSIAN LOAD  = the whole 8-fold-symmetry-reduced ERI tensor + h1 + const.
                        Every entry is a distinct classically-computed number;
                        a Gaussian tensor carries no on-device generation rule.
  (ii) GEOVAC LOAD    = radial-seed instances only. The ANGULAR structure of the
                        GeoVac Coulomb-Sturmian / hydrogenic Hamiltonian is
                        generated ON-DEVICE from (n,l,m) labels via Gaunt/3j
                        (Papers 22/58), so the m-degeneracy (2l+1 orbitals per
                        shell) is NOT loaded. Only the 1D radial reduced
                        integrals — indexed by RADIAL shells (n,l) per center,
                        not by orbitals (n,l,m) — must be tabulated and shipped.
  (iii) GEOVAC GENERATE = size of the label->angular RULE (distinct 3j/Gaunt
                        evaluations). A fixed algorithm; grows only with l_max,
                        independent of center count and of the m-multiplicity.

This is NOT tensor-nonzeros (N3b) and NOT Pauli-term count (QC-1). It is the
independent-scalar SHIP-IN count. A dense tensor GENERATED on-device from a
compact rule + a few radial seeds is still low-I/O.

NOTE ON GRANULARITY. This is Rung 1: the radial-seed count is taken at the
SHELL-QUARTET granularity (the m-collapsed analog of the Gaussian ERI tensor),
which isolates the one mechanism the thesis rests on — angular m-degeneracy is
generated, not loaded. A multipole-resolved refinement is reported as a
secondary column; Rung 2 sharpens the radial count further (distinct
exponent/argument instances of the closed-form {E1, ln, gamma} engine).

Run:  python debug/io_ladder_accounting.py
"""

from __future__ import annotations

import warnings
from math import log

import numpy as np

warnings.filterwarnings("ignore")

from geovac.composed_qubit import _enumerate_states  # noqa: E402


# ---------------------------------------------------------------------------
# Basis reconstruction: replicate build_composed_hamiltonian's sub-block loop
# ---------------------------------------------------------------------------

def subblocks_for_spec(spec):
    """Return [(label, Z, states)] exactly as build_composed_hamiltonian enumerates."""
    subs = []
    for blk in spec.blocks:
        l_min = getattr(blk, "l_min", 0)
        subs.append((blk.label + "_center", blk.Z_center,
                     _enumerate_states(blk.max_n, l_min=l_min)))
        if blk.has_h_partner:
            pm = blk.max_n_partner if blk.max_n_partner > 0 else blk.max_n
            subs.append((blk.label + "_partner", blk.Z_partner,
                         _enumerate_states(pm)))
    return subs


def basis_summary(subs):
    """From sub-blocks, extract global orbitals and global radial shells.

    A global orbital is (subblock_idx, n, l, m).
    A global radial shell is (subblock_idx, n, l): the radial factor depends on
    (n, l) and the center's exponent Z, NOT on m. So the 2l+1 orbitals of a
    shell share ONE radial seed family.
    """
    orbitals, shells = [], []
    shell_l = []  # l value of each global shell (for multipole counting)
    shell_sub = []  # subblock index of each global shell
    for si, (_lab, _Z, states) in enumerate(subs):
        seen = set()
        for (n, l, m) in states:
            orbitals.append((si, n, l, m))
            if (n, l) not in seen:
                seen.add((n, l))
                shells.append((si, n, l))
                shell_l.append(l)
                shell_sub.append(si)
    return orbitals, shells, shell_l, shell_sub


# ---------------------------------------------------------------------------
# Independent-scalar counters
# ---------------------------------------------------------------------------

def unique8(K: int) -> int:
    """Number of 8-fold-permutation-unique (pq|rs) index quartets for K functions.

    npair = K(K+1)/2 ;  n_unique = npair(npair+1)/2 . This is the count of
    distinct real-orbital two-electron integrals a FCIDUMP would carry.
    """
    npair = K * (K + 1) // 2
    return npair * (npair + 1) // 2


def gaussian_load(M: int) -> dict:
    """Independent scalars to ship a Gaussian-encoded Hamiltonian at M orbitals."""
    eri = unique8(M)
    h1 = M * (M + 1) // 2
    return {"eri": eri, "h1": h1, "const": 1, "total": eri + h1 + 1}


def _n_multipoles(la: int, lb: int) -> int:
    """# even-parity multipoles L in [|la-lb|, la+lb]: = min(la,lb)+1.

    (multipole_decomposition keeps L with (la+lb+L) even, |la-lb|<=L<=la+lb.)
    """
    return min(la, lb) + 1


def geovac_load(shells, shell_l, shell_sub) -> dict:
    """Radial-seed instances only (angular generated on-device => 0 loaded).

    Primary (shell-quartet, first cut): the m-collapsed analog of the Gaussian
    ERI tensor — 8-fold-unique quartets of RADIAL shells. One reduced radial
    integral family per quartet.

    Secondary (multipole-resolved proxy): Sum over 8-fold-unique shell quartets
    of n_mult(pair1) * n_mult(pair2), where n_mult(la,lb)=min(la,lb)+1 is the #
    of distinct multipole radial integrals the closed-form engine emits for a
    shell pair. A bounded O(1) l-dependent multiplier on the primary count.

    One-electron: the diagonal h1 = -Z^2/(2 n^2) is a CLOSED-FORM function of the
    labels (n, Z) — generated, 0 loaded. Only cross-center V_ne reduced integrals
    (weight-1 {E1, ln, gamma} closed forms, Paper 58) are shipped: one family per
    cross-subblock shell pair.
    """
    S = len(shells)
    # Primary: shell-quartet count (8-fold unique).
    eri_shellquartet = unique8(S)

    # Secondary: multipole-resolved proxy. Enumerate 8-fold-unique quartets.
    # pair index p=(i,j) with i>=j; quartet (p,q) with p>=q.
    pairs = []
    for i in range(S):
        for j in range(i + 1):
            pairs.append((i, j))
    nmult_pair = [_n_multipoles(shell_l[i], shell_l[j]) for (i, j) in pairs]
    eri_multipole = 0
    for pi in range(len(pairs)):
        for qi in range(pi + 1):
            eri_multipole += nmult_pair[pi] * nmult_pair[qi]

    # One-electron cross-center radial seeds: shell pairs on DIFFERENT subblocks.
    h1_cross = 0
    for i in range(S):
        for j in range(i):
            if shell_sub[i] != shell_sub[j]:
                h1_cross += 1

    return {
        "eri_shellquartet": eri_shellquartet,
        "eri_multipole": eri_multipole,
        "h1_cross": h1_cross,
        "const": 1,
        "total_primary": eri_shellquartet + h1_cross + 1,
        "total_multipole": eri_multipole + h1_cross + 1,
        "S": S,
    }


def generate_rule_size(l_max: int) -> dict:
    """Size of the label->angular RULE: distinct nonzero 3j/Gaunt evaluations.

    The angular tensor is generated from Gaunt coefficients
    <Y_LM | conj(Y_l1m1) Y_l2m2> = product of two Wigner-3j. We count the
    distinct NONZERO Gaunt tuples over the (l,m) label content up to l_max
    (l = n-1). This depends ONLY on l_max — independent of the number of centers
    and of the m-multiplicity of the loaded tensor. The algorithm itself
    (one Gaunt function + bounded multipole loop) is O(1) code, reused verbatim
    across every system.
    """
    from sympy.physics.wigner import gaunt

    seen = set()
    n3j = 0
    for l1 in range(l_max + 1):
        for l2 in range(l_max + 1):
            for L in range(abs(l1 - l2), l1 + l2 + 1):
                if (l1 + l2 + L) % 2 != 0:
                    continue
                for m1 in range(-l1, l1 + 1):
                    for m2 in range(-l2, l2 + 1):
                        M = m2 - m1
                        if abs(M) > L:
                            continue
                        g = gaunt(l1, l2, L, -m1, m2, -M)
                        if g != 0:
                            seen.add((l1, l2, L, m1, m2, M))
                            n3j += 1
    return {"distinct_gaunt": n3j, "l_max": l_max}


# ---------------------------------------------------------------------------
# System table
# ---------------------------------------------------------------------------

def build_specs(max_n=2):
    from geovac.molecular_spec import hydride_spec, MolecularSpec, OrbitalBlock

    he = MolecularSpec(
        name="He",
        blocks=[OrbitalBlock(label="He_core", block_type="atomic",
                             Z_center=2.0, n_electrons=2, max_n=max_n)],
        nuclear_repulsion_constant=0.0,
    )
    h2 = MolecularSpec(
        name="H2",
        blocks=[OrbitalBlock(label="H2_bond", block_type="bond_pair",
                             Z_center=1.0, n_electrons=2, max_n=max_n)],
        nuclear_repulsion_constant=1.0 / 1.4,
    )
    systems = [
        ("He", 1, he),
        ("H2", 2, h2),
        ("LiH", 2, hydride_spec(3, max_n=max_n)),
        ("BeH2", 3, hydride_spec(4, max_n=max_n)),
        ("H2O", 3, hydride_spec(8, max_n=max_n)),
        ("CH4", 5, hydride_spec(6, max_n=max_n)),
    ]
    return systems


def loglog_fit(xs, ys):
    """Return (exponent, prefactor) for y ~ a * x^b via least squares on logs."""
    lx = np.log(np.array(xs, float))
    ly = np.log(np.array(ys, float))
    b, la = np.polyfit(lx, ly, 1)
    return b, float(np.exp(la))


def main():
    print("=" * 92)
    print("I/O LADDER ACCOUNTING — Rung 1 (load vs generate)  [DIAGNOSTIC]")
    print("=" * 92)

    max_n = 2
    l_max = max_n - 1
    systems = build_specs(max_n=max_n)

    rows = []
    print(f"\nBasis: composed hydrogenic, max_n={max_n} (l_max={l_max}) per center.\n")
    hdr = (f"{'system':>6} {'atoms':>5} {'nsub':>4} {'M':>4} {'S':>4} "
           f"{'GAUSS load':>12} {'GEOVAC load':>12} {'GEOVAC(mult)':>12} "
           f"{'ratio G/Ga':>10}")
    print(hdr)
    print("-" * len(hdr))
    for name, natoms, spec in systems:
        subs = subblocks_for_spec(spec)
        orbitals, shells, shell_l, shell_sub = basis_summary(subs)
        M = len(orbitals)
        gl = gaussian_load(M)
        gv = geovac_load(shells, shell_l, shell_sub)
        ratio = gv["total_primary"] / gl["total"]
        rows.append(dict(name=name, natoms=natoms, nsub=len(subs), M=M,
                         S=gv["S"], gauss=gl["total"], gauss_eri=gl["eri"],
                         geovac=gv["total_primary"], geovac_mult=gv["total_multipole"],
                         geovac_eri=gv["eri_shellquartet"], h1_cross=gv["h1_cross"]))
        print(f"{name:>6} {natoms:>5} {len(subs):>4} {M:>4} {gv['S']:>4} "
              f"{gl['total']:>12,} {gv['total_primary']:>12,} "
              f"{gv['total_multipole']:>12,} {ratio:>10.4f}")

    # Generate-rule size (fixed; depends only on l_max)
    gr = generate_rule_size(l_max)
    print(f"\nGEOVAC GENERATE (rule size): {gr['distinct_gaunt']} distinct nonzero "
          f"Gaunt/3j evaluations at l_max={l_max}.")
    print("  -> This is a FIXED algorithm: same code + same 3j table for EVERY")
    print("     system, every center, every m. It does NOT scale with M or")
    print("     with the number of centers; it grows only with l_max.")

    # -------- scaling fits across the 6 systems (vs M) --------
    Ms = [r["M"] for r in rows]
    ga = [r["gauss"] for r in rows]
    gv = [r["geovac"] for r in rows]
    gvm = [r["geovac_mult"] for r in rows]
    bG, aG = loglog_fit(Ms, ga)
    bV, aV = loglog_fit(Ms, gv)
    bVm, aVm = loglog_fit(Ms, gvm)
    print("\n" + "=" * 92)
    print("SCALING vs basis size M  (log-log fit y ~ a*M^b, 6 systems, max_n=2)")
    print("=" * 92)
    print(f"  Gaussian LOAD          : b = {bG:.3f}   (a={aG:.3g})")
    print(f"  GeoVac  LOAD (shell)   : b = {bV:.3f}   (a={aV:.3g})")
    print(f"  GeoVac  LOAD (multipole): b = {bVm:.3f}  (a={aVm:.3g})")
    print(f"  ratio GeoVac/Gaussian at M={Ms[-1]}: {gv[-1]/ga[-1]:.4f} "
          f"(= constant-factor {ga[-1]/gv[-1]:.1f}x fewer scalars)")

    # -------- scaling vs center count (at fixed max_n) --------
    nsubs = [r["nsub"] for r in rows]
    bGc, _ = loglog_fit([max(n, 1) for n in nsubs], ga)
    bVc, _ = loglog_fit([max(n, 1) for n in nsubs], gv)
    print("\nSCALING vs # basis sub-blocks (centers), fixed max_n=2:")
    print(f"  Gaussian LOAD exponent in nsub : {bGc:.3f}")
    print(f"  GeoVac   LOAD exponent in nsub : {bVc:.3f}")
    print("  -> same exponent; at fixed per-center basis the win is a CONSTANT")
    print("     factor (M/S)^4, not a changed exponent.")

    # -------- basis-richness axis: max_n sweep on a single center --------
    print("\n" + "=" * 92)
    print("SCALING vs basis RICHNESS (max_n sweep, SINGLE center)")
    print("  This is the axis where the m-collapse win GROWS.")
    print("=" * 92)
    from geovac.molecular_spec import MolecularSpec, OrbitalBlock
    Ns, Msi, Ssi, gasi, gvsi = [], [], [], [], []
    print(f"{'max_n':>6} {'l_max':>6} {'M(orb)':>7} {'S(shell)':>9} "
          f"{'M/S':>6} {'GAUSS':>14} {'GEOVAC':>12} {'ratio':>9} {'(M/S)^4':>9}")
    for N in range(2, 7):
        spec = MolecularSpec(
            name=f"1c_n{N}",
            blocks=[OrbitalBlock(label="c", block_type="atomic",
                                 Z_center=1.0, n_electrons=2, max_n=N)],
            nuclear_repulsion_constant=0.0)
        subs = subblocks_for_spec(spec)
        orbitals, shells, shell_l, shell_sub = basis_summary(subs)
        M = len(orbitals)
        S = len(shells)
        gl = gaussian_load(M)["total"]
        gvv = geovac_load(shells, shell_l, shell_sub)["total_primary"]
        Ns.append(N); Msi.append(M); Ssi.append(S); gasi.append(gl); gvsi.append(gvv)
        print(f"{N:>6} {N-1:>6} {M:>7} {S:>9} {M/S:>6.2f} "
              f"{gl:>14,} {gvv:>12,} {gl/gvv:>9.1f} {(M/S)**4:>9.1f}")

    bGr, _ = loglog_fit(Ns, gasi)
    bVr, _ = loglog_fit(Ns, gvsi)
    ratio_fit_b, ratio_fit_a = loglog_fit(Ns, [g / v for g, v in zip(gasi, gvsi)])
    print(f"\n  Gaussian LOAD ~ max_n^{bGr:.2f}")
    print(f"  GeoVac   LOAD ~ max_n^{bVr:.2f}")
    print(f"  LOAD RATIO (Gauss/GeoVac) ~ {ratio_fit_a:.3g} * max_n^{ratio_fit_b:.2f}"
          f"   <-- the win GROWS polynomially with angular richness")

    print("\n" + "=" * 92)
    print("VERDICT")
    print("=" * 92)
    print("""GO (with a sharp caveat on WHICH axis).

 * The m-degeneracy collapse is real and exact: GeoVac loads radial seeds indexed
   by shells (n,l), not orbitals (n,l,m). GeoVac LOAD = (S/M)^4 x Gaussian LOAD.
 * vs CENTER count at fixed per-center basis (max_n=2): SAME scaling exponent
   (~nsub^4 for both); GeoVac wins by a CONSTANT factor ~ (M/S)^4 ~ 7-8x.
 * vs BASIS RICHNESS (max_n / l_max): GeoVac LOAD is asymptotically SUBLINEAR in
   the Gaussian tensor size — the ratio grows ~ (2*max_n/3)^4. THIS is the
   'markedly slower' axis the gate asks for.
 * GENERATE is a fixed O(1)-code rule (a handful of distinct 3j values, l_max only)
   — it does not scale with M or center count at all.
 HONESTY CAVEAT: the Gaussian baseline is dense M^4 (fair for these small systems).
   Large systems get Gaussian integral SCREENING (distance sparsity, 'two kinds of
   sparsity' memo) which erodes the dense-M^4 baseline; the comparison is fairest
   in the small-molecule regime, which is GeoVac's regime anyway.""")


if __name__ == "__main__":
    main()
