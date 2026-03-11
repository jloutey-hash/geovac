"""
MOSturmianFCI — Full CI in molecular Sturmian basis (v0.9.34)
=============================================================

Full Configuration Interaction solver using molecular Sturmian orbitals
computed from prolate spheroidal coordinate separation. This is a
fundamentally different architecture from the atom-centered LCAO approach
in lattice_index.py and direct_ci.py.

v0.9.34 changes:
  - Dual-p0 mode: combines Li-scale (p0_Li, nmax_Li) and H-scale
    (p0_H, nmax_H) MO Sturmian sets into a single FCI. Canonical
    orthogonalization of the combined basis eliminates near-linear
    dependence. Two-parameter self-consistency loop updates both
    p0_Li and p0_H from projected 1-RDM subblocks.

v0.9.33 changes:
  - Exact exchange integrals K_kl via Poisson solve.
  - Energy decomposition method decompose_energy().

Self-consistent p0 loop (single-p0 mode):
  1. compute_molecular_sturmian_betas(p0)
  2. compute_h1_matrix(..., orth_method='canonical')
  3. compute_exact_j_integrals(...) + compute_exact_k_integrals(...)
  4. build SD basis, assemble FCI Hamiltonian
  5. diagonalize → E_mol
  6. p0_new = sqrt(-2*E_mol)
  7. iterate to convergence

Author: GeoVac Development Team
Date: March 2026
"""

from __future__ import annotations

import time
from itertools import combinations
from typing import List, Optional, Tuple

import numpy as np
from scipy.sparse import coo_matrix, csr_matrix
from scipy.sparse.linalg import eigsh

from .molecular_sturmian import (
    compute_molecular_sturmian_betas,
    compute_h1_matrix,
    compute_eri_matrix,
    compute_exact_j_integrals,
    compute_exact_k_integrals,
    compute_cross_set_integrals,
    compute_combined_jk_integrals,
)


class MOSturmianFCI:
    """Full CI in molecular Sturmian basis.

    Parameters
    ----------
    Z_A, Z_B : float
        Nuclear charges (Z_A >= Z_B by convention).
    R : float
        Internuclear distance in bohr.
    nmax : int
        Maximum principal quantum number for basis (Li-scale in dual mode).
    n_electrons : int
        Number of electrons.
    n_grid : int
        Gauss-Legendre quadrature points per dimension.
    xi_max : float
        Upper radial integration bound.
    dual_p0 : bool
        If True, use dual-p0 basis combining Li-scale and H-scale MO sets.
    nmax_H : int
        Maximum principal quantum number for H-scale set (only used when
        dual_p0=True). Default 2 to avoid near-linear dependence.
    """

    def __init__(
        self,
        Z_A: float,
        Z_B: float,
        R: float,
        nmax: int = 3,
        n_electrons: int = 4,
        n_grid: int = 40,
        xi_max: float = 12.0,
        dual_p0: bool = False,
        nmax_H: int = 2,
    ) -> None:
        self.Z_A = Z_A
        self.Z_B = Z_B
        self.R = R
        self.nmax = nmax
        self.n_electrons = n_electrons
        self.n_grid = n_grid
        self.xi_max = xi_max
        self.dual_p0 = dual_p0
        self.nmax_H = nmax_H

        # Will be set during solve()
        self.n_mo: int = 0
        self.n_spinorb: int = 0
        self.n_sd: int = 0
        self.mo_results: List[Tuple[int, int, int, int, float]] = []
        self.sd_basis: List[Tuple[int, ...]] = []
        self._sd_index: dict = {}
        self.H1: Optional[np.ndarray] = None
        self.ERI: Optional[np.ndarray] = None
        self.energy: Optional[float] = None
        self.p0_converged: Optional[float] = None
        self.civec: Optional[np.ndarray] = None

        # Dual-p0 specific state
        self.mo_results_Li: List[Tuple[int, int, int, int, float]] = []
        self.mo_results_H: List[Tuple[int, int, int, int, float]] = []
        self.n_mo_Li: int = 0
        self.n_mo_H: int = 0
        self.p0_Li_converged: Optional[float] = None
        self.p0_H_converged: Optional[float] = None

    def _build_basis(self, p0: float) -> bool:
        """Compute MO betas and build SD basis for given p0.

        Returns True if all MOs found, False otherwise.
        """
        self.mo_results = compute_molecular_sturmian_betas(
            Z_A=self.Z_A, Z_B=self.Z_B, R=self.R,
            p0=p0, nmax=self.nmax,
        )

        # Filter valid MOs
        valid = [mo for mo in self.mo_results if np.isfinite(mo[4])]
        self.n_mo = len(valid)
        self.n_spinorb = 2 * self.n_mo
        self._valid_mo_indices = [i for i, mo in enumerate(self.mo_results)
                                  if np.isfinite(mo[4])]

        if self.n_mo < self.n_electrons // 2:
            print(f"  [MO-FCI] Only {self.n_mo} MOs found, need >= {self.n_electrons // 2}")
            return False

        # Build Slater determinant basis
        # Spin-orbital indexing: 2*spatial + spin (0=alpha, 1=beta)
        self.sd_basis = list(combinations(range(self.n_spinorb), self.n_electrons))
        self.n_sd = len(self.sd_basis)
        self._sd_index = {sd: i for i, sd in enumerate(self.sd_basis)}

        return True

    def _compute_integrals(self, p0: float, verbose: bool = False) -> None:
        """Compute orthogonalized H1 and ERI matrices for current p0.

        Uses canonical orthogonalization (drops near-LD eigenvectors),
        exact direct Coulomb integrals J via Poisson solve, and
        exact exchange integrals K via Poisson solve.
        """
        self.H1, X, n_dropped, S, H1_raw = compute_h1_matrix(
            self.mo_results, self.Z_A, self.Z_B, self.R, p0,
            n_grid=self.n_grid, xi_max=self.xi_max,
            orth_method='canonical', orth_threshold=1e-4,
        )
        self._X = X  # store for ERI transform
        self._S = S  # store for diagnostics
        self._H1_raw = H1_raw
        self._n_dropped = n_dropped

        if n_dropped > 0 and verbose:
            print(f"  [MO-FCI] Canonical orth dropped {n_dropped} near-LD vectors")

        # Update dimensions after canonical orthogonalization
        n_orth = self.H1.shape[0]
        self.n_mo = n_orth
        self.n_spinorb = 2 * n_orth
        # Rebuild SD basis with reduced MO count
        self.sd_basis = list(combinations(range(self.n_spinorb), self.n_electrons))
        self.n_sd = len(self.sd_basis)
        self._sd_index = {sd: i for i, sd in enumerate(self.sd_basis)}

        # --- Two-electron integrals ---
        # Exact J via Poisson solve
        J_exact = compute_exact_j_integrals(
            self.mo_results, self.Z_A, self.Z_B, self.R, p0,
            n_grid=self.n_grid, xi_max=self.xi_max,
        )
        self._J_raw = J_exact  # store for diagnostics (raw MO basis)

        # Exact K via Poisson solve (replaces Ohno-Klopman)
        K_exact = compute_exact_k_integrals(
            self.mo_results, self.Z_A, self.Z_B, self.R, p0,
            n_grid=self.n_grid, xi_max=self.xi_max,
        )
        self._K_raw = K_exact  # store for diagnostics (raw MO basis)

        # Build ERI tensor with exact J and K
        valid_mos = [(i, mo) for i, mo in enumerate(self.mo_results)
                     if np.isfinite(mo[4])]
        n_mo_raw = len(valid_mos)
        ERI_raw = np.zeros((n_mo_raw, n_mo_raw, n_mo_raw, n_mo_raw))

        for k in range(n_mo_raw):
            for l in range(n_mo_raw):
                # Direct Coulomb: <kl|kl> = J_kl
                ERI_raw[k, l, k, l] = J_exact[k, l]
                # Exchange: <kl|lk> = K_kl
                ERI_raw[k, l, l, k] = K_exact[k, l]

        # Transform ERI to orthogonal basis using X
        tmp = np.einsum('ai,abcd->ibcd', X, ERI_raw)
        tmp = np.einsum('bj,ibcd->ijcd', X, tmp)
        tmp = np.einsum('ck,ijcd->ijkd', X, tmp)
        self.ERI = np.einsum('dl,ijkd->ijkl', X, tmp)

    def diagnose_orthogonalization(self, p0: float, R: float) -> dict:
        """Diagnose overlap matrix and H1 before/after orthogonalization.

        Parameters
        ----------
        p0 : float
            Momentum parameter.
        R : float
            Internuclear distance (uses self.R if different, but stored).

        Returns
        -------
        diag : dict with keys:
            'S': overlap matrix
            'S_evals': eigenvalues of S
            'H1_raw_diag': H1 diagonal before orthogonalization
            'H1_lowdin_diag': H1 diagonal after Löwdin
            'H1_canonical_diag': H1 diagonal after canonical
            'n_dropped': number dropped by canonical
            'top_offdiag': list of (i, j, S_ij) for largest off-diagonal elements
        """
        mo_results = compute_molecular_sturmian_betas(
            Z_A=self.Z_A, Z_B=self.Z_B, R=R,
            p0=p0, nmax=self.nmax,
        )

        # Löwdin
        H1_low, X_low, _, S, H1_raw = compute_h1_matrix(
            mo_results, self.Z_A, self.Z_B, R, p0,
            n_grid=self.n_grid, xi_max=self.xi_max,
            orth_method='lowdin',
        )

        # Canonical
        H1_can, X_can, n_dropped, _, _ = compute_h1_matrix(
            mo_results, self.Z_A, self.Z_B, R, p0,
            n_grid=self.n_grid, xi_max=self.xi_max,
            orth_method='canonical', orth_threshold=1e-4,
        )

        evals_S = np.sort(np.linalg.eigvalsh(S))

        # Top off-diagonal elements
        n = S.shape[0]
        offdiag = []
        for i in range(n):
            for j in range(i + 1, n):
                offdiag.append((i, j, S[i, j]))
        offdiag.sort(key=lambda x: abs(x[2]), reverse=True)

        return {
            'S': S,
            'S_evals': evals_S,
            'H1_raw_diag': np.diag(H1_raw),
            'H1_lowdin_diag': np.diag(H1_low),
            'H1_canonical_diag': np.diag(H1_can),
            'n_dropped': n_dropped,
            'top_offdiag': offdiag[:5],
            'mo_results': mo_results,
        }

    def _assemble_hamiltonian(self) -> csr_matrix:
        """Build FCI Hamiltonian using Slater-Condon rules.

        For 4845 SDs this is fast enough with a direct loop.
        Uses dense H1 and ERI arrays.
        """
        n_sd = self.n_sd
        n_el = self.n_electrons
        H1 = self.H1
        eri = self.ERI
        sd_basis = self.sd_basis
        sd_index = self._sd_index
        V_NN = self.Z_A * self.Z_B / self.R

        rows: List[int] = []
        cols: List[int] = []
        vals: List[float] = []

        # --- Diagonal elements ---
        for I, sd_I in enumerate(sd_basis):
            h_diag = 0.0
            # One-electron part
            for p in sd_I:
                sp_p = p >> 1  # spatial index
                h_diag += H1[sp_p, sp_p]

            # Two-electron part: Coulomb - Exchange
            for i in range(n_el):
                pi = sd_I[i]
                sp_i = pi >> 1
                sig_i = pi & 1
                for j in range(i + 1, n_el):
                    pj = sd_I[j]
                    sp_j = pj >> 1
                    sig_j = pj & 1
                    # Coulomb: <ij|ij>
                    h_diag += eri[sp_i, sp_j, sp_i, sp_j]
                    # Exchange: -delta(sigma_i, sigma_j) <ij|ji>
                    if sig_i == sig_j:
                        h_diag -= eri[sp_i, sp_j, sp_j, sp_i]

            # Nuclear repulsion (constant shift)
            h_diag += V_NN

            rows.append(I)
            cols.append(I)
            vals.append(h_diag)

        # --- Off-diagonal: singles and doubles ---
        for I in range(n_sd):
            sd_I = sd_basis[I]
            occ_I = frozenset(sd_I)

            # Singles: p -> r
            for kp in range(n_el):
                p = sd_I[kp]
                sp_p = p >> 1
                sig_p = p & 1

                for r in range(self.n_spinorb):
                    if r in occ_I:
                        continue
                    if (r & 1) != sig_p:
                        continue  # spin conservation

                    sp_r = r >> 1

                    # Slater-Condon single excitation matrix element
                    me = H1[sp_p, sp_r]
                    for q in sd_I:
                        if q == p:
                            continue
                        sp_q = q >> 1
                        me += eri[sp_p, sp_q, sp_r, sp_q]
                        if (q & 1) == sig_p:
                            me -= eri[sp_p, sp_q, sp_q, sp_r]

                    if abs(me) < 1e-12:
                        continue

                    # Build target SD and get phase
                    new_sd = list(sd_I)
                    new_sd[kp] = r
                    new_sd_sorted = tuple(sorted(new_sd))
                    J = sd_index.get(new_sd_sorted)
                    if J is None or J <= I:
                        continue  # upper triangle only

                    phase = self._compute_phase(sd_I, kp, r)
                    rows.append(I)
                    cols.append(J)
                    vals.append(phase * me)
                    rows.append(J)
                    cols.append(I)
                    vals.append(phase * me)

            # Doubles: (p,q) -> (r,s)
            for kp in range(n_el):
                p = sd_I[kp]
                sp_p = p >> 1
                sig_p = p & 1

                for kq in range(kp + 1, n_el):
                    q = sd_I[kq]
                    sp_q = q >> 1
                    sig_q = q & 1

                    for r in range(self.n_spinorb):
                        if r in occ_I or (r & 1) != sig_p:
                            continue
                        sp_r = r >> 1

                        for s in range(r + 1, self.n_spinorb):
                            if s in occ_I or s == r:
                                continue
                            if (s & 1) != sig_q:
                                continue
                            sp_s = s >> 1

                            # <pq|rs> - <pq|sr>
                            me = eri[sp_p, sp_q, sp_r, sp_s]
                            if sig_p == sig_q:
                                me -= eri[sp_p, sp_q, sp_s, sp_r]

                            if abs(me) < 1e-12:
                                continue

                            # Build target SD
                            new_sd = list(sd_I)
                            new_sd[kp] = r
                            new_sd[kq] = s
                            new_sd_sorted = tuple(sorted(new_sd))
                            J = sd_index.get(new_sd_sorted)
                            if J is None or J <= I:
                                continue

                            phase = self._compute_double_phase(sd_I, kp, kq, r, s)
                            rows.append(I)
                            cols.append(J)
                            vals.append(phase * me)
                            rows.append(J)
                            cols.append(I)
                            vals.append(phase * me)

        H = coo_matrix((vals, (rows, cols)), shape=(n_sd, n_sd))
        return H.tocsr()

    @staticmethod
    def _compute_phase(sd: Tuple[int, ...], kp: int, r: int) -> int:
        """Phase from single excitation: sd[kp] -> r.

        Count how many occupied orbitals lie between sd[kp] and r.
        """
        p = sd[kp]
        # Count inversions: number of elements in sd between p and r
        count = 0
        for i, occ in enumerate(sd):
            if i == kp:
                continue
            if (p < occ < r) or (r < occ < p):
                count += 1
        return (-1)**count

    @staticmethod
    def _compute_double_phase(
        sd: Tuple[int, ...], kp: int, kq: int, r: int, s: int,
    ) -> int:
        """Phase from double excitation: sd[kp]->r, sd[kq]->s.

        Apply two single excitations and count total inversions.
        """
        p = sd[kp]
        q = sd[kq]

        # First excitation: p -> r
        count1 = 0
        for i, occ in enumerate(sd):
            if i == kp:
                continue
            if (p < occ < r) or (r < occ < p):
                count1 += 1

        # After first excitation, build intermediate SD
        inter = list(sd)
        inter[kp] = r
        inter_sorted = sorted(inter)

        # Second excitation: q -> s in intermediate
        # Find position of q in intermediate
        count2 = 0
        for occ in inter_sorted:
            if occ == q:
                continue
            if (q < occ < s) or (s < occ < q):
                count2 += 1

        return (-1)**(count1 + count2)

    def decompose_energy(self) -> dict:
        """Decompose FCI energy into H1, J, K, V_NN contributions.

        Uses diagonal-weighted CI decomposition: for each Slater determinant I,
        weight by |c_I|^2 and sum the one-electron, Coulomb, exchange, and
        nuclear repulsion contributions.

        Returns
        -------
        decomp : dict with keys 'E_H1', 'E_J', 'E_K', 'V_NN', 'E_sum',
            'E_fci', 'residual' (E_fci - E_sum, from off-diagonal CI mixing).
        """
        if self.civec is None or self.energy is None:
            raise RuntimeError("Call solve() before decompose_energy()")

        c = self.civec
        H1 = self.H1
        eri = self.ERI
        V_NN = self.Z_A * self.Z_B / self.R
        n_el = self.n_electrons

        E_H1 = 0.0
        E_J = 0.0
        E_K = 0.0

        for I, sd_I in enumerate(self.sd_basis):
            w = c[I] ** 2

            # One-electron contribution
            for p in sd_I:
                sp_p = p >> 1
                E_H1 += w * H1[sp_p, sp_p]

            # Two-electron contribution
            for i in range(n_el):
                pi = sd_I[i]
                sp_i = pi >> 1
                sig_i = pi & 1
                for j in range(i + 1, n_el):
                    pj = sd_I[j]
                    sp_j = pj >> 1
                    sig_j = pj & 1
                    # Coulomb
                    E_J += w * eri[sp_i, sp_j, sp_i, sp_j]
                    # Exchange (same spin only)
                    if sig_i == sig_j:
                        E_K -= w * eri[sp_i, sp_j, sp_j, sp_i]

        E_sum = E_H1 + E_J + E_K + V_NN
        return {
            'E_H1': E_H1,
            'E_J': E_J,
            'E_K': E_K,
            'V_NN': V_NN,
            'E_sum': E_sum,
            'E_fci': self.energy,
            'residual': self.energy - E_sum,
        }

    # ------------------------------------------------------------------
    # Dual-p0 methods (v0.9.34)
    # ------------------------------------------------------------------

    def _build_basis_dual(
        self, p0_Li: float, p0_H: float,
    ) -> bool:
        """Build combined Li-scale + H-scale MO basis.

        Returns True if enough MOs found, False otherwise.
        """
        self.mo_results_Li = compute_molecular_sturmian_betas(
            Z_A=self.Z_A, Z_B=self.Z_B, R=self.R,
            p0=p0_Li, nmax=self.nmax,
        )
        self.mo_results_H = compute_molecular_sturmian_betas(
            Z_A=self.Z_A, Z_B=self.Z_B, R=self.R,
            p0=p0_H, nmax=self.nmax_H,
        )

        valid_Li = [mo for mo in self.mo_results_Li if np.isfinite(mo[4])]
        valid_H = [mo for mo in self.mo_results_H if np.isfinite(mo[4])]
        self.n_mo_Li = len(valid_Li)
        self.n_mo_H = len(valid_H)
        n_mo_raw = self.n_mo_Li + self.n_mo_H

        if n_mo_raw < self.n_electrons // 2:
            print(f"  [DUAL] Only {n_mo_raw} MOs found, need >= "
                  f"{self.n_electrons // 2}")
            return False

        self.n_mo = n_mo_raw
        self.n_spinorb = 2 * n_mo_raw
        self.sd_basis = list(
            combinations(range(self.n_spinorb), self.n_electrons)
        )
        self.n_sd = len(self.sd_basis)
        self._sd_index = {sd: i for i, sd in enumerate(self.sd_basis)}

        return True

    def _compute_integrals_dual(
        self, p0_Li: float, p0_H: float, verbose: bool = False,
    ) -> None:
        """Compute combined H1, overlap, and ERI for dual-p0 basis.

        Block structure:
          S = [S_LL  S_LH]    H1 = [H1_LL  H1_LH]
              [S_LH' S_HH]         [H1_LH' H1_HH]
        """
        # --- Within-set integrals ---
        H1_LL, X_LL, nd_LL, S_LL, H1_LL_raw = compute_h1_matrix(
            self.mo_results_Li, self.Z_A, self.Z_B, self.R, p0_Li,
            n_grid=self.n_grid, xi_max=self.xi_max,
            orth_method='none_raw',  # we do orth on combined basis
        )
        H1_HH, X_HH, nd_HH, S_HH, H1_HH_raw = compute_h1_matrix(
            self.mo_results_H, self.Z_A, self.Z_B, self.R, p0_H,
            n_grid=self.n_grid, xi_max=self.xi_max,
            orth_method='none_raw',
        )

        # --- Cross-set integrals ---
        S_LH, H1_LH = compute_cross_set_integrals(
            self.mo_results_Li, self.mo_results_H,
            self.Z_A, self.Z_B, self.R, p0_Li, p0_H,
            n_grid=self.n_grid, xi_max=self.xi_max,
        )

        # --- Assemble full combined matrices ---
        nL = S_LL.shape[0]
        nH = S_HH.shape[0]
        n_total = nL + nH

        S_combined = np.zeros((n_total, n_total))
        S_combined[:nL, :nL] = S_LL
        S_combined[nL:, nL:] = S_HH
        S_combined[:nL, nL:] = S_LH
        S_combined[nL:, :nL] = S_LH.T

        H1_combined = np.zeros((n_total, n_total))
        H1_combined[:nL, :nL] = H1_LL_raw
        H1_combined[nL:, nL:] = H1_HH_raw
        H1_combined[:nL, nL:] = H1_LH
        H1_combined[nL:, :nL] = H1_LH.T

        self._S_combined = S_combined
        self._H1_combined_raw = H1_combined
        self._n_mo_Li_raw = nL
        self._n_mo_H_raw = nH

        # --- Canonical orthogonalization ---
        s_evals, s_evecs = np.linalg.eigh(S_combined)
        orth_threshold = 1e-3
        keep = s_evals > orth_threshold
        n_dropped = int(np.sum(~keep))
        self._n_dropped = n_dropped
        self._S_evals = s_evals

        if verbose:
            print(f"  [DUAL] Combined basis: {nL} Li + {nH} H = "
                  f"{n_total} raw MOs")
            print(f"  [DUAL] S eigenvalues: min={s_evals[0]:.4f}, "
                  f"max={s_evals[-1]:.4f}")
            if n_dropped > 0:
                print(f"  [DUAL] Dropped {n_dropped} near-LD vectors "
                      f"(threshold={orth_threshold})")
            if n_dropped > 4:
                print(f"  [DUAL] WARNING: n_dropped={n_dropped} > 4, "
                      f"reduce nmax_H")

        U_keep = s_evecs[:, keep]
        s_keep = s_evals[keep]
        X = U_keep @ np.diag(s_keep**(-0.5))

        H1_orth = X.T @ H1_combined @ X
        self.H1 = H1_orth
        self._X = X

        # Update basis dimensions after canonical orth
        n_orth = H1_orth.shape[0]
        self.n_mo = n_orth
        self.n_spinorb = 2 * n_orth
        self.sd_basis = list(
            combinations(range(self.n_spinorb), self.n_electrons)
        )
        self.n_sd = len(self.sd_basis)
        self._sd_index = {sd: i for i, sd in enumerate(self.sd_basis)}

        if verbose:
            print(f"  [DUAL] After orth: {n_orth} MOs, {self.n_spinorb} "
                  f"spinorbs, {self.n_sd} SDs")

        # --- Two-electron integrals (combined J, K) ---
        J_combined, K_combined = compute_combined_jk_integrals(
            self.mo_results_Li, self.mo_results_H,
            self.Z_A, self.Z_B, self.R, p0_Li, p0_H,
            n_grid=self.n_grid, xi_max=self.xi_max,
        )
        self._J_raw = J_combined
        self._K_raw = K_combined

        # Build ERI tensor
        ERI_raw = np.zeros((n_total, n_total, n_total, n_total))
        for k in range(n_total):
            for l in range(n_total):
                ERI_raw[k, l, k, l] = J_combined[k, l]
                ERI_raw[k, l, l, k] = K_combined[k, l]

        # Transform to orthogonal basis
        tmp = np.einsum('ai,abcd->ibcd', X, ERI_raw)
        tmp = np.einsum('bj,ibcd->ijcd', X, tmp)
        tmp = np.einsum('ck,ijcd->ijkd', X, tmp)
        self.ERI = np.einsum('dl,ijkd->ijkl', X, tmp)

    def solve_dual(
        self,
        p0_Li_init: Optional[float] = None,
        p0_H_init: float = 1.0,
        damping: float = 0.3,
        max_iter: int = 20,
        tol: float = 1e-4,
        verbose: bool = True,
    ) -> Tuple[float, float, float]:
        """Two-parameter self-consistent p0 loop for dual-p0 FCI.

        Parameters
        ----------
        p0_Li_init : float or None
            Initial Li-scale momentum. Default: sqrt(2 * 7.392).
        p0_H_init : float
            Initial H-scale momentum. Default: 1.0 (exact hydrogen).
        damping : float
            Damping factor for p0 updates.
        max_iter : int
            Maximum iterations.
        tol : float
            Convergence on |dp0_Li| + |dp0_H|.
        verbose : bool
            Print iteration details.

        Returns
        -------
        E_mol : float
            Converged molecular energy.
        p0_Li_star : float
            Converged Li-scale momentum.
        p0_H_star : float
            Converged H-scale momentum.
        """
        if p0_Li_init is None:
            p0_Li_init = np.sqrt(2.0 * 7.392)  # Li nmax=3 atomic

        p0_Li = p0_Li_init
        p0_H = p0_H_init

        if verbose:
            print(f"[DUAL-FCI] Z_A={self.Z_A}, Z_B={self.Z_B}, "
                  f"R={self.R:.3f}")
            print(f"[DUAL-FCI] nmax_Li={self.nmax}, nmax_H={self.nmax_H}, "
                  f"n_el={self.n_electrons}")
            print(f"[DUAL-FCI] Starting p0_Li={p0_Li:.4f}, "
                  f"p0_H={p0_H:.4f}")

        t_total = time.perf_counter()

        for iteration in range(max_iter):
            t0 = time.perf_counter()

            # Step 1: Build combined basis
            if not self._build_basis_dual(p0_Li, p0_H):
                print(f"  [DUAL] iter {iteration}: basis build failed")
                break

            # Step 2: Compute all integrals
            self._compute_integrals_dual(p0_Li, p0_H, verbose=(
                verbose and iteration == 0
            ))

            # Step 3: Assemble and diagonalize
            H_fci = self._assemble_hamiltonian()

            if self.n_sd <= 200:
                H_dense = H_fci.toarray()
                evals, evecs = np.linalg.eigh(H_dense)
                E_mol = evals[0]
                self.civec = evecs[:, 0]
            else:
                evals, evecs = eigsh(H_fci, k=1, which='SA')
                E_mol = evals[0]
                self.civec = evecs[:, 0]

            self.energy = E_mol

            # Step 4: Project 1-RDM onto Li and H subblocks
            nL = self._n_mo_Li_raw
            nH = self._n_mo_H_raw
            X = self._X
            S_comb = self._S_combined

            # Build 1-RDM in orthogonal basis
            gamma_orth = np.zeros((self.n_mo, self.n_mo))
            for I, sd_I in enumerate(self.sd_basis):
                for p in sd_I:
                    sp_p = p >> 1
                    gamma_orth[sp_p, sp_p] += self.civec[I]**2

            # Transform to raw basis: gamma_raw = X @ gamma_orth @ X.T
            gamma_raw = X @ gamma_orth @ X.T

            # Li subblock: trace with H1 and S
            H1_raw = self._H1_combined_raw
            gamma_Li = gamma_raw[:nL, :nL]
            S_Li = S_comb[:nL, :nL]
            E_Li_proj = np.trace(gamma_Li @ H1_raw[:nL, :nL])
            N_Li_el = np.trace(gamma_Li @ S_Li)

            # H subblock
            gamma_H = gamma_raw[nL:, nL:]
            S_H = S_comb[nL:, nL:]
            E_H_proj = np.trace(gamma_H @ H1_raw[nL:, nL:])
            N_H_el = np.trace(gamma_H @ S_H)

            # Step 5: Update p0_Li and p0_H
            if E_mol >= 0:
                if verbose:
                    print(f"  iter {iteration}: E={E_mol:.6f} Ha "
                          f"(POSITIVE)")
                self.p0_Li_converged = p0_Li
                self.p0_H_converged = p0_H
                break

            # p0 from projected energy
            if N_Li_el > 0.1 and E_Li_proj < 0:
                p0_Li_new = np.sqrt(
                    -2.0 * E_Li_proj / max(N_Li_el, 0.5)
                )
            else:
                p0_Li_new = p0_Li  # no update if projection fails

            if N_H_el > 0.1 and E_H_proj < 0:
                p0_H_new = np.sqrt(
                    -2.0 * E_H_proj / max(N_H_el, 0.5)
                )
            else:
                p0_H_new = p0_H

            dp0_Li = abs(p0_Li_new - p0_Li)
            dp0_H = abs(p0_H_new - p0_H)

            p0_Li_old = p0_Li
            p0_H_old = p0_H
            p0_Li = (1.0 - damping) * p0_Li + damping * p0_Li_new
            p0_H = (1.0 - damping) * p0_H + damping * p0_H_new

            if verbose:
                print(f"  iter {iteration}: E={E_mol:.6f} Ha, "
                      f"p0_Li={p0_Li_old:.4f}->{p0_Li:.4f} "
                      f"(dp={dp0_Li:.4f}), "
                      f"p0_H={p0_H_old:.4f}->{p0_H:.4f} "
                      f"(dp={dp0_H:.4f}), "
                      f"N_Li={N_Li_el:.2f}, N_H={N_H_el:.2f}, "
                      f"n_sd={self.n_sd}, "
                      f"t={time.perf_counter()-t0:.1f}s")

            self.p0_Li_converged = p0_Li
            self.p0_H_converged = p0_H
            self.p0_converged = p0_Li  # backward compat

            if dp0_Li + dp0_H < tol:
                if verbose:
                    print(f"  [DUAL] Converged at iter {iteration}, "
                          f"dp_sum={dp0_Li+dp0_H:.2e}")
                break

        t_elapsed = time.perf_counter() - t_total
        if verbose:
            print(f"[DUAL-FCI] Done: E={self.energy:.6f} Ha, "
                  f"p0_Li*={self.p0_Li_converged:.4f}, "
                  f"p0_H*={self.p0_H_converged:.4f}, "
                  f"total={t_elapsed:.1f}s")

        return self.energy, self.p0_Li_converged, self.p0_H_converged

    def solve(
        self,
        p0_init: Optional[float] = None,
        damping: float = 0.5,
        max_iter: int = 20,
        tol: float = 1e-5,
        verbose: bool = True,
    ) -> Tuple[float, float]:
        """Self-consistent p0 loop for MO Sturmian FCI.

        Parameters
        ----------
        p0_init : float or None
            Initial momentum parameter. Default: sqrt(-2*E_atoms).
        damping : float
            Damping factor for p0 update (0=no update, 1=full update).
        max_iter : int
            Maximum SCF iterations.
        tol : float
            Convergence threshold on p0.
        verbose : bool
            Print iteration details.

        Returns
        -------
        E_mol : float
            Converged molecular energy in Hartree.
        p0_star : float
            Converged momentum parameter.
        """
        if p0_init is None:
            # Start from sum of atomic energies
            E_atoms = -(self.Z_A**2) / 2.0 - (self.Z_B**2) / 2.0
            p0 = np.sqrt(-2.0 * E_atoms)
        else:
            p0 = p0_init

        if verbose:
            print(f"[MO-FCI] Z_A={self.Z_A}, Z_B={self.Z_B}, R={self.R:.3f}, "
                  f"nmax={self.nmax}, n_el={self.n_electrons}")
            print(f"[MO-FCI] Starting p0={p0:.4f}")

        t_total = time.perf_counter()

        for iteration in range(max_iter):
            t0 = time.perf_counter()

            # Step 1: Compute MO betas (sets n_mo to raw count)
            if not self._build_basis(p0):
                print(f"  [MO-FCI] iter {iteration}: basis build failed at p0={p0:.4f}")
                break

            n_mo_raw = self.n_mo  # before canonical orth

            # Step 2: Compute integrals (may reduce n_mo via canonical orth)
            self._compute_integrals(p0, verbose=verbose)

            t_integrals = time.perf_counter() - t0

            # Step 3: Assemble and diagonalize
            t1 = time.perf_counter()
            H_fci = self._assemble_hamiltonian()
            t_assemble = time.perf_counter() - t1

            t2 = time.perf_counter()
            if self.n_sd <= 100:
                # Small enough for dense diagonalization
                H_dense = H_fci.toarray()
                evals, evecs = np.linalg.eigh(H_dense)
                E_mol = evals[0]
                self.civec = evecs[:, 0]
            else:
                evals, evecs = eigsh(H_fci, k=1, which='SA')
                E_mol = evals[0]
                self.civec = evecs[:, 0]
            t_diag = time.perf_counter() - t2

            # Step 4: Update p0
            if E_mol >= 0:
                if verbose:
                    print(f"  iter {iteration}: E={E_mol:.6f} Ha (POSITIVE), "
                          f"p0={p0:.4f}, n_sd={self.n_sd}")
                # Don't update p0 if energy is positive
                self.energy = E_mol
                self.p0_converged = p0
                break

            p0_new = np.sqrt(-2.0 * E_mol)
            dp0 = abs(p0_new - p0)
            p0_old = p0
            p0 = (1.0 - damping) * p0 + damping * p0_new

            if verbose:
                print(f"  iter {iteration}: E={E_mol:.6f} Ha, p0={p0_old:.4f}->{p0:.4f} "
                      f"(dp0={dp0:.6f}), n_mo={self.n_mo}, n_sd={self.n_sd}, "
                      f"t={time.perf_counter()-t0:.2f}s")

            self.energy = E_mol
            self.p0_converged = p0

            if dp0 < tol:
                if verbose:
                    print(f"  [MO-FCI] Converged at iter {iteration}, "
                          f"dp0={dp0:.2e} < tol={tol:.0e}")
                break

        t_elapsed = time.perf_counter() - t_total
        if verbose:
            print(f"[MO-FCI] Done: E={self.energy:.6f} Ha, p0*={self.p0_converged:.4f}, "
                  f"total={t_elapsed:.2f}s")

        return self.energy, self.p0_converged
