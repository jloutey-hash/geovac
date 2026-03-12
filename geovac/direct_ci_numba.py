"""
Numba-accelerated FCI Hamiltonian assembly kernels
===================================================

Drop-in replacement for the pure-Python loops in DirectCISolver.assemble_hamiltonian().

Strategy:
  - @njit gives ~100x over CPython for tight numerical loops (no parallel needed).
  - SD→index via combinatorial number system (replaces dict lookup).
  - All data as dense NumPy arrays (already the case in DirectCISolver).
  - Single monolithic kernel: diagonal + singles + doubles in one pass.

Author: GeoVac Development Team
Date: March 2026
"""

import time

import numpy as np

try:
    from numba import njit
    NUMBA_AVAILABLE = True
except ImportError:
    NUMBA_AVAILABLE = False


# ======================================================================
# Data preparation (pure Python — runs once)
# ======================================================================

def sd_basis_to_array(sd_basis: list) -> np.ndarray:
    """Convert list of sorted tuples to contiguous 2D int32 array."""
    n_sd = len(sd_basis)
    n_el = len(sd_basis[0])
    arr = np.empty((n_sd, n_el), dtype=np.int32)
    for i, sd in enumerate(sd_basis):
        for j in range(n_el):
            arr[i, j] = sd[j]
    return arr


def build_colex_to_lex(sd_basis: list, n_sp: int) -> np.ndarray:
    """
    Build mapping from combinatorial (colex) index to sequential (lex) index.

    itertools.combinations gives lex order; the combinatorial number system
    C(occ[0],1) + C(occ[1],2) + ... gives colex order. This array maps between
    them: colex_to_lex[colex_idx] = lex_idx.

    Returns int64 array of shape (n_sd,).
    """
    n_sd = len(sd_basis)
    n_el = len(sd_basis[0])
    binom = _precompute_binomial(n_sp, n_el)
    sd_array = sd_basis_to_array(sd_basis)

    colex_to_lex = np.full(n_sd, -1, dtype=np.int64)
    for i in range(n_sd):
        colex_idx = _sd_to_index(sd_array[i], n_el, binom)
        colex_to_lex[colex_idx] = i

    assert np.all(colex_to_lex >= 0), "Colex index mapping has gaps"
    return colex_to_lex


def prepare_spatial_targets(
    spatial_targets: dict,
    n_spatial: int,
) -> tuple:
    """
    Convert spatial_targets dict to Numba-compatible flat arrays.

    Merges forward and reverse keys: for each (a,b), the stored targets
    are the union of spatial_targets[(a,b)] and spatial_targets[(b,a)].
    This eliminates the need for a swap loop in the assembly kernel.

    Returns (offsets, lengths, values) where:
      offsets[a, b] = start index in values for key (a, b), or -1
      lengths[a, b] = number of target pairs for key (a, b)
      values[i] = (sp_lo, sp_hi) target pair
    """
    # Build merged dict: for each (a,b), union with (b,a)
    merged: dict = {}
    for (a, b), targets in spatial_targets.items():
        if 0 <= a < n_spatial and 0 <= b < n_spatial:
            key = (a, b)
            if key not in merged:
                merged[key] = set()
            merged[key] |= targets

            # Also merge into reverse key
            rev = (b, a)
            if rev not in merged:
                merged[rev] = set()
            merged[rev] |= targets

    offsets = np.full((n_spatial, n_spatial), -1, dtype=np.int64)
    lengths = np.zeros((n_spatial, n_spatial), dtype=np.int64)

    total = sum(len(v) for v in merged.values())
    values = np.empty((max(total, 1), 2), dtype=np.int32)

    pos = 0
    for (a, b), targets in merged.items():
        offsets[a, b] = pos
        lengths[a, b] = len(targets)
        for (lo, hi) in targets:
            values[pos, 0] = lo
            values[pos, 1] = hi
            pos += 1

    return offsets, lengths, values[:max(pos, 1)]


# ======================================================================
# Numba JIT kernels
# ======================================================================

if NUMBA_AVAILABLE:

    @njit(cache=True)
    def _precompute_binomial(max_n: int, max_k: int) -> np.ndarray:
        """Precompute C(n, k) table."""
        binom = np.zeros((max_n + 1, max_k + 1), dtype=np.int64)
        for n in range(max_n + 1):
            binom[n, 0] = 1
            for k in range(1, min(n, max_k) + 1):
                binom[n, k] = binom[n - 1, k - 1] + binom[n - 1, k]
        return binom

    @njit(cache=True)
    def _sd_to_index(occ: np.ndarray, n_el: int, binom: np.ndarray) -> int:
        """Combinatorial number system: sorted tuple → sequential index."""
        idx = 0
        for i in range(n_el):
            idx += binom[occ[i], i + 1]
        return idx

    @njit(cache=True)
    def _sort_small(arr: np.ndarray, n: int) -> None:
        """In-place insertion sort for n <= ~10."""
        for i in range(1, n):
            key = arr[i]
            j = i - 1
            while j >= 0 and arr[j] > key:
                arr[j + 1] = arr[j]
                j -= 1
            arr[j + 1] = key

    @njit(cache=True)
    def _phase_single(sd: np.ndarray, n_el: int, kp: int, r: int) -> int:
        """Fermionic sign for sd[kp] → r."""
        p = sd[kp]
        kr = 0
        for i in range(n_el):
            if sd[i] < r and sd[i] != p:
                kr += 1
        return 1 - 2 * ((kp + kr) & 1)

    @njit(cache=True)
    def _phase_double(
        sd: np.ndarray, n_el: int, kp: int, kq: int, r: int, s: int
    ) -> int:
        """Fermionic sign for sd[kp],sd[kq] → r,s via two sequential excitations."""
        p = sd[kp]
        q = sd[kq]

        # First: remove p, insert r
        kr = 0
        for i in range(n_el):
            if sd[i] != p and sd[i] < r:
                kr += 1
        phase1 = 1 - 2 * ((kp + kr) & 1)

        # Build sorted intermediate: sd - {p} + {r}
        rem = np.empty(n_el - 1, dtype=np.int32)
        j = 0
        for i in range(n_el):
            if sd[i] != p:
                rem[j] = sd[i]
                j += 1
        inter = np.empty(n_el, dtype=np.int32)
        inserted = False
        j = 0
        for i in range(n_el - 1):
            if not inserted and r < rem[i]:
                inter[j] = r
                j += 1
                inserted = True
            inter[j] = rem[i]
            j += 1
        if not inserted:
            inter[j] = r

        # Second: remove q from inter, insert s
        kq_new = 0
        for i in range(n_el):
            if inter[i] == q:
                kq_new = i
                break
        ks = 0
        for i in range(n_el):
            if inter[i] != q and inter[i] < s:
                ks += 1
        phase2 = 1 - 2 * ((kq_new + ks) & 1)

        return phase1 * phase2

    @njit(cache=True)
    def _assemble_kernel(
        sd_array: np.ndarray,       # (n_sd, n_el) int32
        H1: np.ndarray,             # (n_spatial, n_spatial) float64
        eri_4d: np.ndarray,         # (n_spatial,)*4 float64
        h1_diag: np.ndarray,        # (n_spatial,) float64
        tgt_off: np.ndarray,        # (n_spatial, n_spatial) int64
        tgt_len: np.ndarray,        # (n_spatial, n_spatial) int64
        tgt_val: np.ndarray,        # (n_targets, 2) int32
        binom: np.ndarray,          # (n_sp+1, n_el+1) int64
        c2l: np.ndarray,            # (n_sd,) colex→lex mapping
        n_sd: int,
        n_el: int,
        n_sp: int,
        n_spatial: int,
        threshold: float,
        max_off: int,
    ) -> tuple:
        """
        Monolithic assembly: diagonal + singles + doubles, all in one SD pass.

        Returns
        -------
        diag_idx  : int64[n_diag]
        diag_val  : float64[n_diag]
        off_row   : int64[n_off]
        off_col   : int64[n_off]
        off_val   : float64[n_off]
        """
        diag_idx = np.empty(n_sd, dtype=np.int64)
        diag_val = np.empty(n_sd, dtype=np.float64)
        nd = 0

        off_row = np.empty(max_off, dtype=np.int64)
        off_col = np.empty(max_off, dtype=np.int64)
        off_val = np.empty(max_off, dtype=np.float64)
        no = 0

        new_sd = np.empty(n_el, dtype=np.int32)

        for I in range(n_sd):
            sd_I = sd_array[I]

            # ============ DIAGONAL ============
            hd = 0.0
            for k in range(n_el):
                hd += h1_diag[sd_I[k] >> 1]
            for i in range(n_el):
                pi = sd_I[i]
                spi = pi >> 1
                sigi = pi & 1
                for j in range(i + 1, n_el):
                    pj = sd_I[j]
                    spj = pj >> 1
                    hd += eri_4d[spi, spj, spi, spj]
                    if (pj & 1) == sigi:
                        hd -= eri_4d[spi, spj, spj, spi]
            if abs(hd) >= threshold:
                diag_idx[nd] = I
                diag_val[nd] = hd
                nd += 1

            # ============ SINGLES ============
            for kp in range(n_el):
                p = sd_I[kp]
                sp_p = p >> 1
                sig_p = p & 1

                for r in range(n_sp):
                    # Spin conservation
                    if (r & 1) != sig_p:
                        continue
                    # Skip occupied
                    occ = False
                    for k in range(n_el):
                        if sd_I[k] == r:
                            occ = True
                            break
                    if occ:
                        continue

                    sp_r = r >> 1

                    # Matrix element: h1 + sum_q [J(pq,rq) - K(pq,qr)]
                    me = H1[sp_p, sp_r]
                    for k in range(n_el):
                        qq = sd_I[k]
                        if qq == p:
                            continue
                        sp_q = qq >> 1
                        me += eri_4d[sp_p, sp_q, sp_r, sp_q]
                        if (qq & 1) == sig_p:
                            me -= eri_4d[sp_p, sp_q, sp_q, sp_r]

                    if abs(me) < threshold:
                        continue

                    # Target SD
                    for k in range(n_el):
                        new_sd[k] = sd_I[k]
                    new_sd[kp] = r
                    _sort_small(new_sd, n_el)

                    colex_J = _sd_to_index(new_sd, n_el, binom)
                    if colex_J < 0 or colex_J >= n_sd:
                        continue
                    J = c2l[colex_J]
                    if J <= I:
                        continue

                    if no >= max_off:
                        # Return what we have; caller will retry
                        return (diag_idx[:nd], diag_val[:nd],
                                off_row[:no], off_col[:no], off_val[:no])

                    phase = _phase_single(sd_I, n_el, kp, r)
                    off_row[no] = I
                    off_col[no] = J
                    off_val[no] = phase * me
                    no += 1

            # ============ DOUBLES ============
            for kp in range(n_el):
                p = sd_I[kp]
                sp_p = p >> 1
                sig_p = p & 1

                for kq in range(kp + 1, n_el):
                    q = sd_I[kq]
                    sp_q = q >> 1
                    sig_q = q & 1

                    # Targets already merged (forward + reverse) in prep
                    start = tgt_off[sp_p, sp_q]
                    if start < 0:
                        continue
                    length = tgt_len[sp_p, sp_q]

                    for t in range(start, start + length):
                        sp_lo = tgt_val[t, 0]
                        sp_hi = tgt_val[t, 1]

                        if sp_lo < sp_hi:
                            n_spin = 4
                        else:
                            n_spin = 1

                        for si in range(n_spin):
                            if sp_lo < sp_hi:
                                r = (sp_lo << 1) | (si >> 1)
                                s = (sp_hi << 1) | (si & 1)
                            else:
                                r = sp_lo << 1
                                s = (sp_lo << 1) | 1

                            # Occupied check
                            r_occ = False
                            s_occ = False
                            for k in range(n_el):
                                if sd_I[k] == r:
                                    r_occ = True
                                if sd_I[k] == s:
                                    s_occ = True
                            if r_occ or s_occ:
                                continue

                            # Spin conservation
                            sig_r = r & 1
                            sig_s = s & 1
                            if sig_p + sig_q != sig_r + sig_s:
                                continue

                            sp_r = r >> 1
                            sp_s = s >> 1

                            me = 0.0
                            if sig_p == sig_r and sig_q == sig_s:
                                me += eri_4d[sp_p, sp_q, sp_r, sp_s]
                            if sig_p == sig_s and sig_q == sig_r:
                                me -= eri_4d[sp_p, sp_q, sp_s, sp_r]

                            if abs(me) < threshold:
                                continue

                            for k in range(n_el):
                                new_sd[k] = sd_I[k]
                            new_sd[kp] = r
                            new_sd[kq] = s
                            _sort_small(new_sd, n_el)

                            colex_J = _sd_to_index(new_sd, n_el, binom)
                            if colex_J < 0 or colex_J >= n_sd:
                                continue
                            J = c2l[colex_J]
                            if J <= I:
                                continue

                            phase = _phase_double(
                                sd_I, n_el, kp, kq, r, s
                            )
                            if no >= max_off:
                                return (diag_idx[:nd], diag_val[:nd],
                                        off_row[:no], off_col[:no], off_val[:no])

                            off_row[no] = I
                            off_col[no] = J
                            off_val[no] = phase * me
                            no += 1

        return (diag_idx[:nd], diag_val[:nd],
                off_row[:no], off_col[:no], off_val[:no])


# ======================================================================
# Public interface
# ======================================================================

def assemble_hamiltonian_numba(
    sd_basis: list,
    H1: np.ndarray,
    eri_4d: np.ndarray,
    h1_diag: np.ndarray,
    n_sp: int,
    threshold: float,
    spatial_targets: dict,
    n_spatial: int,
) -> tuple:
    """
    Full Numba-accelerated Hamiltonian assembly.

    Returns (diag_idx, diag_val, off_row, off_col, off_val).
    """
    if not NUMBA_AVAILABLE:
        raise RuntimeError("Numba is not installed")

    sd_array = sd_basis_to_array(sd_basis)
    n_sd, n_el = sd_array.shape

    H1 = np.ascontiguousarray(H1, dtype=np.float64)
    eri_4d = np.ascontiguousarray(eri_4d, dtype=np.float64)
    h1_diag = np.ascontiguousarray(h1_diag, dtype=np.float64)

    binom = _precompute_binomial(n_sp, n_el)
    c2l = build_colex_to_lex(sd_basis, n_sp)
    tgt_off, tgt_len, tgt_val = prepare_spatial_targets(
        spatial_targets, n_spatial
    )

    # Estimate max off-diagonal NNZ
    # Typical: ~10-30 per SD. Start generous; retry with 2x if overflow.
    max_off = max(int(n_sd * 50), 100000)

    t0 = time.perf_counter()

    for attempt in range(3):
        result = _assemble_kernel(
            sd_array, H1, eri_4d, h1_diag,
            tgt_off, tgt_len, tgt_val, binom, c2l,
            n_sd, n_el, n_sp, n_spatial, threshold, max_off,
        )
        n_off = len(result[2])
        # Check if kernel hit the overflow guard
        # (it returns early without processing all SDs)
        if n_off < max_off:
            break
        # Overflow: double allocation and retry
        max_off *= 2
        print(f"[Numba] Overflow at {n_off:,} entries, retrying with {max_off:,}")

    dt = time.perf_counter() - t0
    print(
        f"[Numba] Assembly: {dt:.3f}s "
        f"(n_sd={n_sd:,}, n_off={n_off:,})"
    )

    return result


def warmup_jit() -> None:
    """
    Trigger Numba compilation on a tiny problem so the JIT cost
    is paid upfront, not during the real calculation.
    """
    if not NUMBA_AVAILABLE:
        return

    # 2 electrons, 4 spin-orbitals → C(4,2) = 6 SDs
    sd_basis = [
        (0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3),
    ]
    H1 = np.eye(2, dtype=np.float64) * -0.5
    eri_4d = np.zeros((2, 2, 2, 2), dtype=np.float64)
    eri_4d[0, 1, 0, 1] = 0.1
    eri_4d[1, 0, 1, 0] = 0.1
    h1_diag = np.array([-0.5, -0.5], dtype=np.float64)
    spatial_targets = {(0, 1): {(0, 1)}}

    assemble_hamiltonian_numba(
        sd_basis, H1, eri_4d, h1_diag,
        n_sp=4, threshold=1e-14,
        spatial_targets=spatial_targets, n_spatial=2,
    )
