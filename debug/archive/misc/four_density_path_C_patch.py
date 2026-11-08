"""helper: scan with nullity at Imax first, then bisect the minimal t-degree (imports step C)."""
from __future__ import annotations
import json, sys, time
from pathlib import Path
HERE = Path(__file__).parent
sys.path.insert(0, str(HERE))
import four_density_path_C_tmodule as C  # noqa: E402
import four_density_path_A_operator as A  # noqa: E402


def scan2(name, M, start, r_list, Imax, lo, hi, p=C.P0):
    basis, A0, A1 = C.module_matrices(M, start)
    mp_ = C.ModP(A0, A1, p, lo, hi)
    t0 = time.time()
    cols = C.columns(mp_, max(r_list), Imax)
    print(f"[{name}] module dim {len(basis)}; columns built ({time.time()-t0:.0f}s)", flush=True)
    for r in r_list:
        nl, M_, _ = C.nullity(cols, r, Imax, p)
        print(f"[{name}] r={r}: nullity {nl} at t-degree {2*Imax}  ({time.time()-t0:.0f}s)", flush=True)
        if nl:
            lo_, hi_ = 0, Imax
            while lo_ < hi_:
                mid = (lo_ + hi_) // 2
                if C.nullity(cols, r, mid, p)[0]:
                    hi_ = mid
                else:
                    lo_ = mid + 1
            print(f"[{name}] r={r}: minimal t-degree {2*lo_}, nullity there {C.nullity(cols, r, lo_, p)[0]}", flush=True)
            return r, lo_
    return None


if __name__ == "__main__":
    rs = list(range(int(sys.argv[1]), int(sys.argv[2]) + 1, int(sys.argv[3])))
    Imax = int(sys.argv[4])
    print(scan2("path", A.path_module(), ("c0", "h0"), rs, Imax, -Imax - 2, Imax + max(rs) + 4))
