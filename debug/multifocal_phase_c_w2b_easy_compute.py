"""Driver for the Phase C-W2b-easy tensor-product propinquity-bound table.

Reproduces §3 of debug/multifocal_phase_c_w2b_easy_memo.md.

Usage:
    python debug/multifocal_phase_c_w2b_easy_compute.py
"""

from __future__ import annotations

import json
import os
import time

from geovac.gh_convergence_tensor import (
    FiveLemmaStatusTensor,
    compute_tensor_propinquity_bound,
    gh_tensor_theorem_statement,
    tensor_convergence_table,
    tensor_L1prime,
    tensor_L2_central_fejer,
    tensor_L4_berezin,
    tensor_L5_assembly,
    verify_joint_convergence_to_zero,
)


def main():
    print("=" * 75)
    print("Phase C-W2b-easy tensor propinquity bound computation")
    print("=" * 75)

    # Theorem statement
    print("\n--- Theorem L5-T (formal statement) ---")
    print(gh_tensor_theorem_statement())

    # L1'-T at small cutoffs (skip prop at large for cost)
    print("\n--- L1'-T: joint operator system dim factorization ---")
    for n_a, n_b in [(2, 2), (2, 3), (3, 3)]:
        out = tensor_L1prime(n_a, n_b)
        print(
            f"  ({n_a},{n_b}): dim(O_a)={out['dim_a']}, dim(O_b)={out['dim_b']}, "
            f"dim_joint={out['dim_joint']}, factorizes={out['dim_factorizes']}, "
            f"prop_a={out['prop_a']}, prop_b={out['prop_b']}, prop_joint={out['prop_joint']}"
        )

    # L2-T closed form
    print("\n--- L2-T: joint cb-norm + Plancherel factorization ---")
    for n_a, n_b in [(2, 2), (2, 3), (3, 3), (3, 4), (4, 4)]:
        out = tensor_L2_central_fejer(n_a, n_b)
        print(
            f"  ({n_a},{n_b}): cb=[{out['cb_norm_a']} * {out['cb_norm_b']} = {out['cb_norm_joint']}],"
            f" cb_factorizes={out['cb_factorizes']}, plancherel_factorizes_at_0={out['plancherel_factorizes_at_0']}"
        )

    # L4-T panel verification at (2,2)
    print("\n--- L4-T: joint Berezin map panel verification ---")
    out_l4 = tensor_L4_berezin(2, 2)
    for k, v in out_l4.items():
        print(f"  {k}: {v}")

    # L5-T table
    print("\n--- L5-T: joint propinquity bound table ---")
    pairs = [(2, 2), (2, 3), (3, 3), (3, 4), (4, 4)]
    t0 = time.time()
    table = tensor_convergence_table(pairs, gamma_prec=15)
    dt = time.time() - t0
    print(f"  computed in {dt:.1f}s")

    print(
        "\n  (n_a, n_b) | gamma_a  | gamma_b  | cb_norm  | bound (fact) | bound (full)"
    )
    print(
        "  -----------+----------+----------+----------+--------------+-------------"
    )
    for k in pairs:
        v = table[k]
        print(
            f"  ({k[0]:>2}, {k[1]:>2})    | {v.gamma_a:.6f} | {v.gamma_b:.6f} |"
            f" {v.cb_norm_joint:.6f} | {v.propinquity_bound:.6f}     |"
            f" {v.propinquity_bound_full:.6f}"
        )

    ok, ratio = verify_joint_convergence_to_zero(table, threshold_ratio=0.7)
    print(
        f"\n  Convergence (ratio Lambda(4,4) / Lambda(2,2)) = {ratio:.4f}; "
        f"passes threshold 0.7 = {ok}"
    )

    # Five-lemma status
    print("\n--- Five-lemma roadmap status (tensor) ---")
    status = FiveLemmaStatusTensor()
    for k, v in status.to_dict().items():
        print(f"  {k}: {v}")

    # Persist data
    out_dir = os.path.join(
        os.path.dirname(os.path.abspath(__file__)), "data"
    )
    os.makedirs(out_dir, exist_ok=True)
    out_path = os.path.join(
        out_dir, "multifocal_phase_c_w2b_easy_table.json"
    )
    payload = {
        "pairs": list(pairs),
        "table": {f"({a},{b})": table[(a, b)].to_dict() for (a, b) in pairs},
        "convergence_ratio": ratio,
        "convergence_passes_threshold_0_7": ok,
        "five_lemma_status": status.to_dict(),
    }
    with open(out_path, "w") as f:
        json.dump(payload, f, indent=2)
    print(f"\nData written to {out_path}")
    print("\n" + "=" * 75)
    print("Done.")


if __name__ == "__main__":
    main()
