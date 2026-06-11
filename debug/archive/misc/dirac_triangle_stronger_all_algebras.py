"""
Test the stronger inequality |sigma+rho|^2 >= |lam - lam' + rho|^2
on all four panels (SU(3), SU(4), Sp(2), G_2).

If this holds universally, the proof goes through cleanly.
"""

import os
import sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from dirac_triangle_extended_verify import (
    build_A, build_C2, build_G2, panel_dominant_weights, tensor_product,
)
from dirac_triangle_dominant_norm_check import check_lam_minus_lp_in_decomp


def main():
    print("Testing |sigma+rho|^2 >= |lam-lam'+rho|^2 on all four panels.")
    print()

    panels = [
        ("SU(3) p+q<=5", build_A(2), panel_dominant_weights(2, 5)),
        ("SU(4) a+b+c<=3", build_A(3), panel_dominant_weights(3, 3)),
        ("Sp(2)=C_2 a+b<=3", build_C2(), panel_dominant_weights(2, 3)),
        ("G_2 a+b<=2", build_G2(), panel_dominant_weights(2, 2)),
    ]

    for label, la, panel in panels:
        total_sigma = 0
        stronger_failures = []
        int_failures = 0
        for lam in panel:
            for lam_pr in panel:
                results = check_lam_minus_lp_in_decomp(la, lam, lam_pr)
                for r in results:
                    total_sigma += 1
                    if not r["int"]:
                        int_failures += 1
                    if not r["stronger"]:
                        stronger_failures.append({
                            "lam": lam, "lam_pr": lam_pr, "result": r,
                        })
        print(f"{label}: {total_sigma} sigmas total")
        print(f"  INT failures: {int_failures}")
        print(f"  Stronger ineq failures: {len(stronger_failures)}")
        if stronger_failures:
            for f in stronger_failures[:3]:
                r = f["result"]
                print(f"    lam={f['lam']}, lam'={f['lam_pr']}, sigma={r['sigma']}: "
                      f"|sigma+rho|^2={r['sigma_rho_sq']} < |lam-lam'+rho|^2={r['lam_diff_rho_sq']}")


if __name__ == "__main__":
    main()
