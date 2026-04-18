"""Track RH-C smoke test — verify the Dirac-S^3 Ihara zeta module runs."""
import sys
sys.stdout.reconfigure(encoding='utf-8')
from geovac.ihara_zeta_dirac import dirac_s3_ihara_zeta

for n_max in [1, 2, 3]:
    for rule in ['A', 'B']:
        try:
            r = dirac_s3_ihara_zeta(n_max=n_max, adjacency_rule=rule, factor=False)
            print(
                f"n_max={n_max} rule={rule}: "
                f"V={r['V']} E={r['E']} c={r['c']} beta1={r['r_betti1']} "
                f"q_max={r['q_max']} "
                f"rho={r['spectral_radius']:.4f} "
                f"max_nt={r['max_abs_nontrivial']:.4f} "
                f"sqrt_q={r['sqrt_q_max']:.4f} "
                f"dev={r['ramanujan_deviation']:+.4f} "
                f"Ram={r['ramanujan_verdict']} "
                f"charpoly_int={r['charpoly_integer_coefficient']}"
            )
            # Print per-kappa sizes
            print(f"  per-kappa: {r['per_kappa_sizes']}")
        except Exception as e:
            import traceback
            print(f"n_max={n_max} rule={rule}: ERROR {type(e).__name__}: {e}")
            traceback.print_exc()
