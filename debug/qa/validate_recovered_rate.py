"""Validate the RECOVERED (git-history) genuine rank-2 rate machinery reproduces
c(SU(3)) ~ 4/pi -- i.e. the Paper 40 universality numbers are reproducible from
the original derivation (non-circular), not just hardcoded literals."""
import sys, math
sys.path.insert(0, "debug/qa/_resurrected")
from dirac_triangle_extended_verify import build_A
from sp2_g2_rate_constant import run_panel, fit_models

a2 = build_A(2)  # SU(3) = A_2
panel = [16, 25, 36, 49, 64, 100, 144, 200, 300]
rows, meta = run_panel(a2, "SU(3)/A_2 (recovered)", panel)
lams = [r["lambda"] for r in rows]
gs = [r["gamma"] for r in rows]
fit = fit_models(lams, gs, exclude_below_L=5.0, label="SU3-recovered")
print("\n=== VERDICT ===")
print(f"  4/pi              = {4/math.pi:.4f}")
print(f"  16/pi^2 (decoy B) = {16/math.pi**2:.4f}")
if "two_param" in fit:
    print(f"  c (2-param SW a)  = {fit['two_param']['a']:.4f}")
print(f"  paper table SU(3) = 1.243")
