"""Generate Track AK detailed analysis report."""
import sys
sys.path.insert(0, 'c:/Users/jlout/OneDrive/Desktop/Project_Geometric')

from geovac.n_electron_spectral import (
    run_full_analysis,
    _count_s4_22_sigma_channels,
    _count_s4_22_full_channels,
    four_electron_spectral_dimensions,
    tractability_verdict,
    so12_free_spectrum_by_channel,
    compression_analysis,
)
from geovac.n_electron_scope import _so_d_degeneracy

# Generate main report
report = run_full_analysis(
    n_basis_per_angle=5,
    output_file='c:/Users/jlout/OneDrive/Desktop/Project_Geometric/debug/track_ak/spectral_compression_report.txt'
)

# Generate detailed supplement
lines = []
lines.append("=" * 72)
lines.append("TRACK AK: Detailed Supplement")
lines.append("=" * 72)

lines.append("\n\nA. S4 CHANNEL ANALYSIS (sigma only, heteronuclear)")
lines.append("-" * 72)
for lm in range(5):
    n_raw = (lm + 1) ** 4
    n_s4 = _count_s4_22_sigma_channels(lm, homonuclear=False)
    ratio = n_s4 / n_raw if n_raw > 0 else 0
    lines.append(f'  l_max={lm}: raw={n_raw:>6}, S4[2,2]={n_s4:>5}, ratio={ratio:.3f}')

lines.append("\n\nB. S4 CHANNEL ANALYSIS (full m, heteronuclear)")
lines.append("-" * 72)
for lm in range(4):
    n_s4 = _count_s4_22_full_channels(lm, homonuclear=False)
    lines.append(f'  l_max={lm}: S4_full={n_s4:>6}')

lines.append("\n\nC. SPECTRAL DIMENSIONS vs LEVEL 4")
lines.append("-" * 72)
for lm in range(5):
    v = tractability_verdict(lm, n_basis_per_angle=5)
    lines.append(
        f'  l_max={lm}: 4e_spec={v["spectral_dim"]:>6}, '
        f'L4_spec={v["l4_spectral_dim"]:>4}, '
        f'ratio={v["ratio_to_l4_dim"]:>6.1f}x, '
        f'time={v["total_time_s"]:>10.1f}s, '
        f'{v["verdict"]}'
    )

lines.append("\n\nD. SO(12) vs SO(6) DEGENERACY GROWTH")
lines.append("-" * 72)
for nu in range(8):
    g6 = _so_d_degeneracy(6, nu)
    g12 = _so_d_degeneracy(12, nu)
    lines.append(f'  nu={nu}: SO(6)={g6:>6}, SO(12)={g12:>8}, ratio={g12/g6 if g6>0 else 0:.1f}x')

lines.append("\n\nE. SPECTRAL MULTIPLICITY vs SO(12) DEGENERACY (l_max=0)")
lines.append("-" * 72)
spec = so12_free_spectrum_by_channel(0, n_basis_per_angle=5)
for nu, mult in sorted(spec['nu_multiplicities'].items())[:8]:
    mu = nu * (nu + 10) / 2.0
    g12 = _so_d_degeneracy(12, nu)
    lines.append(
        f'  nu={nu}: mu={mu:.1f}, SO(12) degen={g12:>8}, '
        f'spectral mult={mult:>4}, compression={g12/mult if mult>0 else 0:.1f}x'
    )

lines.append("\n\nF. n_basis SENSITIVITY AT l_max=2")
lines.append("-" * 72)
for nb in range(2, 9):
    v = tractability_verdict(2, n_basis_per_angle=nb)
    lines.append(
        f'  n_basis={nb}: dim={v["spectral_dim"]:>6}, '
        f'time={v["total_time_s"]:>10.1f}s, '
        f'L4_ratio={v["ratio_to_l4_dim"]:>6.1f}x, '
        f'{v["verdict"]}'
    )

supplement = "\n".join(lines)
with open('c:/Users/jlout/OneDrive/Desktop/Project_Geometric/debug/track_ak/detailed_supplement.txt', 'w') as f:
    f.write(supplement)

print(supplement)
print("\n\nReports written to debug/track_ak/")
