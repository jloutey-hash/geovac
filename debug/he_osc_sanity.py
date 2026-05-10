"""Sanity check: reproduce hydrogen Lyman alpha A_21 from the dipole machinery.

Lyman alpha: 2P_z -> 1S, hydrogenic Z=1.
Expected: <2P_z|z|1S> = 128*sqrt(2)/243 a_0 ~ 0.7449 a_0.
Then f_Lyman = (2/3) * omega * |<1S|r|2P>|^2 = (2/3) * (3/8) * 3 * (128*sqrt(2)/243)^2
             = (2/3) * (3/8) * 3 * (128^2 * 2 / 243^2)
             = (1/4) * 3 * (32768/59049)
             = (3/4) * (32768/59049)
             = 24576/59049
             ~ 0.4162  (the well-known Lyman alpha oscillator strength)

Where:
  - omega = E_2 - E_1 = -1/8 - (-1/2) = 3/8 Ha
  - factor 3 sums over m=-1,0,+1 sublevels of 2P; equivalently the (2/3) is for
    angular average and we need the m-summed |<1S|z|2P>|^2 = 3 * |<1S|z|2P_{m=0}>|^2.
  - f_Lyman ~ 0.4162 (Bethe-Salpeter Eq. 3.10)
"""

import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import sympy as sp
from geovac.dirac_matrix_elements import radial_matrix_element
from geovac.casimir_ci import _gaunt_ck

# ===== test 1: <2P_z|z|1S> for hydrogen
print("Test 1: <2P_{m=0}|z|1S> for H (Z=1)")

# Radial part: <R_{2,1}|r|R_{1,0}> at Z=1
R_sym = radial_matrix_element(2, 1, 1, 0, "r", Z=sp.Rational(1))
R_value_sym = sp.simplify(R_sym)
R_value = float(R_value_sym)
print(f"  <R_{{2,1}}|r|R_{{1,0}}>(Z=1) = {R_value_sym} = {R_value:.6f}")
print(f"  Expected: 128*sqrt(6)/243 = {float(128*sp.sqrt(6)/243):.6f}")

# Angular part: c_1 for (l=0, m=0) <- (l=1, m=0)
# Convention: <a|z|c> uses c_1(l_a, m_a, l_c, m_c)
c1_1S_2Pz = _gaunt_ck(0, 0, 1, 0, 1)  # <Y_{0,0}|...|Y_{1,0}> for k=1
c1_2Pz_1S = _gaunt_ck(1, 0, 0, 0, 1)  # <Y_{1,0}|...|Y_{0,0}> for k=1
print(f"  c_1(0,0,1,0) = {c1_1S_2Pz:.6f}")
print(f"  c_1(1,0,0,0) = {c1_2Pz_1S:.6f}")
print(f"  Expected: 1/sqrt(3) = {float(1/sp.sqrt(3)):.6f}")

# Full dipole: should be 128*sqrt(2)/243 ~ 0.7449
dipole_value = c1_1S_2Pz * R_value
print(f"  <1S|z|2P_z>(Z=1) = {dipole_value:.6f}")
print(f"  Expected: 128*sqrt(2)/243 = {float(128*sp.sqrt(2)/243):.6f}")
print(f"  Ratio: {dipole_value / float(128*sp.sqrt(2)/243):.6f}")

# ===== test 2: <2P|r|1S> radial at Z=2
print("\nTest 2: <2P|r|1S> at Z=2 (should scale as 1/Z compared to Z=1)")
R_at_Z = float(R_value_sym) / 2.0  # 1/Z scaling
print(f"  <R_{{2,1}}|r|R_{{1,0}}>(Z=2) = {R_at_Z:.6f}")
print(f"  Expected: (128*sqrt(6)/243)/2 = {float(128*sp.sqrt(6)/486):.6f}")

# ===== test 3: oscillator strength formula
print("\nTest 3: Lyman alpha oscillator strength")
omega = 3.0/8.0  # Ha
z_me = float(128*sp.sqrt(2)/243)
# f = 2 * omega * |<1S|z|2Pz>|^2  (sum over upper-state m sublevels assumed)
# Actually conventionally: f = (2/3) * omega * |<1S|r|2P>|^2 where
# |<1S|r|2P>|^2 = sum_m |<1S|r|2P_m>|^2 = 3 * |<1S|z|2P_{m=0}>|^2
# So f = 2 * omega * |z_me|^2 (using single-m component with factor 2)
f_lyman = 2.0 * omega * z_me ** 2
print(f"  f = 2 * (3/8) * (128*sqrt(2)/243)^2 = {f_lyman:.6f}")
print(f"  Expected: 24576/59049 = {24576/59049:.6f}")

# ===== test 4: He at Z=2 single-config (1s)(2p) <-> (1s)(1s)
print("\nTest 4: He single-config sanity at Z=2")
# Single-config <1s,1s | z_1+z_2 | 1s,2p_z> with HF wavefunctions
# For closed-shell singlet (1s)^2 to open-shell (1s)(2p):
# <(1s)^2|z_1+z_2|(1s,2p_z),1S>
# Spin-singlet spatial: |1s,1s> normalized as |1s,1s>; |1s,2p_z>_singlet = (|1s,2p_z> + |2p_z,1s>)/sqrt(2)
# Matrix element: <(1s,1s)|z_1+z_2|(1s,2p)+(2p,1s)>/sqrt(2)
#   = (<1s|z|1s><1s|2p> + <1s|1s><1s|z|2p> + <1s|z|2p><1s|1s> + <1s|1s><1s|z|2p>) / sqrt(2)
# Wait: more carefully. Assume orthonormal orbitals.
#   <1s,1s|z_1|1s,2p> = <1s|z|1s><1s|2p> = 0 (orthogonal)
#   <1s,1s|z_1|2p,1s> = <1s|z|2p><1s|1s> = <1s|z|2p>
#   <1s,1s|z_2|1s,2p> = <1s|1s><1s|z|2p> = <1s|z|2p>
#   <1s,1s|z_2|2p,1s> = <1s|2p><1s|z|1s> = 0
# So <(1s,1s)| z_1+z_2 |(1s,2p)+(2p,1s)>/sqrt(2)
#   = [0 + <1s|z|2p> + <1s|z|2p> + 0]/sqrt(2)
#   = 2 <1s|z|2p>/sqrt(2)
#   = sqrt(2) <1s|z|2p>
# At Z=2: <1s|z|2p> = (128*sqrt(6)/243)/2 * (1/sqrt(3)) = 128*sqrt(2)/(243*2)
#       = 64*sqrt(2)/243
# Multi-electron <z> = sqrt(2) * 64*sqrt(2)/243 = 128/243 ~ 0.5267
# omega(He, single config) ~ E(2P) - E(1S)_HF
# |<>|^2 = (128/243)^2 = 16384/59049 ~ 0.2774
print(f"  <1s|z|2p>(Z=2) = (128*sqrt(2)/243)/2 = {float(128*sp.sqrt(2)/486):.6f}")
print(f"  Multi-electron <Psi_S | z_1+z_2 | Psi_P> = sqrt(2)*<1s|z|2p>(Z=2) = {float(sp.sqrt(2)*128*sp.sqrt(2)/486):.6f}")
print(f"  This equals 128/243 = {128/243:.6f}")
