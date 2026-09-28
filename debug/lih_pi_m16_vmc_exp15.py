"""Run lih_pi_m16_vmc.py against the exp gamma=1.5 trunc-1e-3 energy npz (G2 twin for E_R12=-8.049045)."""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import lih_pi_m16_vmc as V
V.ENERGY_NPZ = os.path.join(V.DATA, "lih_pi_m16_energy_exp1.5.npz")
if __name__ == "__main__":
    V.main()
