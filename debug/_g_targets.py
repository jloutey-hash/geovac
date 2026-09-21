"""Per-component VMC targets for g = <G|H|G>, G=(F-Fbar)Phi0, geminal f=exp(-0.5 r).
  g_T   = <1/2 sum_i (chi^2 v_i^2 + 2 chi gradchi_i.v_i + |gradchi_i|^2)>
  g_Vne = <chi^2 V_ne>,   g_Vee = <chi^2 V_ee>,   g = g_T+g_Vne+g_Vee  (no V_NN).
Also breaks g_T into its 3 sub-terms (drift^2 / cross / |gradF|^2)."""
import numpy as np
from lih_r12ci_sigma2_mc import sample_block
from lih_r12ci_vmc import geom_quantities, F_and_gradF

GAM = 0.5
cfg = dict(nw=12000, burn=3000, nsnap=16, thin=45)
up, au = sample_block(seed=31, **cfg); dn, ad = sample_block(seed=32, **cfg)
ns = len(up)
# global Fbar
Fall = []
POS = []
for s in range(ns):
    pos = np.concatenate([up[s], dn[s]], axis=1); POS.append(pos)
    Fall.append(F_and_gradF(pos, GAM, kind='exp')[0])
Fbar = np.concatenate(Fall).mean()

gT1, gT2, gT3, gVne, gVee, sig2 = [], [], [], [], [], []
for s in range(ns):
    pos = POS[s]
    v, v2, vne, vee = geom_quantities(pos)
    F, gradF = F_and_gradF(pos, GAM, kind='exp')
    chi = F - Fbar
    gcv = (gradF * v).sum(axis=(1, 2)); gc2 = (gradF ** 2).sum(axis=(1, 2))
    gT1.append((0.5 * chi ** 2 * v2).mean())
    gT2.append((chi * gcv).mean())
    gT3.append((0.5 * gc2).mean())
    gVne.append((chi ** 2 * vne).mean()); gVee.append((chi ** 2 * vee).mean())
    sig2.append((chi ** 2).mean())
sem = lambda x: np.std(x, ddof=1) / np.sqrt(len(x))
m = lambda x: np.mean(x)
gT = m(gT1) + m(gT2) + m(gT3)
print(f"Fbar={Fbar:.5f}  acc {au:.2f}/{ad:.2f}  sigma2={m(sig2):.5f}+/-{sem(sig2):.1e}")
print(f"g_T   = {gT:+.5f}   [drift2={m(gT1):+.5f}+/-{sem(gT1):.0e}  cross={m(gT2):+.5f}+/-{sem(gT2):.0e}  |gradF|2={m(gT3):+.5f}+/-{sem(gT3):.0e}]")
print(f"g_Vne = {m(gVne):+.5f} +/- {sem(gVne):.1e}")
print(f"g_Vee = {m(gVee):+.5f} +/- {sem(gVee):.1e}")
print(f"g = g_T+g_Vne+g_Vee = {gT+m(gVne)+m(gVee):+.5f}   (target -0.9303)")
