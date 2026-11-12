"""Pre-registration controls for record HP (docs/predictions/2026-11-12_heavy_atom_panel.md), run before the
record is committed:
  overlap   -- each Gaussian basis read correctly for C, N, O, F: unit overlap diagonal; ANO-R1's general
               contraction gives an identity overlap per atom (record SZ checked this on Ne only)
  fitbound  -- Gaussian-fit error of the Slater kit on the C and F atoms: 16 Gaussians per s/p primitive and 12
               per polarization function against the registered 12/8, at the registered exponents
"""
import importlib.util, json, os, sys
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
s = importlib.util.spec_from_file_location('hp', os.path.join(HERE, 'heavy_atom_panel.py'))
hp = importlib.util.module_from_spec(s); s.loader.exec_module(hp)

if sys.argv[1] == 'overlap':
    out = {}
    for el in 'CNOF':
        nuc = [(np.zeros(3), hp.ELEM[el]['Z'])]
        for b in ('ccpvdz', 'anor1'):
            fns, T = hp.w.gaussian_basis(hp.BASIS.format(hp.GAUSS[b]), nuc)
            S_ = hp.w.ao_integrals(fns, nuc, T)[0]
            out['%s_%s' % (el, b)] = dict(nbf=int(S_.shape[0]), max_diag_dev=float(np.abs(np.diag(S_) - 1).max()),
                                         max_offdiag=float(np.abs(S_ - np.diag(np.diag(S_))).max()))
    print(json.dumps(out, indent=1))
elif sys.argv[1] == 'fitbound':
    out = {}
    for el in sys.argv[2:]:
        nuc = [(np.zeros(3), hp.ELEM[el]['Z'])]
        five = json.load(open(hp._path('atom_sp_%s.json' % el)))['ks_1s_1s_2s_2p_3p']
        k3d = json.load(open(hp._path('atom_3d_%s.json' % el)))['k3d']
        for ng, ngp in ((12, 8), (16, 12)):
            fns, T = hp.heavy_kit(nuc, el, five, k3d, [1.0], [0, 0, 0], n_gauss=ng, ng_pol=ngp)
            S_, h, g = hp.w.ao_integrals(fns, nuc, T)
            out['%s_%d_%d' % (el, ng, ngp)] = hp.atom_uhf_ccsdt(el, S_, h, g)
        a, b = out['%s_12_8' % el], out['%s_16_12' % el]
        out['%s_shift_mHa' % el] = {k: 1e3 * (b[k] - a[k]) for k in a}
    print(json.dumps(out, indent=1))
