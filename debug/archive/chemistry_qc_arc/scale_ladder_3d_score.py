"""Scorer for record SL-3d (docs/predictions/2026-11-01_scale_ladder_3d.md).  Phase 2 only: run after
`debug/scale_ladder_3d_water.py R3d R4d R5d R5Hd` have each appended one row to
debug/data/scale_ladder_3d.jsonl.  Reads only:
  debug/data/scale_ladder_3d.jsonl            new water rows (R2d control from Phase 1; scored rows from Phase 2)
  debug/data/scale_ladder_3d_atom_<R>.json    atom-set 3d exponents and atom changes (Phase 1)
  debug/data/scale_ladder.jsonl               the registered ladder's water rows (record SL)
  debug/data/scale_ladder_r3_diag_atom.json, debug/data/scale_ladder_r3_diag.jsonl   control C1 (R3dB)
It recomputes the controls from the files and refuses to score if they differ from the registered values.

    python debug/scale_ladder_3d_score.py
"""
import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, 'data')
KEY = 'CCSD(T)_fc'

# Registered constants (record SL-3d, "Predictions"); the scorer checks the controls reproduce them.
C1_RATIO, C2_RATIO = 0.963, 1.052          # water change / atom change, controls C1 (R3dB) and C2 (R2d)
HALF_WIDTH = 0.156                          # 3 x max |ratio - 1| over the controls (3 x 0.052)
BAND = (1 - HALF_WIDTH, 1 + HALF_WIDTH)     # SL3d-T: water change / atom change inside this band
CCPVDZ_FC = -76.241025                      # registry p62_ccpvdz_fc
DZP_FC = -76.281739                         # registry p62_dzp_fc
SCORED = {'R3d': ('R3', 'R3'), 'R4d': ('R4', 'R4'), 'R5d': ('R5', 'R5'), 'R5Hd': ('R5H', 'R5')}
#          new row   (old ladder row, atom rung)


def rows(path, key='rung'):
    out = {}
    for line in open(path):
        r = json.loads(line)
        if r[key] in out:
            sys.exit('VOID: two rows for %s in %s; the record scores the first run only' % (r[key], path))
        out[r[key]] = r
    return out


def main():
    new = rows(os.path.join(DATA, 'scale_ladder_3d.jsonl'))
    old = rows(os.path.join(DATA, 'scale_ladder.jsonl'))
    old['R2'] = old['R2rep']
    atom = {R: json.load(open(os.path.join(DATA, 'scale_ladder_3d_atom_%s.json' % R)))
            for R in ('R2', 'R3', 'R4', 'R5')}
    missing = [n for n in ['R2d'] + list(SCORED) if n not in new]
    if missing:
        sys.exit('not scorable yet: missing rows %s' % missing)

    # --- controls (must reproduce the registered values) ---
    d = json.load(open(os.path.join(DATA, 'scale_ladder_r3_diag_atom.json')))['energies']
    diag = rows(os.path.join(DATA, 'scale_ladder_r3_diag.jsonl'))
    a1 = 1e3 * (d['R3dB'][KEY] - d['R3'][KEY])
    w1 = 1e3 * (diag['R3dB'][KEY] - old['R3'][KEY])
    a2 = atom['R2']['atom_change_mHa'][KEY]
    w2 = 1e3 * (new['R2d'][KEY] - old['R2'][KEY])
    for got, reg, name in ((w1 / a1, C1_RATIO, 'C1'), (w2 / a2, C2_RATIO, 'C2')):
        if abs(got - reg) > 5e-4:
            sys.exit('VOID: control %s ratio %.4f differs from the registered %.3f' % (name, got, reg))
    print('controls: C1 R3dB atom %.3f water %.3f ratio %.3f | C2 R2d atom %.3f water %.3f ratio %.3f'
          % (a1, w1, w1 / a1, a2, w2, w2 / a2))

    # --- each scored row must use the registered s/p exponents and the atom-set 3d ---
    for n, (o, R) in SCORED.items():
        r = new[n]
        if (max(abs(x - y) for x, y in zip(r['o_k_1s_1s_2s_2p_3p'], atom[R]['s_p_exponents'])) > 1e-12
                or abs(r['k3d'] - atom[R]['k3d_new']) > 1e-12):
            sys.exit('VOID: row %s does not use the registered exponents' % n)

    # --- SL3d-T: transfer of the 3d re-setting, per scored rung ---
    print('\nSL3d-T  band %.3f..%.3f' % BAND)
    print('%-5s %10s %10s %7s %9s %9s %9s  %s' % ('rung', 'atom mHa', 'water mHa', 'ratio', 'RHF part',
                                                  'corr part', 'atom corr', 'verdict'))
    t_held, rhf_neg, corr_lt1 = {}, {}, {}
    for n, (o, R) in SCORED.items():
        A = atom[R]['atom_change_mHa'][KEY]
        Ac = A - atom[R]['atom_change_mHa']['UHF']
        W = 1e3 * (new[n][KEY] - old[o][KEY])
        Wh = 1e3 * (new[n]['RHF'] - old[o]['RHF'])
        ratio = W / A
        t_held[n] = BAND[0] <= ratio <= BAND[1]
        rhf_neg[n] = Wh < 0
        corr_lt1[n] = (W - Wh) / Ac < 1
        print('%-5s %10.3f %10.3f %7.3f %9.3f %9.3f %9.3f  %s' % (n, A, W, ratio, Wh, W - Wh, Ac,
                                                                 'HELD' if t_held[n] else 'MISSED'))
    nt = sum(t_held.values())

    # --- SL3d-M: R2 -> R3 on water under the new rule ---
    step = 1e3 * (new['R3d'][KEY] - new['R2d'][KEY])
    m_held = step > 0
    print('\nSL3d-M  W(R3d) - W(R2d) = %+.3f mHa (predicted > 0, non-monotone): %s'
          % (step, 'HELD' if m_held else 'MISSED'))

    # --- reported ---
    print('\nSL3d-H (reported) water RHF part < 0 on each scored rung: %s' % rhf_neg)
    print('SL3d-C (reported) water corr part / atom corr part < 1 on each scored rung: %s' % corr_lt1)
    for n in ('R5d', 'R5Hd'):
        print('SL3d-R (reported) W(%s) = %.6f: %+.2f mHa vs cc-pVDZ, %+.2f mHa vs Slater DZP' % (
            n, new[n][KEY], 1e3 * (new[n][KEY] - CCPVDZ_FC), 1e3 * (new[n][KEY] - DZP_FC)))
    steps = ['R2d', 'R3d', 'R4d', 'R5d']
    print('new ladder fc-CCSD(T):', ', '.join('%s %.6f' % (s, new[s][KEY]) for s in steps + ['R5Hd']))

    # --- gate ---
    if nt == 4 and m_held:
        gate = 'GO'
    elif nt == 4 or nt == 3:
        gate = 'BORDERLINE'
    else:
        gate = 'STOP'
    print('\nGate: %s  (SL3d-T held on %d of 4 scored rungs; SL3d-M %s)' % (gate, nt, 'held' if m_held else 'missed'))


if __name__ == '__main__':
    main()
