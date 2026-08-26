"""
How much does the same-spin double-excitation sign defect move the *published*
balanced-LiH headlines?  (energy error at the minimum, R_eq, omega_e, tilt)

The E_coupled convention carries an R-independent core double-count; it is
calibrated here against the banked full-BO curve so absolute energy errors can be
quoted in the published convention.
"""
from __future__ import annotations
import json
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import numpy as np
from davidson_pes import fit_shape, R_TRUE

E_EXACT = -8.0706          # LiH BO total (pk_amplification_lih.json 'reference')

# banked full-BO totals used to calibrate the core-double-count offset per max_n
BANKED = {
    2: ([2.0, 2.5, 3.0, 3.015, 3.5, 4.0, 5.0],
        [-7.824551, -7.894006, -7.928962, -7.929426, -7.927324, -7.895162, -7.771569]),
    3: ([2.9, 3.015, 3.1, 3.3, 3.5],
        [-8.0493805415, -8.0545240769, -8.0570773414, -8.0591647662, -8.0562794055]),
}


def offset_for(max_n):
    """E_BO = E_coupled + offset. Calibrated at R_true from the banked curve where
    available; for max_n with no banked curve, extrapolate the (very weakly
    max_n-dependent) offset from the validation sweeps."""
    fn = f'debug/data/davidson_pes_n{max_n}_decider.json'
    d = json.load(open(fn))
    row = [r for r in d['rows'] if abs(r['R'] - R_TRUE) < 1e-9][0]
    Ec = row['faithful']['E']
    if max_n in BANKED:
        Rb, Eb = BANKED[max_n]
        i = Rb.index(R_TRUE) if R_TRUE in Rb else None
        if i is not None:
            return Eb[i] - Ec, 'calibrated vs banked BO curve'
    return None, 'no banked BO curve for this max_n'


if __name__ == '__main__':
    print("=" * 100)
    print("Impact of the same-spin double-excitation sign defect on the published headlines")
    print("=" * 100)
    print(f"E_exact(LiH, BO) = {E_EXACT} Ha\n")
    hdr = (f"{'max_n':>5s} {'mode':>10s} {'E_min(BO)':>13s} {'err Ha':>10s} {'err %':>7s} "
           f"{'R_eq':>7s} {'R_eq%':>7s} {'omega_e':>8s} {'w_e%':>7s} {'tilt':>10s}")
    print(hdr); print('-' * len(hdr))
    for max_n in (2, 3, 4):
        fn = f'debug/data/davidson_pes_n{max_n}_decider.json'
        if not os.path.exists(fn):
            continue
        d = json.load(open(fn))
        off, how = offset_for(max_n)
        for tag in ('faithful', 'corrected'):
            pts = [(r['R'], r[tag]['E']) for r in d['rows'] if tag in r]
            if len(pts) < 5:
                continue
            Rs = [p[0] for p in pts]; Es = [p[1] for p in pts]
            f = fit_shape(Rs, Es, order=4)
            om = f['own_min']
            if not om:
                continue
            # energy at the fitted own minimum
            c = np.polyfit(np.array(Rs) - R_TRUE, Es, 4)
            Emin = float(np.poly1d(c)(om['R_eq'] - R_TRUE))
            if off is not None:
                Ebo = Emin + off
                err = Ebo - E_EXACT
                print(f"{max_n:5d} {tag:>10s} {Ebo:+13.6f} {err:+10.5f} "
                      f"{abs(err)/abs(E_EXACT)*100:7.3f} {om['R_eq']:7.4f} "
                      f"{om['R_eq_err_pct']:+7.2f} {om['omega_e_cm1']:8.0f} "
                      f"{om['omega_e_err_pct']:+7.1f} {f['tilt_at_Rtrue']:+10.6f}")
            else:
                print(f"{max_n:5d} {tag:>10s} {'(no BO offset)':>13s} {'':>10s} {'':>7s} "
                      f"{om['R_eq']:7.4f} {om['R_eq_err_pct']:+7.2f} "
                      f"{om['omega_e_cm1']:8.0f} {om['omega_e_err_pct']:+7.1f} "
                      f"{f['tilt_at_Rtrue']:+10.6f}")
        if off is not None:
            print(f"      (core double-count offset E_BO - E_coupled = {off:+.6f}; {how})")
    print("\nsign-defect shift in E_coupled at R_true:")
    for max_n in (2, 3, 4):
        fn = f'debug/data/davidson_pes_n{max_n}_decider.json'
        if not os.path.exists(fn):
            continue
        d = json.load(open(fn))
        rr = [r for r in d['rows'] if abs(r['R'] - R_TRUE) < 1e-9]
        if rr and 'corrected' in rr[0]:
            r = rr[0]
            print(f"  max_n={max_n}: corrected - faithful = "
                  f"{r['corrected']['E'] - r['faithful']['E']:+.6f} Ha")
