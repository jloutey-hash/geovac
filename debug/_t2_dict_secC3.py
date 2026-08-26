"""C3 (light grid): the FULL physical fibre in position space, both factors as
2D-Euclidean propagator superpositions, with the j0 factor realised as the box smearing.

    J(s,t) = (2/pi) int_0^inf db f_1(b) <f_2>_{|W|}(b),
    <f_2>_{|W|}(b) = (1/2|W|) int_{-|W|}^{|W|} f_2(|b+v|) dv,
    f_i(b) = sqrt(c_i) int_1^inf w(u) u K_1(sqrt(u^2+b^2/c_i))/sqrt(u^2+b^2/c_i) du.
"""
from __future__ import annotations
import sys
import time
import mpmath as mp

sys.path.insert(0, 'debug')
from t2_euclidean_dictionary import J_direct                     # noqa: E402
from _t2_dict_secC import glint, f_prop_fast                     # noqa: E402


def main():
    mp.mp.dps = int(sys.argv[1]) if len(sys.argv) > 1 else 20
    NB = int(sys.argv[2]) if len(sys.argv) > 2 else 16
    NV = int(sys.argv[3]) if len(sys.argv) > 3 else 14
    NU = int(sys.argv[4]) if len(sys.argv) > 4 else 30
    print(f"C3 light: dps={mp.mp.dps}, outer GL={NB}/panel, smear GL={NV}, f-prop GL={NU}/panel")
    for ss, ts in [('0.4', '0.4'), ('0.25', '0.6')]:
        s, t = mp.mpf(ss), mp.mpf(ts)
        bW = s + t
        a1 = 1/mp.sqrt(s*(1 - s))
        Jd = J_direct(s, t)
        Bmax = 26/a1 + 2*bW
        t0 = time.time()

        def smeared(b):
            return glint(lambda v: f_prop_fast(t, abs(b + v), NU), -bW, bW, NV)/(2*bW)

        def integ(b):
            return f_prop_fast(s, b, NU)*smeared(b)
        Jp = (2/mp.pi)*(glint(integ, 0, Bmax/3, NB) + glint(integ, Bmax/3, Bmax, NB))
        rd = abs(Jd - Jp)/abs(Jd)
        print(f"  s={ss} t={ts} |W|={mp.nstr(bW, 6)}   J_momentum={mp.nstr(Jd, 16)}"
              f"   J_position={mp.nstr(Jp, 16)}   rel.d={mp.nstr(rd, 3)}   [{time.time()-t0:.0f}s]",
              flush=True)


if __name__ == '__main__':
    main()
