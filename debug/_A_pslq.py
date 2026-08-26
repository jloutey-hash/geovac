import mpmath as mp
mp.mp.dps = 40
A = mp.mpf('0.07905403211681687671461246470412102690878')
# natural genus-0 constants at arguments 2 and 2sqrt2 (region corners of the 2D K0 Yukawa moment)
s2 = mp.sqrt(2); a1=mp.mpf(2); a2=2*s2
consts = {
 '1':mp.mpf(1), 'pi':mp.pi, 'gamma':mp.euler,
 'e^-2':mp.e**-2, 'e^-2s2':mp.e**(-a2),
 'K0(2)':mp.besselk(0,a1),'K0(2s2)':mp.besselk(0,a2),
 'K1(2)':mp.besselk(1,a1),'K1(2s2)':mp.besselk(1,a2),
 'E1(2)':mp.e1(a1),'E1(2s2)':mp.e1(a2),
 'piE^-2':mp.pi*mp.e**-2,'piE^-2s2':mp.pi*mp.e**(-a2),
}
def guarded(name, keys):
    basis=[A]+[consts[k] for k in keys]
    rel=mp.pslq(basis, tol=mp.mpf(10)**-30, maxcoeff=10**6, maxsteps=10**6)
    # decoy: same-magnitude random-ish target
    decoy=A*mp.mpf('1.0000000000031415926535')+mp.mpf('1e-9')
    reld=mp.pslq([decoy]+[consts[k] for k in keys], tol=mp.mpf(10)**-30, maxcoeff=10**6, maxsteps=10**6)
    ht = max(abs(c) for c in rel) if rel else None
    htd= max(abs(c) for c in reld) if reld else None
    print(f"  basis={keys}")
    print(f"    A   rel={rel}  height={ht}")
    print(f"    dec rel={reld} height={htd}")
    if rel and rel[0]!=0 and (htd is None or ht < htd/100):
        print(f"    >>> CANDIDATE (A-coeff {rel[0]}, decoy {'none' if not reld else 'higher'})")
    print(flush=True)

print("A =",mp.nstr(A,38),flush=True)
for keys in [
  ['pi','e^-2','e^-2s2'],
  ['pi','K0(2)','K0(2s2)'],
  ['K0(2)','K1(2)','K0(2s2)','K1(2s2)'],
  ['pi','e^-2','e^-2s2','E1(2)','E1(2s2)'],
  ['piE^-2','piE^-2s2','e^-2','e^-2s2'],
  ['pi','K0(2s2)','K1(2s2)'],
  ['pi','K0(2)','K1(2)'],
]:
    guarded('A',keys)
