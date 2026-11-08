"""Paper 22 headline diagnostic: does the production ERI angular factor realize
the pair-diagonal-m rule (m_a=m_c & m_b=m_d) or the full-Gaunt rule
(m_a+m_b=m_c+m_d)?

We count nonzero ANGULAR factors (radial-independent) over all (l,m) quartets
up to l_max, exactly as Paper 22's density theorem counts, using the SAME
_gaunt_ck the production two_electron_integral uses.
"""
from itertools import product
from geovac.casimir_ci import _gaunt_ck


def angular_factor_nonzero(a, b, c, d):
    (la, ma), (lb, mb), (lc, mc), (ld, md) = a, b, c, d
    k_min = max(abs(la - lc), abs(lb - ld))
    k_max = min(la + lc, lb + ld)
    tot = 0.0
    for k in range(k_min, k_max + 1):
        if (la + lc + k) % 2 != 0 or (lb + ld + k) % 2 != 0:
            continue
        tot += _gaunt_ck(la, ma, lc, mc, k) * _gaunt_ck(lb, mb, ld, md, k)
    return abs(tot) > 1e-12


for l_max in (2, 3):
    orbs = [(l, m) for l in range(l_max + 1) for m in range(-l, l + 1)]
    Q = len(orbs)
    total = Q ** 4
    nz = 0
    nz_pairdiag = 0   # m_a==m_c and m_b==m_d
    nz_global = 0     # m_a+m_b==m_c+m_d
    nz_violate_global = 0  # nonzero but m_a+m_b != m_c+m_d (would be a bug)
    for a, b, c, d in product(orbs, repeat=4):
        if angular_factor_nonzero(a, b, c, d):
            nz += 1
            pd = (a[1] == c[1]) and (b[1] == d[1])
            gl = (a[1] + b[1]) == (c[1] + d[1])
            if pd:
                nz_pairdiag += 1
            if gl:
                nz_global += 1
            else:
                nz_violate_global += 1
    print(f"l_max={l_max}  Q={Q}  total quartets={total}")
    print(f"  nonzero angular factors : {nz}  ({100*nz/total:.4f}%)")
    print(f"  of which pair-diagonal-m : {nz_pairdiag}  ({100*nz_pairdiag/total:.4f}% of all)")
    print(f"  of which global m-cons   : {nz_global}  ({100*nz_global/total:.4f}% of all)")
    print(f"  nonzero that VIOLATE global m-cons : {nz_violate_global}")
    print(f"  --> pair-diagonal-only? {nz == nz_pairdiag}   full-Gaunt(global)? {nz == nz_global}")
    print()
