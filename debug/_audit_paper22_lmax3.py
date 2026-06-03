"""Independent recomputation of D(l_max=3) angular ERI density for Paper 22 audit."""
from sympy.physics.wigner import wigner_3j

l_max = 3
orbs = []
for l in range(l_max + 1):
    for m in range(-l, l + 1):
        orbs.append((l, m))
N = len(orbs)
print("N_orb =", N, "N_orb^4 =", N**4)

nz = 0
for (l1, m1) in orbs:
    for (l2, m2) in orbs:
        for (l3, m3) in orbs:
            for (l4, m4) in orbs:
                if m1 + m2 != m3 + m4:
                    continue
                kmin = max(abs(l1 - l3), abs(l2 - l4))
                kmax = min(l1 + l3, l2 + l4)
                found = False
                for k in range(kmin, kmax + 1):
                    if (l1 + l3 + k) % 2 != 0:
                        continue
                    if (l2 + l4 + k) % 2 != 0:
                        continue
                    q1 = m1 - m3
                    q2 = m2 - m4
                    a1 = wigner_3j(l1, k, l3, -m1, q1, m3)
                    if a1 == 0:
                        continue
                    a2 = wigner_3j(l2, k, l4, -m2, q2, m4)
                    if a2 == 0:
                        continue
                    b1 = wigner_3j(l1, k, l3, 0, 0, 0)
                    if b1 == 0:
                        continue
                    b2 = wigner_3j(l2, k, l4, 0, 0, 0)
                    if b2 == 0:
                        continue
                    found = True
                    break
                if found:
                    nz += 1

total = N**4
density = nz / total * 100
print("nonzero:", nz)
print("total:", total)
print("density: {:.4f}%".format(density))
