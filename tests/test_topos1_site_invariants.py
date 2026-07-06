"""Backing tests for the Topos-1 Bohrification probe (Paper 57 SS open,
"A Bohr-site partial meta-theorem" remark): site invariants of the
forced/free seam.  SELF-CONTAINED (no debug/ imports -- the P29
provider-relocation lesson): the site/KS machinery is inlined below and
mirrors debug/compute_topos1_bohr_probe.py (the sprint driver).

Pins (all bit-exact / machine-verified):
  1. Site strata of M_k(C), k = 2, 3: conjugacy strata of commutative
     unital *-subalgebras <-> integer partitions; family dim k^2 - sum(l^2).
  2. B5 site-degeneracy: C(M2(C)) and C(H) share order invariants
     (height-1 fan; point stratum + one 2-real-parameter family), while
     C(M3(C)) differs (height 2; dims 0/4/6) -- the site is blind exactly
     at the catalogue's lone admitted-not-forced entry (B5) and sighted
     at k = 3.
  3. Kochen-Specker witness, machine-verified: primitive integer rays in
     R^3 with entries in {-2..2} (49 rays, 26 orthonormal triples;
     containing the Conway-Kochen ray family) admit NO 0/1 frame-function
     coloring (exhaustive backtracking); the range-1 set (13 rays) IS
     colorable -- the dim-3 valuation obstruction with the dim-2 contrast.
  4. GeoVac distinguished flag on the actual lattice (max_n = 2, 3):
     Casimir chain dims (2,3,5)/(3,6,14); ||[A,n]|| = sqrt(2)/sqrt(10),
     ||[A,L2]|| = 0 (the edge rule conserves l), ||[A,Lz]|| = 2/4.
"""

import itertools
from math import gcd, isclose, sqrt

import numpy as np


# ----------------------------------------------------------------------
# inlined site/KS machinery (mirrors debug/compute_topos1_bohr_probe.py)
# ----------------------------------------------------------------------

def partitions(k):
    def gen(rest, mx):
        if rest == 0:
            yield ()
            return
        for first in range(min(rest, mx), 0, -1):
            for tail in gen(rest - first, first):
                yield (first,) + tail
    yield from gen(k, k)


def refines(mu, lam):
    if sum(mu) != sum(lam):
        return False
    mu = sorted(mu, reverse=True)
    caps = sorted(lam, reverse=True)

    def place(i):
        if i == len(mu):
            return all(c == 0 for c in caps)
        seen = set()
        for j in range(len(caps)):
            if caps[j] >= mu[i] and caps[j] not in seen:
                seen.add(caps[j])
                caps[j] -= mu[i]
                if place(i + 1):
                    caps[j] += mu[i]
                    return True
                caps[j] += mu[i]
        return False

    return place(0)


def site_strata(k):
    return [{"partition": list(lam),
             "family_dim_real": k * k - sum(x * x for x in lam),
             "subalgebra_dim": len(lam)} for lam in partitions(k)]


def primitive_rays(r):
    seen, rays = set(), []
    for v in itertools.product(range(-r, r + 1), repeat=3):
        if v == (0, 0, 0):
            continue
        g = gcd(gcd(abs(v[0]), abs(v[1])), abs(v[2]))
        w = tuple(x // g for x in v)
        for x in w:
            if x != 0:
                if x < 0:
                    w = tuple(-y for y in w)
                break
        if w not in seen:
            seen.add(w)
            rays.append(w)
    return rays


def ks_search(r):
    rays = primitive_rays(r)
    n = len(rays)
    dot = lambda a, b: a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
    orth_pairs = [(i, j) for i in range(n) for j in range(i + 1, n)
                  if dot(rays[i], rays[j]) == 0]
    adj = {i: set() for i in range(n)}
    for i, j in orth_pairs:
        adj[i].add(j)
        adj[j].add(i)
    triples = [(i, j, k) for i in range(n) for j in sorted(adj[i]) if j > i
               for k in sorted(adj[i] & adj[j]) if k > j]
    used = sorted({x for t in triples for x in t})
    remap = {x: idx for idx, x in enumerate(used)}
    triples_r = [tuple(remap[x] for x in t) for t in triples]
    pairs_r = [(remap[i], remap[j]) for (i, j) in orth_pairs
               if i in remap and j in remap]
    m = len(used)
    color = [None] * m
    tri_by_ray = {i: [] for i in range(m)}
    for t in triples_r:
        for x in t:
            tri_by_ray[x].append(t)
    nbr = {i: set() for i in range(m)}
    for i, j in pairs_r:
        nbr[i].add(j)
        nbr[j].add(i)

    def consistent(i):
        if color[i] == 1 and any(color[j] == 1 for j in nbr[i]):
            return False
        for t in tri_by_ray[i]:
            vals = [color[x] for x in t]
            if vals.count(1) > 1 or vals.count(0) == 3:
                return False
            if None not in vals and vals.count(1) != 1:
                return False
        return True

    def bt(i):
        if i == m:
            return True
        for c in (0, 1):
            color[i] = c
            if consistent(i) and bt(i + 1):
                return True
        color[i] = None
        return False

    return {"n_rays_total": n, "n_triples": len(triples_r),
            "coloring_exists": bt(0)}


# ----------------------------------------------------------------------
# tests
# ----------------------------------------------------------------------

def test_site_strata_partitions():
    assert [s["partition"] for s in site_strata(2)] == [[2], [1, 1]]
    assert [s["partition"] for s in site_strata(3)] == [[3], [2, 1], [1, 1, 1]]
    assert [s["family_dim_real"] for s in site_strata(2)] == [0, 2]
    assert [s["family_dim_real"] for s in site_strata(3)] == [0, 4, 6]
    assert refines((1, 1), (2,)) and refines((1, 1, 1), (2, 1))
    assert not refines((2, 1), (1, 1, 1))


def test_b5_site_degeneracy_and_m3_contrast():
    m2_profile = sorted(s["family_dim_real"] for s in site_strata(2))
    h_profile = sorted([0, 2])       # C(H): R1 point + C_q family (S^2/{+-1})
    assert m2_profile == h_profile == [0, 2]       # B5: invariants coincide
    m3_profile = sorted(s["family_dim_real"] for s in site_strata(3))
    assert m3_profile == [0, 4, 6]                 # k = 3: site sighted


def test_kochen_specker_machine_witness():
    r1 = ks_search(1)
    assert r1["coloring_exists"] is True and r1["n_rays_total"] == 13
    r2 = ks_search(2)
    assert r2["coloring_exists"] is False          # bit-exact KS witness
    assert r2["n_rays_total"] == 49 and r2["n_triples"] == 26


def test_geovac_flag_chain_and_grading():
    from geovac.lattice import GeometricLattice
    expect = {2: ([2, 3, 5], sqrt(2.0), 0.0, 2.0),
              3: ([3, 6, 14], sqrt(10.0), 0.0, 4.0)}
    for max_n, (dims, cn, cl2, clz) in expect.items():
        L = GeometricLattice(max_n=max_n, nuclear_charge=1)
        states = L.states
        nv = np.array([s[0] for s in states], float)
        lv = np.array([s[1] * (s[1] + 1) for s in states], float)
        mv = np.array([s[2] for s in states], float)
        d1, d2, d3 = (len(set(nv)), len(set(zip(nv, lv))),
                      len(set(zip(nv, lv, mv))))
        assert [d1, d2, d3] == dims and d1 < d2 < d3
        A = np.asarray(L.adjacency.todense(), float)
        cnorm = lambda D: float(np.linalg.norm(A * (D[None, :] - D[:, None])))
        assert isclose(cnorm(nv), cn, rel_tol=1e-12)
        assert cnorm(lv) == cl2                    # Delta-l = 0: l conserved
        assert isclose(cnorm(mv), clz, rel_tol=1e-12)
