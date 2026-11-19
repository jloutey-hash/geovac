"""tree_attain_l: module model of a Coulomb tree with l > 0 densities (axial: P_v polynomial in d_z and lap).

Hat normalisation: non-root density rho_v = P_v(d) K0(r); channel-L potential Phi_{w,L} solves
    (theta - L)(theta + L + 1) Phi_{w,L} = - r^2 W_{w,L},   theta = r d/dr,
W_{w,L} = channel-L part (Legendre P_L(mu)) of rho_w prod_{c child of w} Phi_c,  Phi_c = sum_L Phi_{c,L} P_L(mu).
Root: P_0(d) is moved onto the children by 3-D integration by parts,
    int (P_0(d) K0(t|r|)) Psi d^3r = int K0(t|r|) P_0(-d) Psi d^3r,   Psi = prod_{w child of root} Phi_w,
so the fibre is a constant times t^3 T(G), T(g) = int_0^oo K0(tx) g(x) dx,  G = x^2 <P_0(-d) Psi>_{S^2}.

Model: an element is a dict  key = (k, i, kept, e) -> coefficient polynomial in z = r^2 (list, index m),
meaning  r^(e + 2(m - OFF)) K0^(k-i) D^i prod_{(y,L,st) in kept} (Phi_{y,L} if st = 0 else theta Phi_{y,L}),
D = theta K0, e in {0, 1}.  theta acts as a derivation:
    theta r^p = p r^p;  theta K0 = D;  theta D = r^2 K0;  theta Phi = (theta Phi);
    theta (theta Phi_{y,L}) = -(theta Phi) + L(L+1) Phi - r^2 W_{y,L}.
Coefficients: Fractions (mod=None) or ints mod a prime.
"""
from __future__ import annotations

from fractions import Fraction
from typing import Dict, List, Tuple

OFF = 6          # r-powers down to r^(-12) representable


# ------------------------------------------------------------------------------------------------ Legendre
def legendre(L: int) -> List[Fraction]:
    """Coefficients (ascending in mu) of P_L."""
    p0, p1 = [Fraction(1)], [Fraction(0), Fraction(1)]
    if L == 0:
        return p0
    for n in range(1, L):
        nxt = [Fraction(0)] * (n + 2)
        for i, c in enumerate(p1):
            nxt[i + 1] += Fraction(2 * n + 1, n + 1) * c
        for i, c in enumerate(p0):
            nxt[i] -= Fraction(n, n + 1) * c
        p0, p1 = p1, nxt
    return p1


def pmul(a: List[Fraction], b: List[Fraction]) -> List[Fraction]:
    out = [Fraction(0)] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        if x:
            for j, y in enumerate(b):
                out[i + j] += x * y
    return out


def mu_int(poly: List[Fraction]) -> Fraction:
    """(1/2) int_{-1}^{1} poly dmu (the S^2 average of an axial function)."""
    return sum((c / (k + 1) for k, c in enumerate(poly) if k % 2 == 0), Fraction(0))


def project(poly: List[Fraction]) -> Dict[int, Fraction]:
    """poly(mu) = sum_L c_L P_L(mu)."""
    out = {}
    for L in range(len(poly)):
        PL = legendre(L)
        c = Fraction(2 * L + 1) * mu_int(pmul(poly, PL))
        if c:
            out[L] = c
    return out


# ------------------------------------------------------------------------------------------------ one Bessel factor
# single radial factor: dict (i, q) -> Fraction, meaning r^q K0 (i = 0) or r^q D (i = 1)
def th1(f: Dict) -> Dict:
    out: Dict = {}
    for (i, q), c in f.items():
        out[(i, q)] = out.get((i, q), 0) + q * c
        if i == 0:
            out[(1, q)] = out.get((1, q), 0) + c
        else:
            out[(0, q + 2)] = out.get((0, q + 2), 0) + c
    return {k: v for k, v in out.items() if v}


def shift1(f: Dict, s: int) -> Dict:
    return {(i, q + s): c for (i, q), c in f.items()}


def lin1(*terms) -> Dict:
    out: Dict = {}
    for c, f in terms:
        for k, v in f.items():
            out[k] = out.get(k, 0) + c * v
    return {k: v for k, v in out.items() if v}


def dz_channels(ch: Dict[int, Dict], th, shift, lin) -> Dict[int, Dict]:
    """d_z on sum_L f_L(r) P_L(mu):  d_z [f P_L] = (L+1)/(2L+1) r^-1 (theta - L) f P_{L+1}
                                                + L/(2L+1) r^-1 (theta + L + 1) f P_{L-1}."""
    out: Dict[int, Dict] = {}
    for L, f in ch.items():
        tf = th(f)
        up = shift(lin((1, tf), (-L, f)), -1)
        out[L + 1] = lin((1, out.get(L + 1, {})), (Fraction(L + 1, 2 * L + 1), up))
        if L >= 1:
            dn = shift(lin((1, tf), (L + 1, f)), -1)
            out[L - 1] = lin((1, out.get(L - 1, {})), (Fraction(L, 2 * L + 1), dn))
    return {L: f for L, f in out.items() if f}


def lap_channels(ch: Dict[int, Dict], th, shift, lin) -> Dict[int, Dict]:
    """lap [f P_L] = r^-2 (theta(theta+1) - L(L+1)) f P_L."""
    out = {}
    for L, f in ch.items():
        tf = th(f)
        g = lin((1, th(tf)), (1, tf), (-L * (L + 1), f))
        g = shift(g, -2)
        if g:
            out[L] = g
    return out


def apply_P(Pv: Dict[Tuple[int, int], Fraction], ch: Dict[int, Dict], th, shift, lin, sign: int = 1) -> Dict:
    """sum c (sign d_z)^a lap^b applied to a channel dict."""
    out: Dict[int, Dict] = {}
    for (a, b), c in Pv.items():
        cur = ch
        for _ in range(b):
            cur = lap_channels(cur, th, shift, lin)
        for _ in range(a):
            cur = dz_channels(cur, th, shift, lin)
        for L, f in cur.items():
            out[L] = lin((1, out.get(L, {})), (Fraction(c) * sign ** a, f))
    return {L: f for L, f in out.items() if f}


def density_channels(Pv) -> Dict[int, Dict]:
    return apply_P(Pv, {0: {(0, 0): Fraction(1)}}, th1, shift1, lin1)


# ------------------------------------------------------------------------------------------------ the tree
class LTree:
    def __init__(self, n: int, edges, root: int, P: Dict[int, Dict]) -> None:
        self.n, self.root, self.P = n, root, P
        adj = {v: [] for v in range(n)}
        for a, b in edges:
            adj[a].append(b)
            adj[b].append(a)
        self.children: Dict[int, List[int]] = {}
        self.size: Dict[int, int] = {}
        order: List[int] = []

        def rec(v, p):
            self.children[v] = sorted(w for w in adj[v] if w != p)
            self.size[v] = 1 + sum(rec(w, v) for w in self.children[v])
            order.append(v)
            return self.size[v]
        assert rec(root, -1) == n
        self.post = order
        self.rho = {v: density_channels(P[v]) for v in range(n)}
        # W[w][L] = list of (coef, (i, q), ckept)  with ckept a tuple of (c, Lc, 0)
        self.W: Dict[int, Dict[int, List]] = {}
        self.chan: Dict[int, List[int]] = {}
        for w in self.post:
            if w == root:
                continue
            terms = []   # (angular poly, coef, (i,q), ckept)
            for L0, f in self.rho[w].items():
                for iq, c in f.items():
                    terms.append((legendre(L0), Fraction(c), iq, ()))
            for c_ in self.children[w]:
                new = []
                for poly, cf, iq, ck in terms:
                    for Lc in self.chan[c_]:
                        new.append((pmul(poly, legendre(Lc)), cf, iq, ck + ((c_, Lc, 0),)))
                terms = new
            WL: Dict[int, Dict] = {}
            for poly, cf, iq, ck in terms:
                for L, pc in project(poly).items():
                    key = (iq, tuple(sorted(ck)))
                    WL.setdefault(L, {})
                    WL[L][key] = WL[L].get(key, 0) + cf * pc
            self.W[w] = {L: [(v, k[0], k[1]) for k, v in d.items() if v] for L, d in WL.items()}
            self.W[w] = {L: l for L, l in self.W[w].items() if l}
            self.chan[w] = sorted(self.W[w])

    def psi_channels(self) -> Dict[int, List]:
        """Channel decomposition of Psi = prod_{w child of root} Phi_w: L -> list of (coef, kept)."""
        terms = [([Fraction(1)], ())]
        for w in self.children[self.root]:
            terms = [(pmul(poly, legendre(L)), ck + ((w, L, 0),)) for poly, ck in terms for L in self.chan[w]]
        out: Dict[int, Dict] = {}
        for poly, ck in terms:
            for L, pc in project(poly).items():
                key = tuple(sorted(ck))
                out.setdefault(L, {})
                out[L][key] = out[L].get(key, 0) + pc
        return {L: [(v, k) for k, v in d.items() if v] for L, d in out.items()}


# ------------------------------------------------------------------------------------------------ the module model
class Model:
    def __init__(self, tree: LTree, mod: int | None = None) -> None:
        self.t, self.mod = tree, mod

    def cv(self, c):
        if self.mod is None:
            return Fraction(c)
        c = Fraction(c)
        return (c.numerator % self.mod) * pow(c.denominator % self.mod, -1, self.mod) % self.mod

    def norm(self, x):
        return x % self.mod if self.mod else x

    def padd(self, d: Dict, key, poly, shift: int = 0, scale=1) -> None:
        if shift < 0:
            assert all(c == 0 for c in poly[:-shift]), "r-power below representable range"
            poly, shift = poly[-shift:], 0
        cur = d.get(key)
        if cur is None:
            cur = d[key] = []
        need = len(poly) + shift
        if len(cur) < need:
            cur.extend([0] * (need - len(cur)))
        for m, c in enumerate(poly):
            if c:
                cur[m + shift] = self.norm(cur[m + shift] + scale * c)

    def trim(self, d: Dict) -> Dict:
        out = {}
        for k, v in d.items():
            while v and v[-1] == 0:
                v.pop()
            if v:
                out[k] = v
        return out

    @staticmethod
    def rpow_key(e: int, m: int) -> int:
        return e + 2 * (m - OFF)

    @staticmethod
    def split_power(p: int) -> Tuple[int, int]:
        e = p % 2
        return e, (p - e) // 2 + OFF

    def mono(self, coef, k: int, i: int, p: int, kept: tuple) -> Dict:
        e, m = self.split_power(p)
        assert m >= 0
        poly = [0] * (m + 1)
        poly[m] = self.cv(coef)
        return {(k, i, tuple(sorted(kept)), e): poly}

    def theta(self, vec: Dict) -> Dict:
        out: Dict = {}
        W = self.t.W
        for (k, i, kept, e), c in vec.items():
            self.padd(out, (k, i, kept, e), [self.norm((e + 2 * (m - OFF)) * cm) for m, cm in enumerate(c)])
            if k - i > 0:
                self.padd(out, (k, i + 1, kept, e), c, 0, k - i)
            if i > 0:
                self.padd(out, (k, i - 1, kept, e), c, 1, i)
            for idx, (y, L, st) in enumerate(kept):
                rest = kept[:idx] + kept[idx + 1:]
                if st == 0:
                    self.padd(out, (k, i, tuple(sorted(rest + ((y, L, 1),))), e), c)
                else:
                    self.padd(out, (k, i, kept, e), c, 0, -1)
                    if L:
                        self.padd(out, (k, i, tuple(sorted(rest + ((y, L, 0),))), e), c, 0, L * (L + 1))
                    for cw, (j, q), ck in W[y][L]:
                        p_add = 2 + q
                        e2 = (e + p_add) % 2
                        sh = (e + p_add - e2) // 2
                        key = (k + 1, i + j, tuple(sorted(rest + ck)), e2)
                        self.padd(out, key, c, sh, self.norm(-self.cv(cw)))
        return self.trim(out)

    def shift(self, vec: Dict, s: int) -> Dict:
        """multiply by r^s."""
        out: Dict = {}
        for (k, i, kept, e), c in vec.items():
            e2 = (e + s) % 2
            sh = (e + s - e2) // 2
            self.padd(out, (k, i, kept, e2), c, sh)
        return self.trim(out)

    def lin(self, *terms) -> Dict:
        out: Dict = {}
        for c, v in terms:
            cc = self.cv(c)
            for k, poly in v.items():
                self.padd(out, k, poly, 0, cc)
        return self.trim(out)

    def G(self) -> Dict:
        """G = x^2 <P_0(-d) Psi>: channel 0 of P_0(-d) applied to the channel dict of Psi, times r^2."""
        tr = self.t
        ch = {}
        for L, lst in tr.psi_channels().items():
            v: Dict = {}
            for coef, kept in lst:
                v = self.lin((1, v), (1, self.mono(coef, 0, 0, 0, kept)))
            if v:
                ch[L] = v
        out = apply_P(tr.P[tr.root], ch, self.theta, self.shift, self.lin, sign=-1)
        return self.shift(out.get(0, {}), 2)

    def pure(self, key) -> bool:
        return key[0] == self.t.n - 1 and key[2] == ()


# ------------------------------------------------------------------------------------------------ densities
def Pdict(name: str) -> Dict[Tuple[int, int], Fraction]:
    """axial P(d) as {(a, b): c} = sum c d_z^a lap^b."""
    return {
        "s": {(0, 0): Fraction(1)},
        "pz": {(1, 0): Fraction(1)},
        "dz2": {(2, 0): Fraction(1)},                                     # d_z^2 (= (2/3) h_2 + lap/3)
        "h2": {(2, 0): Fraction(3, 2), (0, 1): Fraction(-1, 2)},          # h_2(d) = (3 d_z^2 - lap)/2
        "lap": {(0, 1): Fraction(1)},
        "cancel": {(0, 1): Fraction(9, 4), (2, 0): Fraction(-15, 4)},            # lap - (5/2) h_2
    }[name]


def symbol_poly(name: str) -> List[Fraction]:
    """P(-omega) restricted to the axis direction as a polynomial in mu (d_z -> -mu, lap -> 1)."""
    out = [Fraction(0)] * 3
    for (a, b), c in Pdict(name).items():
        out[a] += c * (-1) ** a
    return out
