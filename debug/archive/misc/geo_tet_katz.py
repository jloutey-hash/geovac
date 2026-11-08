"""Blind geometric route, part 2: Katz middle convolution of local monodromy tuples.

Y/t^3 = const * int_0^inf sqrt(z) (z + x)^(-3/2) Phi(z) dz  (z = L^2, Phi the kite), an
Euler transform with kernel of monodromy -1.  So, up to rank-one pieces created at the
boundary, each irreducible local system W of g(v) = sqrt(z) Phi(z) (v = -z) goes to Katz's
middle convolution MC_{-1}(W), which is irreducible (W irreducible, not a Kummer sheaf).

This script implements the Dettweiler-Reiter algorithm for MC_lambda on tuples
(A_1..A_r) of monodromy matrices (exact rationals, sympy), validates it on the Legendre
case and on MC_{-1} o MC_{-1} = id, and applies it to geometric inputs only:

  W_ab : the D = 3 density bubble (b = 2, a = 1) times sqrt(z).  The bubble is
         P'_{-1/2}(z_L), z_L = (a^2+b^2-v)/(2ab) (Paper 59, function L(a,b;c) at c = L),
         so its local system is Legendre's (unipotent at v = (a-b)^2, (a+b)^2, -J2 at
         infinity, trivial at v = 0); sqrt(z) adds -1 at v = 0.
  slice: the same at a = b = 1 (control against Paper 59's order-3 K3 piece).
  E_A  : rank one, -1 at v = 0 (sqrt z), 4, 16 (the density-geminal three-line cut).
  E_G  : the kite top-sector cut curve (geminal pair, I_1 at 9/4 and 81/4, I_4^* at oo),
         with and without a sqrt(v) twist at 0: shows what decides the top sector's image.

No differential operator for Y or for any sector is computed.
"""
import json
import sympy as sp

I2 = sp.eye(2)


def mc(As, lam):
    """Dettweiler-Reiter middle convolution MC_lam of the tuple As (list of sympy matrices)."""
    r = len(As)
    n = As[0].shape[0]
    N = n * r
    Bs = []
    for k in range(r):
        B = sp.eye(N)
        for j in range(r):
            if j < k:
                blk = lam * (As[j] - sp.eye(n))
            elif j == k:
                blk = lam * As[k]
            else:
                blk = As[j] - sp.eye(n)
            B[k * n:(k + 1) * n, j * n:(j + 1) * n] = blk
        Bs.append(B)
    # K = sum of ker(A_i - 1) placed in block i; L = intersection of ker(B_k - 1)
    Kvecs = []
    for i, A in enumerate(As):
        for v in (A - sp.eye(n)).nullspace():
            w = sp.zeros(N, 1)
            w[i * n:(i + 1) * n, 0] = v
            Kvecs.append(w)
    stack = sp.Matrix.vstack(*[B - sp.eye(N) for B in Bs])
    Lvecs = stack.nullspace()
    sub = Kvecs + Lvecs
    if sub:
        S = sp.Matrix.hstack(*sub)
        # basis of the span
        S = sp.Matrix.hstack(*[S[:, i] for i in S.T.rref()[1]]) if S.rank() < S.shape[1] else S
        rS = S.rank()
    else:
        S = sp.zeros(N, 0)
        rS = 0
    # complement: extend S to a basis of Q^N
    cols = [S[:, i] for i in range(S.shape[1])]
    basis = list(cols)
    for i in range(N):
        e = sp.zeros(N, 1)
        e[i] = 1
        if sp.Matrix.hstack(*(basis + [e])).rank() > len(basis):
            basis.append(e)
    Bm = sp.Matrix.hstack(*basis)
    Bi = Bm.inv()
    m = N - rS
    out = []
    for B in Bs:
        C = Bi * B * Bm
        # invariant subspace = first rS coordinates; quotient block = lower-right
        assert all(C[rS + i, j] == 0 for i in range(m) for j in range(rS)), 'subspace not invariant'
        out.append(sp.Matrix(C[rS:, rS:]))
    return out


def jordan(A):
    """Jordan block structure {eigenvalue: [sizes]} of A."""
    P, J = A.jordan_form()
    blocks = {}
    i = 0
    n = J.shape[0]
    while i < n:
        lam = J[i, i]
        s = 1
        while i + s < n and J[i + s - 1, i + s] == 1 and J[i + s, i + s] == lam:
            s += 1
        blocks.setdefault(str(sp.nsimplify(lam)), []).append(s)
        i += s
    return blocks


def describe(As, names):
    prod = sp.eye(As[0].shape[0])
    for A in As:
        prod = prod * A
    Tinf = prod.inv()
    res = {}
    for nm, A in list(zip(names, As)) + [('oo', Tinf)]:
        res[nm] = {'jordan': jordan(A), 'rank(T-1)': (A - sp.eye(A.shape[0])).rank()}
    return res


def invariant_forms(As):
    """Solve A^T F A = F for all A; return (dim, dim symmetric, dim antisymmetric)."""
    n = As[0].shape[0]
    syms = sp.symbols('f0:%d' % (n * n))
    F = sp.Matrix(n, n, syms)
    eqs = []
    for A in As:
        eqs += list(A.T * F * A - F)
    sol = sp.linsolve(eqs, syms)
    sol = list(sol)[0]
    free = sorted(set().union(*[s.free_symbols for s in sol]), key=str)
    Fg = sp.Matrix(n, n, sol)
    sym_dim = len(free)
    # split into symmetric / antisymmetric parts
    basis = []
    for f in free:
        Fi = Fg.subs({g: (1 if g == f else 0) for g in free})
        basis.append(Fi)
    nsym = sum(1 for Fi in basis if Fi == Fi.T)
    nanti = sum(1 for Fi in basis if Fi == -Fi.T)
    return {'dim': len(free), 'symmetric': nsym, 'antisymmetric': nanti,
            'forms': [str(Fi.tolist()) for Fi in basis]}


out = {}
M0 = sp.Matrix([[1, 2], [0, 1]])
M1 = sp.Matrix([[1, 0], [-2, 1]])

# ------------------------------------------------ validation 1: Legendre from x^(-1/2)(x-1)^(-1/2)
leg = mc([sp.Matrix([[-1]]), sp.Matrix([[-1]])], -1)
out['val_legendre'] = {'rank': leg[0].shape[0], 'local': describe(leg, ['0', '1'])}
# ------------------------------------------------ validation 2: MC_{-1} o MC_{-1} = id
back = mc(leg, -1)
out['val_involution'] = {'rank': back[0].shape[0], 'local': describe(back, ['0', '1'])}

# ------------------------------------------------ W_ab off the slice: points v = 0, 1, 9
W = [-I2, M0, M1]
out['W_ab_input'] = describe(W, ['v=0', 'v=1', 'v=9'])
K3 = mc(W, -1)
out['K3_offslice'] = {'rank': K3[0].shape[0], 'local': describe(K3, ['x=0', 'x=1', 'x=9']),
                      'invariant_forms': invariant_forms(K3)}
# ------------------------------------------------ slice control: a = b = 1, points v = 0, 4
Ws = [-M0, M1]
K3s = mc(Ws, -1)
out['K3_slice'] = {'rank': K3s[0].shape[0], 'local': describe(K3s, ['x=0', 'x=4']),
                   'invariant_forms': invariant_forms(K3s)}
# ------------------------------------------------ E_A from the rank-one cut character
EA = mc([sp.Matrix([[-1]])] * 3, -1)
out['E_A'] = {'rank': EA[0].shape[0], 'local': describe(EA, ['x=0', 'x=4', 'x=16']),
              'invariant_forms': invariant_forms(EA)}
# ------------------------------------------------ E_G (geminal-pair cut curve), two twists
Da = sp.Matrix([[1, 1], [0, 1]])
# second vanishing cycle chosen so that (Da*Db) has trace -2 (I_4^* at infinity)
Db = sp.Matrix([[1, 0], [-4, 1]])
assert (Da * Db).trace() == -2
for tw, A0 in [('untwisted_at_0', I2), ('sqrt_twisted_at_0', -I2)]:
    try:
        EG = mc([A0, Da, Db], -1)
        out['E_G_' + tw] = {'rank': EG[0].shape[0],
                            'local': describe(EG, ['x=0', 'x=9/4', 'x=81/4']) if EG[0].shape[0] else {}}
    except Exception as exc:  # noqa
        out['E_G_' + tw] = 'error: ' + str(exc)

with open('debug/data/geo_tet_katz.json', 'w') as fh:
    json.dump(out, fh, indent=1)
for k, v in out.items():
    print(k, ':', json.dumps(v))
