"""T2 evaluator with per-corner treatment.  Diagnosis (see sprint_T2_outer_corner_memo.md):

  corner   b=s+t   J behaviour              treatment              fibre
  (0,0)    ->0     rho^{3/2}+rho^{5/2}, sing Duffy rho=sigma^2      Jdecay (non-osc, spectral)
  (1,1)    ->2     analytic (rho^2,rho^3)    plain GL (moderate N)  quadosc (osc-robust)
  (0,1)    ->1     analytic                 plain GL               quadosc
  (1,0)    ->1     = (0,1) by s<->t          -- (folded)            --
  bulk / edges     analytic                 plain GL               Jdecay (edge k-cutoff)

T2 = (8/pi) [ I_00 + I_11 + 2*I_01 + I_bulk ].
The "outer wall" the corpus reported was FIBRE error at the 3 oscillatory corners
(the memo's tail/tensor fibres gave ~12-19 dig there), NOT outer-quadrature convergence.
"""
import sys, time
import mpmath as mp
sys.path.insert(0, 'debug')
from _fastgl import fast_gl
from _t2_decayfibre import Jdecay
from _t2_fibre_stress import J_quadosc

# ---- fibres -------------------------------------------------------------
OSC_THRESH = 3.5   # osc_real above this -> quadosc; below -> decay-map GL

def Jsmart(s, t):
    """Uniform fibre. osc_real = b/(sqrt cs+sqrt ct) = # of j0 periods across the
    envelope (envelope decays like e^{-(sqrt cs+sqrt ct)k}). Few periods -> decay-map
    GL is spectral; many -> route to quadosc (oscillation-robust)."""
    cs = s*(1-s); ct = t*(1-t)
    if cs == 0 or ct == 0:
        return mp.mpf(0)
    b = s+t
    osc = float(b/(mp.sqrt(cs)+mp.sqrt(ct)))
    if osc <= OSC_THRESH:
        Nk = min(max(int(600 + 250*osc), 400), 1600)   # calibrated for >=42 dig
        Nk = int(round(Nk/100.0))*100                  # quantize -> fast_gl cache reuse
        return Jdecay(s, t, Nk)
    return J_quadosc(s, t)

# back-compat aliases
Jfast = Jsmart
Josc = J_quadosc

# ---- (0,0) singular corner: Duffy rho=sigma^2 + angular sin^2 -----------
def I_00(delta, Nsig, Nal, fibre):
    H = mp.pi/2
    xs, ws = fast_gl(Nsig); xa, wa = fast_gl(Nal)
    smax = mp.sqrt(delta); tot = mp.mpf(0)
    for xi, wi in zip(xs, ws):
        sig = smax*(xi+1)/2; wsig = smax*wi/2
        for xj, wj in zip(xa, wa):
            psi = H*(xj+1)/2; wpsi = H*wj/2
            al = mp.sin(psi)**2; jal = mp.sin(2*psi)
            s = sig*sig*al; t = sig*sig*(1-al)
            tot += wsig*wpsi*jal*2*sig**3*fibre(s, t)
    return tot

# ---- analytic corner patch (plain GL on the delta-square) ---------------
def I_patch(corner, delta, N, fibre):
    xs, ws = fast_gl(N)
    (cx, cy) = corner  # each 0 or 1
    tot = mp.mpf(0)
    for xi, wi in zip(xs, ws):
        s = (delta*(xi+1)/2) if cx == 0 else (1-delta+delta*(xi+1)/2)
        wsi = delta*wi/2
        for xj, wj in zip(xs, ws):
            t = (delta*(xj+1)/2) if cy == 0 else (1-delta+delta*(xj+1)/2)
            wtj = delta*wj/2
            tot += wsi*wtj*fibre(s, t)
    return tot

# ---- bulk: R1=[d,1-d]x[0,1], R2=[0,d]x[d,1-d], R3=[1-d,1]x[d,1-d] --------
def _rect(ax, bx, ay, by, Nx, Ny, fibre):
    xs, wx = fast_gl(Nx); ys, wy = fast_gl(Ny)
    tot = mp.mpf(0)
    for xi, wi in zip(xs, wx):
        s = ax+(bx-ax)*(xi+1)/2; wsi = (bx-ax)*wi/2
        for xj, wj in zip(ys, wy):
            t = ay+(by-ay)*(xj+1)/2; wtj = (by-ay)*wj/2
            tot += wsi*wtj*fibre(s, t)
    return tot

def I_bulk(delta, N, fibre):
    d = delta; one = mp.mpf(1)
    r1 = _rect(d, 1-d, mp.mpf(0), one, N, N, fibre)          # middle vertical strip
    r2 = _rect(mp.mpf(0), d, d, 1-d, N, N, fibre)            # left middle
    r3 = _rect(1-d, one, d, 1-d, N, N, fibre)                # right middle
    return r1+r2+r3

def I_sliver(delta, N, fibre):
    """Analytic sliver triangle {s in[0,d], t in[d-s,d]} (the part of the square [0,d]^2
    ABOVE the Duffy triangle s+t<=d). b=s+t in [d,2d] small -> non-oscillatory -> fast.
    Map s=d*u, t=(d-s)+s*v; area element = d^2 * u du dv."""
    d = delta; xs, ws = fast_gl(N)
    tot = mp.mpf(0)
    for xi, wi in zip(xs, ws):
        u = (xi+1)/2; wu = wi/2; s = d*u
        for xj, wj in zip(xs, ws):
            v = (xj+1)/2; wv = wj/2
            t = (d - s) + s*v
            tot += wu*wv*d*d*u*fibre(s, t)
    return tot

def T2_Lshape(delta, N00, Nrect, fibre=None, verbose=True):
    """CLEAN 3-piece decomposition (provably-correct, matches direct GL):
       T2 = (8/pi)[ Duffy([0,d]^2) + GL([d,1]x[0,1]) + GL([0,d]x[d,1]) ].
    Only (0,0) (singular rho^{3/2}) gets Duffy; the 3 analytic oscillatory corners live
    inside the two big rectangles (plain GL, Jsmart auto-routes quadosc at their deep nodes)."""
    if fibre is None:
        fibre = Jsmart
    d = delta; one = mp.mpf(1); comp = {}
    t0 = time.time(); comp['duffy'] = I_00(d, N00, N00, fibre)      # triangle {s+t<=d}
    if verbose: print(f"    duffy tri{{s+t<=d}} = {mp.nstr(comp['duffy'],30)}  ({time.time()-t0:.0f}s)", flush=True)
    t0 = time.time(); comp['sliver'] = I_sliver(d, max(12, N00//2), fibre)  # square-minus-triangle
    if verbose: print(f"    sliver           = {mp.nstr(comp['sliver'],30)}  ({time.time()-t0:.0f}s)", flush=True)
    t0 = time.time(); comp['rectR'] = _rect(d, one, mp.mpf(0), one, Nrect, Nrect, fibre)
    if verbose: print(f"    rect[d,1]x[0,1]  = {mp.nstr(comp['rectR'],30)}  ({time.time()-t0:.0f}s)", flush=True)
    t0 = time.time(); comp['rectT'] = _rect(mp.mpf(0), d, d, one, Nrect, Nrect, fibre)
    if verbose: print(f"    rect[0,d]x[d,1]  = {mp.nstr(comp['rectT'],30)}  ({time.time()-t0:.0f}s)", flush=True)
    v = (8/mp.pi)*(comp['duffy'] + comp['sliver'] + comp['rectR'] + comp['rectT'])
    return v, comp

def T2(delta=mp.mpf('0.1'), N00=(40, 40), Nc=20, Nbulk=48, verbose=True):
    t0 = time.time()
    i00 = I_00(delta, N00[0], N00[1], Jfast)
    if verbose: print(f"    I_00   = {mp.nstr(i00,28)}  ({time.time()-t0:.0f}s)", flush=True)
    t0 = time.time()
    i11 = I_patch((1, 1), delta, Nc, Josc)
    if verbose: print(f"    I_11   = {mp.nstr(i11,28)}  ({time.time()-t0:.0f}s)", flush=True)
    t0 = time.time()
    i01 = I_patch((0, 1), delta, Nc, Josc)
    if verbose: print(f"    I_01   = {mp.nstr(i01,28)}  ({time.time()-t0:.0f}s)", flush=True)
    t0 = time.time()
    ibulk = I_bulk(delta, Nbulk, Jfast)
    if verbose: print(f"    I_bulk = {mp.nstr(ibulk,28)}  ({time.time()-t0:.0f}s)", flush=True)
    return (8/mp.pi)*(i00 + i11 + 2*i01 + ibulk)

def full_T2(delta, N00, Nc, Nbulk, verbose=True):
    """Return (T2, components dict)."""
    comp = {}
    t0 = time.time(); comp['I_00'] = I_00(delta, N00, N00, Jsmart)
    if verbose: print(f"    I_00   = {mp.nstr(comp['I_00'],30)}  ({time.time()-t0:.0f}s)", flush=True)
    t0 = time.time(); comp['I_11'] = I_patch((1, 1), delta, Nc, Jsmart)
    if verbose: print(f"    I_11   = {mp.nstr(comp['I_11'],30)}  ({time.time()-t0:.0f}s)", flush=True)
    t0 = time.time(); comp['I_01'] = I_patch((0, 1), delta, Nc, Jsmart)
    if verbose: print(f"    I_01   = {mp.nstr(comp['I_01'],30)}  ({time.time()-t0:.0f}s)", flush=True)
    t0 = time.time(); comp['I_bulk'] = I_bulk(delta, Nbulk, Jsmart)
    if verbose: print(f"    I_bulk = {mp.nstr(comp['I_bulk'],30)}  ({time.time()-t0:.0f}s)", flush=True)
    v = (8/mp.pi)*(comp['I_00'] + comp['I_11'] + 2*comp['I_01'] + comp['I_bulk'])
    return v, comp

if __name__ == '__main__':
    import json
    mp.mp.dps = int(sys.argv[1]) if len(sys.argv) > 1 else 30
    delta = mp.mpf(sys.argv[2]) if len(sys.argv) > 2 else mp.mpf('0.1')
    Nc = int(sys.argv[3]) if len(sys.argv) > 3 else 12
    Nbulk = int(sys.argv[4]) if len(sys.argv) > 4 else 40
    N00 = int(sys.argv[5]) if len(sys.argv) > 5 else 32
    ref = mp.mpf('0.3953557659017139641')
    print(f"RUN dps={mp.mp.dps} delta={delta} Nc={Nc} Nbulk={Nbulk} N00={N00}", flush=True)
    t0 = time.time()
    v, comp = full_T2(delta, N00, Nc, Nbulk)
    print(f"  T2 = {mp.nstr(v, mp.mp.dps-4)}   ({time.time()-t0:.0f}s total)", flush=True)
    print(f"  |T2 - anchor19| = {mp.nstr(abs(v-ref),4)}", flush=True)
    out = {'dps': mp.mp.dps, 'delta': str(delta), 'Nc': Nc, 'Nbulk': Nbulk, 'N00': N00,
           'T2': mp.nstr(v, mp.mp.dps), 'components': {k: mp.nstr(val, mp.mp.dps) for k, val in comp.items()},
           'abs_vs_anchor19': mp.nstr(abs(v-ref), 6)}
    with open(f'debug/data/t2_corner_run_dps{mp.mp.dps}_Nc{Nc}_Nb{Nbulk}.json', 'w') as f:
        json.dump(out, f, indent=1)
    print("  wrote debug/data/...", flush=True)
