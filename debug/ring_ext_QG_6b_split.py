"""QG step 6b (gamma = 1): does the K3 piece of G4 = R3 o R1 split off (f3 a right factor of G4)?
Bounded mod-p search on side R.  Positive control first: the reversed composition R1 o R3 has f3 as a
right factor; the search bounds are fixed where the control succeeds, then G4 is tested with those and
larger bounds."""
import sys, json; sys.path.insert(0,'debug')
import sympy as sp
import ring_ext_QG_ops as O, ring_ext_QG_6_factor as F, density_loop_D_gauge as D
from ring_ext_QG_ops import tt as t
from fractions import Fraction
k3=json.load(open('debug/data/ring_ext_QG_7_k3map_g1.json')); d=json.load(open('debug/data/ring_ext_QG_6_gcrd_g1.json'))
G4=[sp.sympify(x,locals={'t':t}) for x in d['G4']]; R3=[sp.sympify(x,locals={'t':t}) for x in k3['R3_monic']]
r=sp.sympify(k3['R1_r'],locals={'t':t}); R1=[-r,sp.Integer(1)]
f3 = D.hyper_t([Fraction(1, 2)] * 3, [1, 1], NT=F.NT, scale=Fraction(1, 4))
ap=9*t**4+312*t**2-12560; q6=33*t**6-612*t**4+5232*t**2-49088
QR3=O.coeffs_to_Q(R3); Qrev=O.coeffs_to_Q(O.comp(R1,R3)); QG4=O.coeffs_to_Q(G4)
found=None
for den in [sp.Integer(1),(t**2-4),(t**2-4)**2,ap,(t**2-4)*ap,(t**2-4)**2*ap,q6,(t**2-4)*q6,ap*q6,(t**2-4)*ap*q6,(t**2-4)**2*ap*q6, ap**2*q6, (t**2-4)**2*ap**2*q6]:
    for Dg in (10,20,30,40):
        n=F.test(QR3,f3,3,Dg,den)
        if n: found=(den,Dg); print('R3 <- L3 map: den',den,'Dg',Dg,'nullity',n, flush=True); break
    if found: break
den,Dg=found
print('control R1 o R3 with same bounds:', F.test(Qrev,f3,3,Dg,den), flush=True)
for den2 in [den, den*(t**2-4), den*(t**2-4)**2*ap]:
    for Dg2 in (Dg, Dg+10, Dg+20):
        print('G4 (R) den',sp.factor(den2),'Dg',Dg2,':', F.test(QG4,f3,3,Dg2,den2), flush=True)
