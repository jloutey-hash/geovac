import sys, time
sys.path.insert(0, r'C:\Users\jlout\Desktop\Project_Geometric\debug')
import mpmath as mp
import aadx_ap0 as A
import aapc_exact as P

for case in P.CASES:
    for integral, L in (('DX', 12), ('DD', 12), ('DX', 20), ('DD', 20)):
        t = time.time()
        q32 = A.quad_GL(case, L, integral, dps=32)
        q44 = A.quad_GL(case, L, integral, dps=44)
        with mp.workdps(60):
            print(case, integral, L, 'q44', mp.nstr(q44, 12), 'rel |q32/q44 - 1|', mp.nstr(abs(q32 / q44 - 1), 3),
                  round(time.time() - t), 's', flush=True)
