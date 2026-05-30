import numpy as np
from scipy.special import eval_gegenbauer, factorial
trap = np.trapezoid  # numpy>=2.0 (np.trapz removed)

# EXACT target: hydrogen 1s position entropy S_x = 3 + ln(pi) = 4.144729...
r=np.linspace(1e-9,60,2000000)
R=2*np.exp(-r)
rho=R**2
print("norm_x =", trap(rho*r**2, r), " (target 1)")
radial = trap(np.where(rho>1e-300,-rho*np.log(rho)*r**2,0.0), r)
Sx = radial + np.log(4*np.pi)
print("S_x(1s) =", round(Sx,5), "  EXACT 3+ln(pi) =", round(3+np.log(np.pi),5))

eta=np.linspace(1e-9,300,3000000)
pref=np.sqrt(2/np.pi*factorial(0)/factorial(1))*1.0*4*1
G=pref*eval_gegenbauer(0,1,(eta**2-1)/(eta**2+1))/(eta**2+1)**2
print("norm_p =", trap(G**2*eta**2, eta), " (target 1)")
radp = trap(np.where(G**2>1e-300,-G**2*np.log(G**2)*eta**2,0.0), eta)
Sp = radp + np.log(4*np.pi)
print("S_p(1s) =", round(Sp,5))
print("SUM(1s) =", round(Sx+Sp,5), "  BBM floor 3(1+ln pi) =", round(3*(1+np.log(np.pi)),5))
