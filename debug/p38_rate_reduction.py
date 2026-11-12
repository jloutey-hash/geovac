"""Paper 38 open item: O(1/n) rate. Reduction gamma(a) = (4/pi) int_0^pi theta |sum a_k sin k theta|^2 / sum a_k^2,
k = 2j+1. Matrix A = pi I - (4/pi) T. Calibrate on the paper's Cesaro numbers, then optimum + explicit tapers."""
import sys
import numpy as np
from scipy.linalg import eigh_tridiagonal, eigh

def A_matrix(ks):
    k = np.asarray(ks, float)[:, None]; l = np.asarray(ks, float)[None, :]
    d = k - l
    odd = (np.abs(np.rint(d)) % 2 == 1)
    with np.errstate(divide="ignore"):
        T = np.where(odd, 1.0 / d**2 - 1.0 / (k + l)**2, 0.0)
    return np.pi * np.eye(len(ks)) - (4 / np.pi) * T

def gamma(a, ks):
    a = np.asarray(a, float); A = A_matrix(ks)
    return a @ A @ a / (a @ a)

if __name__ == "__main__":
    # calibration: paper's Cesaro taper a_j = sqrt(2j+1)(1 - 2j/(n+1)), n*gamma = 6.35 (25), 6.60 (50)
    for n in (25, 50):
        for label, js in (("j half-int", np.arange(0, (n + 1) / 2 + 1e-9, 0.5)), ("j int", np.arange(0, (n + 1) // 2 + 1))):
            js = js[1 - 2 * js / (n + 1) >= 0]
            a = np.sqrt(2 * js + 1) * (1 - 2 * js / (n + 1))
            print(n, label, "N=", len(js), "n*gamma=%.4f" % (n * gamma(a, 2 * js + 1)))

def run():
    print("calibration (k = 1..n):")
    for n, ref in ((25, 6.35), (50, 6.60), (100, 6.75)):
        k = np.arange(1, n + 1); j = (k - 1) / 2
        a = np.sqrt(2 * j + 1) * (1 - 2 * j / (n + 1))
        print("  n=%d n*gamma=%.4f (paper %.2f)" % (n, n * gamma(a, k), ref))
    print("optimum n*lambda_min(A) and the sine taper a_k = sin(pi k/(n+1)):")
    for n in (25, 50, 100, 200, 400, 800, 1600, 3200):
        k = np.arange(1, n + 1); A = A_matrix(k)
        lam = eigh(A, eigvals_only=True, subset_by_index=[0, 0])[0]
        s = np.sin(np.pi * k / (n + 1))
        print("  n=%5d opt %.6f  sine %.6f" % (n, n * lam, n * gamma(s, k)), flush=True)
run()
