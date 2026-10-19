"""Largest Gram eigenvalue for blocks of n zeros with every gap >= delta (multistart over gaps)."""
import sys, numpy as np
sys.path.insert(0, sys.argv[1]); from blocks import k
from scipy.optimize import minimize
n = int(sys.argv[2]); rng = np.random.default_rng(0)
def lam(g):
    y = np.concatenate([[0], np.cumsum(g)]); return np.linalg.eigvalsh(k(y[:, None]-y[None, :]))[-1]
for delta in [0.4, 0.5, 0.55, 0.6, 0.7, 0.8]:
    best = 0
    for _ in range(60):
        g0 = delta + rng.exponential(0.3, n-1)
        r = minimize(lambda g: -lam(g), g0, method="L-BFGS-B", bounds=[(delta, 6)]*(n-1))
        best = max(best, -r.fun)
    print(f"n={n} delta={delta} max lambda={best:.4f}  (k(delta)={k(delta):.3f})", flush=True)
