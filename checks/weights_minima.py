"""Local minima of the weighted window function for a rational weights file.
usage: python3 -I weights_minima.py R weights.txt   (weights_minima8.txt: R = 8, verify/weights8.txt)"""
import itertools, sys, numpy as np
from scipy.optimize import minimize
s2 = np.sqrt(2.0); K0 = s2*np.sin(1/s2); t0 = 1/(s2*np.pi)
def dsinc(t):
    t = np.asarray(t, float); out = np.zeros_like(t); nz = np.abs(t) > 1e-12
    out[nz] = (np.cos(np.pi*t[nz]) - np.sinc(t[nz]))/t[nz]; return out
def k(x): return 0.5*(np.sinc(x - t0) + np.sinc(x + t0))/K0
def dk(x): return 0.5*(dsinc(x - t0) + dsinc(x + t0))/K0
R = int(sys.argv[1]); n = R - 1
pairs = [(r, i) for r in range(1, R) for i in range(R - r)]
x = np.zeros(n + len(pairs))
for line in open(sys.argv[2]):
    f = line.split()
    if f[0] == "b": x[int(f[1]) - 1] = int(f[2])/int(f[3])
    else: x[n + pairs.index((int(f[1]), int(f[2])))] = int(f[3])/int(f[4])
M = np.zeros((len(pairs), n))
for j, (r, i) in enumerate(pairs): M[j, i:i+r] = 1.0
def fg(g):
    D = M @ g; kd = k(D); return x[:n] @ g + x[n:] @ kd**2, x[:n] + M.T @ (x[n:]*2*kd*dk(D))
starts = [np.array(wd) for wd in itertools.product([1.04, 1.98, 3.0], repeat=n)]
rng = np.random.default_rng(2); starts += [rng.uniform(0.6, 4, n) for _ in range(3000)]
best = {}
for s in starts:
    r = minimize(fg, s, jac=True, method="L-BFGS-B", bounds=[(0.3, 60)]*n, options={"ftol": 1e-15, "gtol": 1e-11})
    key = tuple(np.round(r.x, 3)); best[key] = min(best.get(key, 1e9), r.fun)
out = sorted((f, g) for g, f in best.items())
for f, g in out[:25]: print(f"{f:.7f}", g)
