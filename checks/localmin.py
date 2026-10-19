"""Local minima of F_R(g) = mu0*sum g + sum_{r<R} (2/(R-r)) sum_i w(y_{i+r}-y_i) over gaps g_i >= 0.55,
started from every word in {1.04, 1.98, 3.0}^(R-1) (the minimisers are subsets of a near-integer lattice).
usage: python3 -I localmin.py R mu0_den[,mu0_den...]"""
import sys, itertools, numpy as np
from scipy.optimize import minimize
s2 = np.sqrt(2.0); H = 1.5 - (np.cos(1/s2)/np.sin(1/s2))/s2; K0 = s2*np.sin(1/s2)
def k(x): return 0.5*(np.sinc((np.pi*x - 1/s2)/np.pi) + np.sinc((np.pi*x + 1/s2)/np.pi))/K0
def F(g, R, mu0):
    y = np.concatenate([[0.0], np.cumsum(g)]); v = mu0*np.sum(g)
    for r in range(1, R): v += 2.0/(R-r)*np.sum(k(y[r:]-y[:-r])**2)
    return v
R = int(sys.argv[1])
for den in [int(d) for d in sys.argv[2].split(',')]:
    mu0 = 1.0/den; best = (np.inf, None)
    for word in itertools.product([1.04, 1.98, 3.0], repeat=R-1):
        r = minimize(F, np.array(word), args=(R, mu0), method="L-BFGS-B", bounds=[(0.55, 40)]*(R-1))
        if r.fun < best[0]: best = (r.fun, r.x)
    c0 = best[0]; mu = (R-1)*mu0
    print(f"R={R} mu0=1/{den} c0~{c0:.7f} bound~{(H-mu)/(1-c0):.7f} gaps={np.round(best[1],4)}", flush=True)
