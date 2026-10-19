"""Upper bounds for the chain ground state e(mu) = inf over infinite gap sequences (gaps >= 0.55) of the
per-zero value of mu*gap + 2*sum_{j>i} w(y_j - y_i), from periodic gap patterns (period p <= 6).
Gives an upper limit for the uncapped method: max_mu (H - mu)/(1 - e(mu))."""
import itertools, numpy as np
from scipy.optimize import minimize
s2 = np.sqrt(2.0); H = 1.5 - (np.cos(1/s2)/np.sin(1/s2))/s2; K0 = s2*np.sin(1/s2)
def k(x): return 0.5*(np.sinc((np.pi*x - 1/s2)/np.pi) + np.sinc((np.pi*x + 1/s2)/np.pi))/K0
def energy(g, mu, reach=400):
    p = len(g); y = np.concatenate([[0.0], np.cumsum(np.tile(g, reach//p + 2))])
    e = mu*np.mean(g)
    for i in range(p):
        d = y[i+1:i+1+reach] - y[i]; e += 2*np.sum(k(d)**2)/p
    return e
for mu in [0.002, 0.0025, 0.003, 0.0035, 0.004, 0.005, 0.006]:
    best = (np.inf, None)
    for p in range(1, 7):
        for word in itertools.product([1.04, 1.98, 3.0], repeat=p):
            r = minimize(energy, np.array(word), args=(mu,), method="L-BFGS-B", bounds=[(0.55, 20)]*p)
            if r.fun < best[0]: best = (r.fun, np.round(r.x, 4))
    e = best[0]
    print(f"mu={mu}: e_upper={e:.7f} limit<= {(H-mu)/(1-e):.7f} pattern={best[1]}", flush=True)
