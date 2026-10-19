"""Chain energy with pairs at most m steps apart: e_m(g,mu) = mu*gbar + (2/p) sum_i sum_{r<=m} w(y_{i+r}-y_i).
Minimised over periodic patterns (period <= P) with gaps >= 11/20; prints (H-mu)/(1-e_m), the limit for a window of m+1 points.
usage: python3 -I chain_window.py P mu1,mu2,... [m1,m2,...]   (default m = 6,7,8; chain_window8.txt: P = 4;
chain_window9.txt: P = 4, m = 8,9,10,11)"""
import itertools, sys, numpy as np
from scipy.optimize import minimize
s2 = np.sqrt(2.0); H = 1.5 - (np.cos(1/s2)/np.sin(1/s2))/s2; K0 = s2*np.sin(1/s2)
def k(x): return 0.5*(np.sinc((np.pi*x - 1/s2)/np.pi) + np.sinc((np.pi*x + 1/s2)/np.pi))/K0
def Em(g, m):
    p = len(g); y = np.concatenate([[0.0], np.cumsum(np.tile(g, m//p + 2))]); e = 0.0
    for i in range(p):
        e += 2*np.sum(k(y[i+1:i+1+m] - y[i])**2)/p
    return e
P = int(sys.argv[1]) if len(sys.argv) > 1 else 4
for mu in [float(x) for x in sys.argv[2].split(",")]:
    row = []
    for m in ([int(v) for v in sys.argv[3].split(",")] if len(sys.argv) > 3 else [6, 7, 8]):
        b = np.inf
        for p in range(1, P+1):
            for word in itertools.product([1.04, 1.98, 3.0], repeat=p):
                r = minimize(lambda g: mu*np.mean(g) + Em(g, m), np.array(word), method="L-BFGS-B", bounds=[(0.55, 60)]*p)
                b = min(b, r.fun)
        row.append(f"m={m}:{(H-mu)/(1-b):.7f}")
    print(f"mu={mu}", " ".join(row), flush=True)
