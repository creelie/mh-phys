"""Uncapped block defect: D(M) >= (1/m) sum over all windows of m consecutive simple zeros of trPsi(G_window).
Find c(mu) = min_g [ trPsi(G(g))/m + mu*span(g)/(m-1) ] and the bound (H - mu)/(1 - c)."""
import sys, numpy as np
from scipy.optimize import minimize
s2 = np.sqrt(2.0); H = 1.5 - np.cot(1/s2)/s2 if hasattr(np, 'cot') else 1.5 - (np.cos(1/s2)/np.sin(1/s2))/s2
def K(x):  # closed form of int_{-1/2}^{1/2} cos(sqrt2 t) cos(2 pi x t) dt
    a = np.pi*x - 1/s2; b = np.pi*x + 1/s2
    return 0.5*(np.sinc(a/np.pi) + np.sinc(b/np.pi))
K0 = K(0.0)
def k(x): return K(x)/K0
def Psi(l): return np.where(l <= 2, (l-1)**2, 2*l-3)
def F(g):
    y = np.concatenate([[0], np.cumsum(np.abs(g))])
    G = k(y[:, None]-y[None, :])
    return Psi(np.linalg.eigvalsh(G)).sum()
def cmin(m, mu, starts=300, rng=np.random.default_rng(0)):
    best = np.inf; arg = None
    obj = lambda g: F(g)/m + mu*np.abs(g).sum()/(m-1)
    for _ in range(starts):
        g0 = rng.uniform(0, 3, m-1) * rng.choice([0.3, 1, 1.5])
        r = minimize(obj, g0, method="Nelder-Mead", options=dict(xatol=1e-7, fatol=1e-10, maxiter=4000))
        if r.fun < best: best, arg = r.fun, np.abs(r.x)
    return best, arg
if __name__ == "__main__":
    m = int(sys.argv[1])
    print("H =", H, " k(1.0573) =", k(1.0572782919432706))
    out = []
    for mu in [float(v) for v in sys.argv[2].split(',')] if len(sys.argv) > 2 else [0.002, 0.005, 0.01, 0.02, 0.03, 0.05]:
        c, g = cmin(m, mu)
        out.append(((H-mu)/(1-c), mu, c, g))
        print(f"m={m} mu={mu} c={c:.6f} bound={(H-mu)/(1-c):.7f} argmin gaps={np.round(g,3)}", flush=True)
    print("best", max(out, key=lambda o: o[0])[:3])
