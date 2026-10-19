"""Window tilt: W_a(t) = cos(a t) on [-1/2,1/2]. H(a) = 2 - [int W^2 + iint W(s)W(t)|s-t|]/(int W)^2 (exact),
c0(a, mu0) = min of the R-window functional F_R (structured local search), uncapped bound (H - (R-1)mu0)/(1 - c0).
usage: python3 -I tilt.py R a1,a2,... den1,den2,..."""
import sys, itertools, numpy as np
from scipy.integrate import quad, dblquad
from scipy.optimize import minimize, brentq
def Hval(a):
    W = lambda t: np.cos(a*t)
    I1 = quad(W, -0.5, 0.5)[0]; I2 = quad(lambda t: W(t)**2, -0.5, 0.5)[0]
    # iint W(s)W(t)|s-t| = 2 int_{-1/2}^{1/2} W(s) int_{-1/2}^{s} W(t)(s-t) dt ds
    inner = lambda s: quad(lambda t: W(t)*(s - t), -0.5, s)[0]
    I3 = 2*quad(lambda s: W(s)*inner(s), -0.5, 0.5, epsabs=1e-13, epsrel=1e-13)[0]
    return 2 - (I2 + I3)/I1**2, I1
R = int(sys.argv[1])
for a in [float(v) for v in sys.argv[2].split(',')]:
    H, K0 = Hval(a)
    def k(x): return 0.5*(np.sinc((np.pi*x - a/2)/np.pi) + np.sinc((np.pi*x + a/2)/np.pi))/K0
    z = [brentq(k, n - 0.3 + 1e-9, n + 0.9) for n in (1, 2, 3)] if k(1.0) > 0 else [np.nan]*3
    def F(g, mu0):
        y = np.concatenate([[0.0], np.cumsum(g)]); v = mu0*np.sum(g)
        for r in range(1, R): v += 2.0/(R-r)*np.sum(k(y[r:]-y[:-r])**2)
        return v
    starts = [z[0]-0.02, z[1]-0.05, z[2]]
    for den in [int(d) for d in sys.argv[3].split(',')]:
        mu0 = 1.0/den; best = (np.inf, None)
        for word in itertools.product(starts, repeat=R-1):
            r = minimize(F, np.array(word), args=(mu0,), method="L-BFGS-B", bounds=[(0.55, 40)]*(R-1))
            if r.fun < best[0]: best = (r.fun, r.x)
        c0 = best[0]; mu = (R-1)*mu0
        print(f"R={R} a={a:.4f} H={H:.7f} zeros={np.round(z,4)} mu0=1/{den} c0={c0:.7f} bound={(H-mu)/(1-c0):.7f} gaps={np.round(best[1],3)}", flush=True)
