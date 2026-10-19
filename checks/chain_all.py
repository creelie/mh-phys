"""Upper limit of the uncapped method over every weight mu >= 0.
For a periodic gap pattern g (period p, mean gbar) and mu >= 0, e(g, mu) = mu*gbar + E(g), with
E(g) = (2/p) sum_{i<=p} sum_{j>i} w(y_j - y_i) (all pairs). Any pattern gives c0 <= e(g, mu), so the
bound of the method at weight mu is at most (H - mu)/(1 - e(g, mu)).
For a fixed pattern, f(mu) = (H - mu)/(1 - mu*gbar - E) has derivative of constant sign (that of gbar*H + E - 1),
so on [mu_a, mu_b] each pattern with e(g, mu_b) < 1 gives bound <= max(f(mu_a), f(mu_b)), with no loss; the best of the patterns
found at the grid points is used on each interval (patterns rounded to 4 decimals, E recomputed).
E(g) is summed to N neighbours and the rest is bounded by |k(x)| <= 1/(K0 (pi x - 1/sqrt2)) for x > 1,
so the printed E is an upper bound up to floating-point rounding.
usage: python3 -I chain_all.py"""
import itertools, numpy as np
from scipy.optimize import minimize
s2 = np.sqrt(2.0); H = 1.5 - (np.cos(1/s2)/np.sin(1/s2))/s2; K0 = s2*np.sin(1/s2)
def k(x): return 0.5*(np.sinc((np.pi*x - 1/s2)/np.pi) + np.sinc((np.pi*x + 1/s2)/np.pi))/K0
def E_trunc(g, reach):
    p = len(g); y = np.concatenate([[0.0], np.cumsum(np.tile(g, reach//p + 2))]); e = 0.0
    for i in range(p):
        d = y[i+1:i+1+reach] - y[i]; e += 2*np.sum(k(d)**2)/p
    return e
def E_upper(g, N=200000):
    p = len(g); gbar = np.mean(g); e = E_trunc(g, N)
    a = np.pi*(N - p)*gbar - 1/s2                      # distance of the (N+1)-th neighbour is >= (N-p+1)*gbar
    return e + 2/(K0**2*np.pi*gbar*a)                  # sum_{m>N} 1/(K0 (pi (m-p) gbar - 1/sqrt2))^2 <= integral
def best(mu):
    b = (np.inf, None)
    for p in range(1, 5):
        for word in itertools.product([1.04, 1.98, 3.0], repeat=p):
            r = minimize(lambda g: mu*np.mean(g) + E_trunc(g, 400), np.array(word), method="L-BFGS-B", bounds=[(0.3, 60)]*p)
            if r.fun < b[0]: b = (r.fun, r.x)
    return b[1]
grid = list(np.round(np.arange(0.0, 0.0201, 0.0005), 6)) + list(np.round(np.arange(0.025, 0.1001, 0.005), 6)) + [0.15, 0.2, 0.3, 0.4, 0.5, 0.6, round(H, 6)]
pats = {}
for m in [m for m in grid[1:] if m <= 0.3]:          # beyond 0.3 the lattice of gap 1 already makes the bound decrease
    g = tuple(float(x) for x in np.round(best(m), 4)); pats[g] = (np.mean(g), E_upper(np.array(g)))
    print(f"mu={m:.4f}: pattern={g} gbar={pats[g][0]:.4f} E<={pats[g][1]:.7f}", flush=True)
worst = 0.0
f = lambda m, gbar, E: (H - m)/(1 - m*gbar - E)
for a, b in zip(grid[:-1], grid[1:]):
    ok = [(max(f(a, *v), f(b, *v)), g) for g, v in pats.items() if 1 - b*v[0] - v[1] > 0]   # c0 <= e < 1 needed
    lim, g = min(ok) if ok else (1.0, None)
    worst = max(worst, lim)
    print(f"mu in [{a:.4f},{b:.4f}]: bound<={lim:.7f} with pattern {g}", flush=True)
print(f"max over all mu: {worst:.7f}")
