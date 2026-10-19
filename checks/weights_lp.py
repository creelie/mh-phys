"""Weighted windows with analytic gradients and mirror-symmetric weights.
F(g) = sum_i b_i g_i + sum_{r,i} a[r][i] w(y_{i+r}-y_i), sum_i a[r][i] <= 2, sum_i b_i = mu, a, b >= 0.
Cutting planes: LP over the weights against local minimisers of the current F.
usage: python3 -I weights_lp.py R mu iterations bmin_frac procs   (Table I: weights_lp.py 8 0.0025 45 0.2 4)
Every iteration saves the best weights and the cuts to ww_R<R>_mu<mu>.npz (and resumes from it).
b_i >= bmin_frac*mu/n keeps every linear weight positive, which the verifier needs."""
import itertools, os, sys, time, numpy as np
from multiprocessing import Pool
from scipy.optimize import minimize, linprog
s2 = np.sqrt(2.0); H = 1.5 - (np.cos(1/s2)/np.sin(1/s2))/s2; K0 = s2*np.sin(1/s2); t0 = 1/(s2*np.pi)


def dsinc(t):
    t = np.asarray(t, float); out = np.zeros_like(t); nz = np.abs(t) > 1e-12
    out[nz] = (np.cos(np.pi*t[nz]) - np.sinc(t[nz]))/t[nz]; return out


def k(x): return 0.5*(np.sinc(x - t0) + np.sinc(x + t0))/K0
def dk(x): return 0.5*(dsinc(x - t0) + dsinc(x + t0))/K0


R = int(sys.argv[1]); mu = float(sys.argv[2]); iters = int(sys.argv[3]); bmin = float(sys.argv[4]); procs = int(sys.argv[5]); n = R - 1
out = f"ww_R{R}_mu{mu}.npz"
pairs = [(r, i) for r in range(1, R) for i in range(R - r)]
P = len(pairs); pr = np.array([p[0] for p in pairs])
M = np.zeros((P, n))                      # pair distances D = M @ g
for j, (r, i) in enumerate(pairs): M[j, i:i+r] = 1.0


def feats(g): D = M @ g; return np.concatenate([g, k(D)**2])


def fg(g, x):
    D = M @ g; kd = k(D); f = x[:n] @ g + x[n:] @ kd**2
    return f, x[:n] + M.T @ (x[n:] * 2*kd*dk(D))


perm = list(range(n))[::-1] + [n + pairs.index((r, R-1-r-i)) for r, i in pairs]   # mirror image of a variable
starts = [np.array(wd) for wd in itertools.product([1.04, 1.98], repeat=n)]
starts += [np.array(wd) for wd in itertools.product([1.04, 1.98, 3.0], repeat=n) if sum(np.array(wd) == 3.0) == 1]
rng = np.random.default_rng(1)


def one(args):
    s, x = args
    r = minimize(fg, s, args=(x,), jac=True, method="L-BFGS-B", bounds=[(0.3, 60)]*n,
                 options={"ftol": 1e-14, "gtol": 1e-10})
    return r.fun, r.x


def search(x, extra=()):
    best = {}
    for f, g in pool.map(one, [(s, x) for s in list(starts) + list(extra)], chunksize=8):
        key = tuple(np.round(g, 4)); best[key] = min(best.get(key, 1e9), f)
    return sorted((f, g) for g, f in best.items())


pool = Pool(procs)


x0 = np.concatenate([np.full(n, mu/n), 2.0/(R - pr)])
nv = n + P
if os.path.exists(out):
    z = np.load(out); cuts = list(z["cuts"]); best = (float(z["c"]), z["x"]); print(f"resumed: {len(cuts)} cuts, best c={best[0]:.7f}", flush=True)
else:
    b = search(x0); print(f"uniform: c0={b[0][0]:.7f} bound={(H-mu)/(1-b[0][0]):.7f} at {b[0][1]}", flush=True)
    cuts = [np.array(g) for _, g in b[:60]]; cuts += [g[::-1] for g in cuts]; best = (b[0][0], x0)
for it in range(iters):
    A = [np.concatenate([-feats(g), [1.0]]) for g in cuts]
    for r in range(1, R): A.append(np.concatenate([np.zeros(n), (pr == r).astype(float), [0.0]]))
    bub = [0.0]*len(cuts) + [2.0]*(R-1)
    Aeq = [np.concatenate([np.ones(n), np.zeros(P), [0.0]])]; beq = [mu]
    for j in range(nv):                       # symmetric weights
        if perm[j] > j:
            e = np.zeros(nv + 1); e[j] = 1; e[perm[j]] = -1; Aeq.append(e); beq.append(0.0)
    res = linprog(np.concatenate([np.zeros(nv), [-1.0]]), A_ub=A, b_ub=bub, A_eq=Aeq, b_eq=beq,
                  bounds=[(bmin*mu/n, None)]*n + [(0, None)]*P + [(None, None)], method="highs")
    x = res.x[:nv]; c_lp = res.x[-1]
    extra = [np.clip(g + rng.normal(0, 0.08, n), 0.4, 50) for g in cuts[-30:]]
    b = search(x, extra)
    if b[0][0] > best[0]: best = (b[0][0], x.copy())
    print(f"it {it}: LP c={c_lp:.7f} true min={b[0][0]:.7f} bound={(H-mu)/(1-b[0][0]):.7f}"
          f"  best={(H-mu)/(1-best[0]):.7f}", flush=True)
    new = [np.array(g) for f, g in b[:15] if f < c_lp - 1e-10]
    cuts += new + [g[::-1] for g in new]
    np.savez(out, cuts=np.array(cuts), c=best[0], x=best[1])
    if not new: break
c, x = best
print(f"best: c={c:.8f} bound={(H-mu)/(1-c):.8f}")
print("b/mu =", np.round(x[:n]/mu, 4))
for r in range(1, R): print(f"a_{r} =", np.round(x[n:][pr == r], 4))
