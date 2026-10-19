"""Ceiling LP with off-line pairs of positive depth.
A pair 1/2 +- delta + i*gamma enters frequency j with weight 2*cosh(j*eta), eta = depth in window units.
Depths are capped at H (Selberg's density theorem makes deeper pairs exponentially rare)."""
import sys, numpy as np
from scipy.optimize import linprog, minimize
rng = np.random.default_rng(1)

def feats(x, eta, d, N, triples):
    j = np.arange(0, N)
    W = np.ones((N, len(x)))
    if d: W[:, :d] = 2*np.cosh(np.outer(j, eta))
    P = np.sum(W*np.exp(2j*np.pi*np.outer(j, x)/N), axis=1)
    f = [np.abs(P[1:])**2]
    if triples:
        t = np.array([P[a]*P[b]*np.conj(P[a+b]) for a, b in triples]); f += [t.real, t.imag]
    return np.concatenate(f)

def solve(N, use_triples, H, iters=400, starts=60):
    triples = [(a, b) for a in range(1, N) for b in range(a, N) if a+b < N] if use_triples else []
    target = np.concatenate([np.arange(1, N), np.zeros(2*len(triples))])
    cols, costs = [], []
    def add(x, eta, d):
        cols.append(feats(x, eta, d, N, triples)); costs.append((N-2*d)/N)
    for _ in range(200):
        A = rng.normal(size=(N, N))+1j*rng.normal(size=(N, N))
        Q, R = np.linalg.qr(A); Q = Q*(np.diag(R)/abs(np.diag(R)))
        add(np.angle(np.linalg.eigvals(Q)) % (2*np.pi)*N/(2*np.pi), np.zeros(0), 0)
    for d in range(1, N//2+1):
        for _ in range(20): add(rng.uniform(0, N, N-d), rng.uniform(0, H, d), d)
    for it in range(iters):
        A = np.array(cols).T
        Aeq = np.vstack([np.ones(len(cols)), A]); beq = np.concatenate([[1], target])
        m_ = Aeq.shape[0]; S = np.hstack([np.eye(m_), -np.eye(m_)]); S[0] = 0
        res = linprog(np.concatenate([costs, 100*np.ones(2*m_)]), A_eq=np.hstack([Aeq, S]), b_eq=beq, bounds=(0, None), method="highs")
        if res.status != 0: print("LP status", res.status, res.message); return
        slack = res.x[len(cols):].sum(); y = res.eqlin.marginals; val = res.fun; y0, yv = y[0], y[1:]
        found = 0; minrc = 0.0
        for d in range(0, N//2+1):
            k = N-d
            def obj(z): return (N-2*d)/N - y0 - yv@feats(z[:k], z[k:], d, N, triples)
            for _ in range(max(2, starts//(N//2+1))):
                z0 = np.concatenate([rng.uniform(0, N, k), rng.uniform(0, H, d)])
                r = minimize(obj, z0, method="L-BFGS-B", bounds=[(None, None)]*k+[(0, H)]*d)
                minrc = min(minrc, r.fun)
                if r.fun < -1e-7: add(r.x[:k] % N, r.x[k:], d); found += 1
        print(f"N={N} H={H} triples={bool(triples)} it={it} value={val:.6f} lb~{val+minrc:.6f} cols={len(cols)} new={found} slack={slack:.1e}", flush=True)
        np.savez(f"dual_depth_N{N}_H{H}_t{int(use_triples)}.npz", y=y, val=val, minrc=minrc)
        if found == 0: break

if __name__ == "__main__":
    solve(int(sys.argv[1]), sys.argv[2] == "1", float(sys.argv[3]))
