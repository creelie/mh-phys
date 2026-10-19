"""Ceiling programme: min expected simple fraction over N-periodic marked
configurations (marks 1 = simple, 2 = double / depth-0 off-line pair) whose
form factor matches E|p_j|^2 = j (0<j<N), optionally also the three-point data
E[p_a p_b conj(p_{a+b})] = 0 for a+b<N (CUE / Rudnick-Sarnak range).
Column generation; pricing by multistart local optimisation."""
import numpy as np, sys
from scipy.optimize import linprog, minimize
rng = np.random.default_rng(1)

def feats(x, m, N, triples):
    j = np.arange(1, N)
    E = np.exp(2j * np.pi * np.outer(j, x) / N) @ m           # p_j, j=1..N-1
    f = [np.abs(E) ** 2]
    if triples:
        P = np.concatenate([[m.sum()], E])                       # P[k] = p_k
        t = np.array([P[a] * P[b] * np.conj(P[a + b]) for a, b in triples])
        f += [t.real, t.imag]
    return np.concatenate(f)

def simple_frac(m, N): return np.sum(m == 1) / N

def solve(N, use_triples, iters=400, starts=60, frac=1.0):
    triples = [(a, b) for a in range(1, N) for b in range(a, N) if a + b < frac * N] if use_triples else []
    target = np.concatenate([np.arange(1, N), np.zeros(2 * len(triples))])
    cols, costs = [], []
    def add(x, m):
        cols.append(feats(x, m, N, triples)); costs.append(simple_frac(m, N))
    # seed: CUE samples (all simple) and random configurations with doubles
    for _ in range(200):
        A = rng.normal(size=(N, N)) + 1j * rng.normal(size=(N, N))
        Q, R = np.linalg.qr(A); Q = Q * (np.diag(R) / abs(np.diag(R)))
        add(np.angle(np.linalg.eigvals(Q)) % (2 * np.pi) * N / (2 * np.pi), np.ones(N))
    for d in range(1, N // 2 + 1):
        for _ in range(20):
            k = N - d; m = np.ones(k); m[:d] = 2
            add(rng.uniform(0, N, k), m)
    best = None
    for it in range(iters):
        A = np.array(cols).T
        Aeq = np.vstack([np.ones(len(cols)), A]); beq = np.concatenate([[1], target])
        m_ = Aeq.shape[0]; S = np.hstack([np.eye(m_), -np.eye(m_)]); S[0] = 0
        res = linprog(np.concatenate([costs, 100 * np.ones(2 * m_)]), A_eq=np.hstack([Aeq, S]),
                      b_eq=beq, bounds=(0, None), method="highs")
        slack = res.x[len(cols):].sum() if res.status == 0 else None
        if res.status != 0:
            print("LP status", res.status, res.message); return None
        y = res.eqlin.marginals; val = res.fun
        y0, yv = y[0], y[1:]
        # pricing: minimise reduced cost s(c) - y0 - yv.f(c)
        found = 0; minrc = 0.0
        for d in range(0, N // 2 + 1):
            k = N - d; m = np.ones(k); m[:d] = 2
            obj = lambda x: simple_frac(m, N) - y0 - yv @ feats(x, m, N, triples)
            for _ in range(max(2, starts // (N // 2 + 1))):
                r = minimize(obj, rng.uniform(0, N, k), method="L-BFGS-B")
                minrc = min(minrc, r.fun)
                if r.fun < -1e-7:
                    add(r.x % N, m); found += 1
        print(f"N={N} triples={bool(triples)} it={it} value={val:.6f} lagrangian_lb~{val+minrc:.6f} cols={len(cols)} new={found} slack={slack:.2e}", flush=True)
        np.savez(f"dual_N{N}_t{int(bool(triples))}_f{frac}.npz", y=y, val=val, minrc=minrc, triples=np.array(triples))
        if found == 0: break
        best = val
    return val

if __name__ == "__main__":
    N = int(sys.argv[1]); t = sys.argv[2] == "1"
    solve(N, t, frac=float(sys.argv[3]) if len(sys.argv) > 3 else 1.0)
