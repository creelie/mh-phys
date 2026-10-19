"""Upper limit of the block-defect method: lattice configurations are allowed, so c(mu) <= min_g 2W(g) + mu g."""
import sys, numpy as np
sys.path.insert(0, sys.argv[1] if len(sys.argv) > 1 else ".")
from blocks import k, H
g = np.linspace(0.3, 40, 200001)
W = np.zeros_like(g)
for n in range(1, 3000): W += k(n*g)**2
loc = (g > 0.9) & (g < 1.3); i = np.argmin(np.where(loc, W, np.inf))
print("local lattice min near g =", g[i], " 2W =", 2*W[i])
best = (0, 0, 0)
for mu in np.linspace(1e-5, 0.03, 3000):
    v = 2*W + mu*g; j = np.argmin(v); c = v[j]; b = (H-mu)/(1-c)
    if b > best[0]: best = (b, mu, c, g[j])
print("lattice upper limit: bound =", best[0], " mu =", best[1], " c =", best[2], " at g =", best[3])
