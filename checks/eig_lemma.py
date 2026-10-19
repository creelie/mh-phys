"""Eigenvalue lemma: for real x_1..x_n with gaps >= d1, spec(k(x_i-x_j)) lies in [0, 1/(d1*K0)].
Majorant M(u) = sinc^2(pi d1 u) >= cos(sqrt2 u) on |u| <= 1/2; Mhat = (1/d1) max(0, 1-|xi|/d1).
Checks: (a) the Taylor-polynomial proof of the majorant inequality (exact rationals + Arb for pi),
(b) a floating-point cross-check, (c) numerical lambda_max of random and lattice configurations."""
import sys, numpy as np
from flint import arb, ctx, fmpq
ctx.prec = 128
d1 = fmpq(11, 20)
pi = arb.pi(); s2 = arb(2).sqrt(); K0 = s2*(1/s2).sin()
print("K0 =", K0, " bound 1/(d1 K0) =", 1/(arb(d1)*K0), " critical d0 = 1/(2K0) =", 1/(2*K0))
# (a) cos x <= 1 - x^2/2 + x^4/24 and sinc^2 t >= 1 - t^2/3 + 2t^4/45 - t^6/315 (alternating, decreasing terms for t <= 1).
a = (pi*d1)**2/3; b = 2*(pi*d1)**4/45; c = (pi*d1)**6/315
# need (a-1) + (1/6 - b) s + c s^2 <= 0 for s = u^2 in [0, 1/4]; convex in s, so check s = 0 and s = 1/4
e0 = a - 1; e1 = (a - 1) + (arb(1)/6 - b)/4 + c/16
print("(a) endpoint values (must be < 0):", e0, e1, " t_max = pi d1/2 =", pi*d1/2)
assert e0 < 0 and e1 < 0 and pi*d1/2 < 1
# (b) floating-point cross-check of (a): (sinc^2(pi d1 u) - cos(sqrt2 u))/u^2 on a grid of (0, 1/2]
u = np.linspace(1e-3, 0.5, 100000); dd = float(d1.p)/float(d1.q)
r = (np.sinc(dd*u)**2 - np.cos(np.sqrt(2)*u))/u**2
print("(b) min of (sinc^2 - cos)/u^2 on (0, 1/2]:", r.min(), "(positive means the majorant holds)")
# (c) numerics
def k(x):
    x = np.asarray(x, float); A = np.pi*x - 1/np.sqrt(2); B = np.pi*x + 1/np.sqrt(2)
    return 0.5*(np.sinc(A/np.pi) + np.sinc(B/np.pi))/(np.sqrt(2)*np.sin(1/np.sqrt(2)))
rng = np.random.default_rng(0); best = 0
for trial in range(3000):
    m = rng.integers(2, 60); g = 0.55 + rng.exponential(rng.choice([0.01, 0.1, 0.5]), m-1)
    y = np.concatenate([[0], np.cumsum(g)]); lam = np.linalg.eigvalsh(k(y[:, None]-y[None, :]))[-1]
    best = max(best, lam)
for m in [10, 100, 1000]:
    y = 0.55*np.arange(m); print(f"(c) lattice 0.55, {m} points: lambda_max = {np.linalg.eigvalsh(k(y[:,None]-y[None,:]))[-1]:.6f}")
print(f"(c) random 0.55-separated configurations (3000): max lambda_max = {best:.6f}")
print("k(0.55) =", k(0.55), " w(0.55) =", k(0.55)**2)
