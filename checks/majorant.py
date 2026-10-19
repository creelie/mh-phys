"""Bound lambda_max(G) for delta-separated points via a band-limited majorant.
f(t) = cos(sqrt2 t)/K0 on [-1/2,1/2].  M = (ftilde 1_J) * Fejer_Delta / m_eps >= f 1_I, Mhat supported in [-Delta, Delta].
lambda_max <= Mhat(0) + 2 sum_{n>=1} e(n delta), e(x) = sup_{|u|>=x} |Mhat(u)|."""
import numpy as np
s2 = np.sqrt(2.0); K0 = s2*np.sin(1/s2)
def fej(s, D): return D*np.sinc(D*s)**2           # Fejer kernel, integral 1, FT = triangle on [-D, D]
def bound(delta, D, eps, n_t=6001, n_x=4001):
    t = np.linspace(-0.5-eps, 0.5+eps, n_t); dt = t[1]-t[0]
    ft = np.cos(s2*np.clip(np.abs(t)-eps, 0, None))/K0            # local max envelope of f (f decreasing in |t|)
    s = np.linspace(-eps, eps, 4001); m_eps = np.trapezoid(fej(s, D), s)
    x = np.linspace(0, D, n_x)
    hat = np.array([np.trapezoid(ft*np.cos(2*np.pi*xx*t), dx=dt) for xx in x]) * np.clip(1-x/D, 0, None) / m_eps
    e = np.maximum.accumulate(np.abs(hat)[::-1])[::-1]
    tail = sum(2*np.interp(n*delta, x, e) for n in range(1, int(D/delta)+1))
    return hat[0] + tail, hat[0], tail
best = {}
for delta in [0.6, 0.7, 0.8]:
    for D in [1.5, 2, 3, 4, 6, 8]:
        for eps in [0.05, 0.1, 0.15, 0.2, 0.3]:
            b = bound(delta, D, eps)
            if delta not in best or b[0] < best[delta][0]: best[delta] = (b[0], b[1], b[2], D, eps)
    print(f"delta={delta}: lambda_max <= {best[delta][0]:.4f} (Mhat0={best[delta][1]:.4f}, tail={best[delta][2]:.4f}, Delta={best[delta][3]}, eps={best[delta][4]})", flush=True)
