"""H(W) = 1 - int w(x)(1 - sinc^2(pi x)) dx, w = (What/What(0))^2, for W(t) = cos(a t) on [-1/2,1/2] (a = sqrt2 is Montgomery-Taylor).
Computed in Fourier space: int w (1-sinc^2) = int (V(al) - V*Tri(al)) with V = (W*W)/K0^2 supported on [-1,1]."""
import numpy as np
from scipy.integrate import quad
def H(a, M=4000):
    K0 = quad(lambda t: np.cos(a*t), -0.5, 0.5)[0]
    k = lambda x: quad(lambda t: np.cos(a*t)*np.cos(2*np.pi*x*t), -0.5, 0.5)[0]/K0
    # x-space integral with tail: w ~ O(1/x^2) and (1 - sinc^2) -> 1; integrate to X and add tail estimate by large X
    X = 400.0; xs = np.linspace(0, X, 400001)
    if abs(a) < 1e-12:
        kv = np.sinc(xs)
    else:
        A = np.pi*xs - a/2; B = np.pi*xs + a/2
        kv = 0.5*(np.sinc(A/np.pi) + np.sinc(B/np.pi))/(K0)
    f = kv**2*(1 - np.sinc(xs)**2)
    return 1 - 2*np.trapezoid(f, xs)
for a in [0.0, 1.0, np.sqrt(2), 1.6, 2.0, 2.5]:
    print(f"a={a:.4f} H={H(a):.7f}")
