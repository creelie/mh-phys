"""Rigorous binary64 tables for the second-order bound of bbp_cf, by Arb ball arithmetic at 128 bits.
w(x) = k(x)^2, k(x) = K(x)/K(0), K(x) = int_{-1/2}^{1/2} cos(sqrt2 t) cos(2 pi x t) dt = [s(u_-) + s(u_+)]/2 with
s(u) = sin(u)/u and u_(+-) = (2 pi x +- sqrt2)/2, so K' = (pi/2)[s'(u_-) + s'(u_+)] and K'' = (pi^2/2)[s''(u_-) + s''(u_+)];
s' and s'' are summed from their Taylor series when |u| <= 1.  Written as little-endian doubles, m = 0..L, c = 0..L-1:
  lower bounds of w(m/G); lower and upper bounds of w'(m/G)/G; lower bounds of w''/G^2 on the cell [c/G, (c+1)/G],
the last as the least of the bounds on S equal sub-cells.
Usage: python3 dtable.py G L S out.bin"""
import sys, math, struct, hashlib
from flint import arb, ctx, fmpq
ctx.prec = 128
G, L, S, out = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3]), sys.argv[4]
s2 = arb(2).sqrt(); pi = arb.pi(); K0 = s2*(1/s2).sin()
def series(u, d):               # d-th derivative of sin(u)/u for |u| <= 1: sum_m (-1)^m u^(2m-d) (2m)!/((2m-d)! (2m+1)!)
    s = arb(0)
    for m in range(1, 32):
        c = math.factorial(2*m)//math.factorial(2*m - d)
        s += (-1)**m * c * u**(2*m - d) / math.factorial(2*m + 1)
    return s + arb(0, fmpq(70*70, math.factorial(64)))   # alternating, decreasing terms: the tail is below the first omitted term
def sd(u, d):
    if u.abs_upper() <= 1: return series(u, d)
    sn, cs = u.sin(), u.cos()
    if d == 1: return (u*cs - sn)/(u*u)
    return ((2 - u*u)*sn - 2*u*cs)/(u*u*u)
def ders(x):
    f = 2*pi*x; um, up_ = (f - s2)/2, (f + s2)/2
    K = (um.sinc() + up_.sinc())/2; K1 = (pi/2)*(sd(um, 1) + sd(up_, 1)); K2 = (pi*pi/2)*(sd(um, 2) + sd(up_, 2))
    return K, 2*K*K1/(K0*K0)/G, 2*(K1*K1 + K*K2)/(K0*K0)/(G*G)
def down(x): return math.nextafter(float(x), -math.inf)
def upr(x): return math.nextafter(float(x), math.inf)
wp, gl, gh, hl = [], [], [], []
for m in range(L + 1):
    K, d1, d2 = ders(arb(fmpq(m, G)))
    a = float((K/K0).abs_lower()); a = 0.0 if a <= 0 else math.nextafter(a, -math.inf)
    wp.append(0.0 if a == 0 else max(0.0, math.nextafter(a*a, -math.inf)))
    gl.append(down(d1.lower())); gh.append(upr(d1.upper()))
for c in range(L):
    lo2 = math.inf
    for t in range(S):
        K, d1, d2 = ders(arb(fmpq(2*S*c + 2*t + 1, 2*S*G), fmpq(1, 2*S*G)))
        lo2 = min(lo2, down(d2.lower()))
    hl.append(lo2)
data = b"".join(struct.pack("<d", v) for v in wp + gl + gh + hl)
open(out, "wb").write(data)
print("points", L + 1, "cells", L, "grid", G, "subcells", S, "sha256", hashlib.sha256(data).hexdigest())
