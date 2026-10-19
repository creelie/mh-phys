"""Rigorous binary64 lower bounds for w(x) = k(x)^2 on the cells [i/G, (i+1)/G], 0 <= i < n, by Arb ball arithmetic.
k(x) = K(x)/K(0), K(x) = int_{-1/2}^{1/2} cos(sqrt2 t) cos(2 pi x t) dt = (1/2)[sinc((sqrt2 - 2pi x)/2) + sinc((sqrt2 + 2pi x)/2)].
Usage: python3 table.py G n out.bin"""
import sys, math, struct, hashlib
from flint import arb, ctx, fmpq
ctx.prec = 128
G, n, out = int(sys.argv[1]), int(sys.argv[2]), sys.argv[3]
s2 = arb(2).sqrt(); pi = arb.pi(); k0 = s2*(1/s2).sin()
def cell(i): return arb(fmpq(2*i+1, 2*G), fmpq(1, 2*G))
vals = []
for i in range(n):
    x = cell(i); f = 2*pi*x
    k = (((s2 - f)/2).sinc() + ((s2 + f)/2).sinc())/2/k0
    lo = float(k.abs_lower())
    lo = 0.0 if lo <= 0 else math.nextafter(lo, -math.inf)
    v = 0.0 if lo == 0 else max(0.0, math.nextafter(lo*lo, -math.inf))
    vals.append(v)
data = b"".join(struct.pack("<d", v) for v in vals)
open(out, "wb").write(data)
print("cells", n, "grid", G, "sha256", hashlib.sha256(data).hexdigest())
