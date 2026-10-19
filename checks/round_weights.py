"""Round the weights of a weight search to rationals for the verifier bbp_ww.
usage: python3 -I round_weights.py ww_R<R>_mu<mu>.npz R out.txt da db   (denominators 10^da, 10^db; Table I: da = 5, db = 7)
a_{r,i} = floor(a*10^da)/10^da, b_i = floor(b*10^db)/10^db, so sum_i a_{r,i} <= 2 and sum_i b_i <= mu.
Prints mu_eff = sum b_i as a fraction."""
import sys
from fractions import Fraction
import numpy as np
z = np.load(sys.argv[1]); R = int(sys.argv[2]); n = R - 1
x = z["x"]
pairs = [(r, i) for r in range(1, R) for i in range(R - r)]
DA, DB = 10**int(sys.argv[4]), 10**int(sys.argv[5])
lines, mu = [], Fraction(0)
for i in range(n):
    num = int(np.floor(x[i]*DB)); lines.append(f"b {i+1} {num} {DB}"); mu += Fraction(num, DB)
tot = {}
for j, (r, i) in enumerate(pairs):
    num = int(np.floor(x[n+j]*DA)); lines.append(f"a {r} {i} {num} {DA}"); tot[r] = tot.get(r, 0) + num
assert all(v <= 2*DA for v in tot.values()), tot
open(sys.argv[3], "w").write("\n".join(lines) + "\n")
print("mu_eff =", mu, "=", float(mu), "  sums a_r:", {r: v/DA for r, v in tot.items()}, "  c =", float(z["c"]))
