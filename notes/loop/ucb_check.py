"""Round 2: test the Uniform Collar Bound E_2(X,H) << X^{2+eps}/H of the
Liouville-correlation preprints against exact values of E_2, computed by FFT
autocorrelation of the Liouville function."""
import numpy as np, sys

def liouville(N):
    om = np.zeros(N + 1, dtype=np.int8)          # Omega(n) mod 2
    rem = np.arange(N + 1, dtype=np.int64)
    sieve = np.ones(N + 1, dtype=bool); sieve[:2] = False
    for p in range(2, int(N ** 0.5) + 1):
        if sieve[p]:
            sieve[p*p::p] = False
    for p in np.nonzero(sieve)[0]:
        pk = p
        while pk <= N:
            om[pk::pk] ^= 1
            pk *= p
    lam = np.where(om == 1, -1, 1).astype(np.float64); lam[0] = 0
    return lam

def correlations(lam, X):
    a = lam[1:X + 1]
    n = 1 << (2 * X - 1).bit_length()
    F = np.fft.rfft(a, n)
    ac = np.fft.irfft(F * np.conj(F), n)[:X]      # ac[h] = sum_{n<=X-h} lam(n)lam(n+h)
    return np.rint(ac)

for X in [10**5, 10**6, 4 * 10**6]:
    lam = liouville(X)
    C = correlations(lam, X)
    L = lam[1:X + 1].sum()
    assert C[0] == X and L * L == X + 2 * C[1:].sum()       # collar identity check
    for a in [0.5, 0.6, 0.7, 0.8, 0.9]:
        H = int(X ** a)
        E2 = float((C[H + 1:2 * H + 1] ** 2).sum())
        print(f"X={X:>8} H=X^{a}: E2={E2:.3e}  X^2/H={X*X/H:.3e}  "
              f"E2/(X^2/H)={E2/(X*X/H):8.2f}  E2/(XH)={E2/(X*H):.3f}")
    sys.stdout.flush()
