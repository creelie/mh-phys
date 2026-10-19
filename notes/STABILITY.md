# Raising the proved bound: uncapped block defects

The proved bounds:
- 0.6725007, Alpöge–Furman Theorem D;
- 0.6730085, creelie/zeta-simple-zeros (unrefereed), which uses a stability refinement of the rank–trace lemma, a 7-point local inequality and blocks of 269 zeros.

That argument loses most of its gain in Lemma 4.3, where the defect of a whole 269-zero block is capped at 1.

Uncapped version (`blocks.py`): pinch the Gram matrix of the simple zeros into all windows of m consecutive zeros and use tr Ψ(G) exactly. This gives D(M°) ≥ (1/m) Σ_windows tr Ψ(G_window). If every window satisfies tr Ψ/m + μ·span/(m−1) ≥ c, then the proportion is at least (H − μ)/(1 − c).

| window m | best bound found (multistart, not certified) |
|---|---|
| 3 | 0.672277 (no gain) |
| 4 | 0.672557 |

The minimisers are near-lattice configurations with spacing about 1.04. The kernel k vanishes at 1.057, 2.03 and 3.02, close to the integers.

Upper limit of this method (`lattice.py`): lattice configurations are admissible, and the infinite lattice at spacing 1.0001 has per-zero defect 2W = 0.006136. This caps the method at 0.674762 (μ = 0.00578).

Uncapped 7-point inequality (`f6.py`, same weights 2/(7−r) as zeta-simple-zeros, cap removed): the best is μ0 = 0.0005, c0 = 0.00534, bound 0.673095 (multistart, not certified). The minimisers are subsets of a near-integer lattice (gaps 1.04 and 1.98), because the kernel vanishes near 1.057, 2.03 and 3.02. Dropping the cap needs λ_max(G_B) ≤ 2 for the blocks, which is not yet proved.

Paper rule (Tubai, 15:30 UTC 19 Oct): any paper from this thread must name the BLLM engine (Anaxia, creelie/eBLLM) where it helped the constructions, in a detailed footnote that says exactly what it did. Cite engine release v3.0.5 (10.5281/zenodo.23469389; concept 10.5281/zenodo.23286641), and the BLLM paper where it fits, and name the engine in the data-availability paragraph.

## Longer windows (`fR.py`, multistart, not certified)

| R points | best bound | at μ |
|---|---|---|
| 7 | 0.673095 | 0.003 |
| 10 | 0.673306 | 0.003 |
| 13 | 0.673393 | 0.004 (still rising) |

Inside the capped framework these gains vanish. Blocks must have size m ≤ 1/c0, so each block loses (R−1)/m of its windows, and the 7-point framework with optimised μ0 already sits at 0.67301. The gain therefore needs blocks of unlimited size, which in turn needs λ_max(G_B) ≤ 2.

## The eigenvalue lemma (numerical evidence, `lammax.py`)

The largest eigenvalue of the Gram matrix of 12 zeros with every gap ≥ δ:
- 2.67 at δ = 0.4;
- 2.15 at δ = 0.5;
- 1.96 at δ = 0.55;
- 1.79 at δ = 0.6;
- 1.54 at δ = 0.7.

This matches the lattice value 1.087/δ, so δ* ≈ 0.54.

Plan:
- Runs with all gaps ≥ 0.6 have λ ≤ 2, and their defect is their full energy, with no cap.
- A gap below 0.6 has k ≥ 0.54, so w ≥ 0.29. Clusters split into blocks of 2 or 3 zeros, each with D ≥ 0.58, which pays for the cut windows (at most (R−1)c0 ≈ 0.09 per cut) many times over.

Still to prove: λ ≤ 2 for runs of any length with gaps ≥ 0.6, by a weighted large-sieve or majorant argument plus a certified finite check.

## The eigenvalue lemma by a majorant (`majorant.py`): not enough

A Fejér-type band-limited majorant combined with a Gershgorin tail gives λ_max ≤ 2.50 for δ = 0.6, 2.33 for δ = 0.7 and 2.28 for δ = 0.8. None of these reaches 2. The numerical truth is about 1.087/δ (1.52 at δ = 0.8). The classical large sieve constant 1 + 1/δ is too weak when δ < 1. The lemma stays open.

## Without the lemma: the capped framework with longer windows

Blocks of m = R − 1 + ⌊1/c0⌋ zeros keep A0 = c0(m − R + 1) ≤ 1, so the cap never bites, and the bound becomes (H − μ(m−1)/m)/(1 − A0/m). Using the multistart minima:
- R = 7 at μ = 0.002: 0.673025;
- R = 10 at μ = 0.002: 0.673132;
- R = 13 at μ = 0.002: 0.673152.

These need no new lemma, only a certified R-point inequality. That is the verifier of zeta-simple-zeros extended from 6 to R − 1 dimensions.

## Certification run (`verify/bb.c`, `verify/table.py`)

The C branch-and-bound uses a rigorous Arb table of lower bounds for w on cells of width 1/4000 and downward rounding. It reproduces the zeta-simple-zeros 7-point inequality F ≥ 19/5000 at μ0 = 1/3000 in 7.7 s (89 million boxes; the same three gap components [3809,4778], [7221,9363], [10572,44827]).

For R = 8 at μ0 = 1/2800 and target 0.0047, it fails at gaps (1.046, 1.977, 1.045, 1.996, 1.996, 1.991, 1.048). A local minimisation from there reaches F = 0.004684, below the multistart value 0.004746. So the multistart minima overstate c0. With the true value the capped bound for R = 8 is about 0.673043, a gain of 0.00003. **The capped framework is saturated at about 67.30%.**

Without the cap (which needs the eigenvalue lemma), R = 8 would give about 0.67315, and larger R somewhat more, below the 67.48% lattice limit.
