# Uncapped block defects: 67.311% of zeros simple and on the line

Result (proved, one computer-checked inequality): liminf N_0^s(T,2T)/N(T,2T) ≥ (3125/3112)(H − 7/3200) = 0.6731133552 (R = 8, μ0 = 1/3200, c0 = 13/3125). A first certified run at μ0 = 1/2000 gave 0.6730728.
Previous records: 0.6730085 (creelie/zeta-simple-zeros, unrefereed) and 0.6725007 (Alpöge–Furman, Theorem D).
It uses the framework of Alpöge–Furman (arXiv:2608.13637) and the stability lemma of zeta-simple-zeros.
The stability lemma is reproved in full in the paper.

## 1. Eigenvalue lemma (proved, `eig_lemma.py`)

Let x_1, ..., x_n be real with |x_i − x_j| ≥ δ1 = 11/20 for i ≠ j. Then the matrix (k(x_i − x_j)) has spectrum in [0, 1/(δ1 K0)] ⊂ [0, 1.97903].
Here K0 = √2 sin(1/√2) = 0.918725.

Proof. Write c*Gc = K0^{-1} ∫_{-1/2}^{1/2} cos(√2 u) |Σ c_i e(u x_i)|² du. Take M(u) = sinc²(π δ1 u). Its Fourier transform is (1/δ1)(1 − |ξ|/δ1)_+, which is supported in [−δ1, δ1], so ∫ M |S|² = (1/δ1) Σ|c_i|². It remains to show M ≥ cos(√2 u) on |u| ≤ 1/2:
- cos x ≤ 1 − x²/2 + x⁴/24;
- sinc² t ≥ 1 − t²/3 + 2t⁴/45 − t⁶/315 for 0 ≤ t ≤ 1 (an alternating series with decreasing terms);
- the difference, divided by u², is a convex quadratic in u². It is negative at u² = 0 (−0.00481) and at u² = 1/4 (−0.0569).

The bound is sharp: δ1-spaced lattices reach 1.97903 (1000 points: 1.979025). The critical spacing is 1/(2K0) = 0.54423.

## 2. Clusters

k is decreasing on [0, 1], because K'(x) = −∫ 2πt cos(√2 t) sin(2πxt) dt < 0. So a gap g < 0.55 has w(g) ≥ w(0.55) = 0.35956.

Pinch every block of m consecutive simple zeros into three kinds of piece:
- runs: maximal pieces with all gaps ≥ 0.55;
- pairs of consecutive zeros joined by a gap < 0.55;
- singletons.

A run has eigenvalues ≤ 1.979 + o(1) < 2, so its defect is its full energy 2Σ_{i<j} w, with no cap. A pair has defect 2w ≥ 0.719.

A cluster of s ≥ 2 zeros gives ⌊s/2⌋ pairs. It costs at most c0·s for its own zeros and c0(R − 1) for the run it cuts. Since 2w⌊s/2⌋ − c0 s ≥ 4w/3 − 2c0 ≥ (R − 1)c0 whenever (R + 1)c0 ≤ 0.479, every block satisfies

  D(G_B) + μ span(B) ≥ c0 (m − R + 1) − o(1),   with μ = (R − 1) μ0.

## 3. Window inequality (certified, `verify/bbp.c`, `../stability/verify/bb.c`, Arb table `w4000.bin`)

F_8(g) = μ0 Σ g + Σ_{r<8} (2/(8−r)) Σ_i w(y_{i+r} − y_i) ≥ c0 for all g ∈ [0, ∞)^7:
- μ0 = 1/3200, c0 = 13/3125 = 0.00416: checked by `bbp` (parallel, fail-safe queue), 10,039,560,518 nodes, 15 min on 2 threads, VERIFIED (`verify/bbp8_3200.log`, rerun with the final code in `verify/bbp8_3200_final.log`). This is the theorem.
- μ0 = 1/2000, c0 = 605/100000: checked by the serial `bb`, 5,440,924,378 nodes, 9 min 57 s, VERIFIED (`verify/bb8_2000.log`).
- Table sha256 8ea9e020… (46000 cells, grid 1/4000).
- The true minimum is about 0.0060703, at gaps (1.044, 1.975, 1.040, 1.972, 1.040, 1.975, 1.044).

## 4. Assembly

Blocks of fixed size m. If span(B) ≥ S := c0(m − R + 1)/μ, the inequality is trivial. Otherwise all pairs are within the fixed distance S, so the kernel limit holds uniformly. Average over the m offsets as in zeta-simple-zeros:

  N_0^s ≥ H N + D(M°) − o(N),   D(M°) ≥ (A0/m) N_0^s − μ (m−1)/m · N − o(N),   A0 = c0(m − R + 1).

As m → ∞ the bound tends to (H − μ)/(1 − c0) = (3125/3112)(H − 7/3200) = 0.6731133552 (R = 8, μ0 = 1/3200, c0 = 13/3125). The block condition (R + 1)c0 = 0.03744 ≤ (4/3)w(0.55) = 0.47942 holds.

## 5. How far this goes

- R = 8, μ0 = 1/3200: true minimum c0 ≈ 0.0041773 (bound 0.673125); certified 0.00416 (bound 0.6731134).
- R = 9 (`lm9.txt`): μ0 = 1/2400 gives c0 ≈ 0.0059432 and bound ≈ 0.673168; certification needs about 40 CPU-hours of branch and bound (planned on a GitHub Actions farm).
- Chain ground state from periodic patterns gives upper limits for the method. `chain.py` truncates the pair sum at 400 neighbours, which lowers the energy by about 5e-5; `chain_all.py` sums 2e5 neighbours and bounds the rest, and uses the monotonicity of (H − μ)/(1 − μḡ − E) in μ for a fixed pattern to cover every μ ≥ 0: 0.67340 at μ = 0.002, 0.67347 at μ = 0.0025, 0.67352 at μ = 0.003, largest 0.673556 near μ = 0.0061 (alternation (1.036, 1.964)), rigorous interval bound below 0.67357, so the method stays below 67.36%.
- Tilting the window to cos(a t) (`tilt.py`): a = 1.5 gains only +0.00004 at R = 7.

## 6. Paper

`paper/main.tex` (REVTeX, 10 pages): "More Than 67.31 Percent of the Zeta Zeros Are Simple and on the Critical Line". Builds clean on TeX Live 2026 (no bad boxes). Repository creelie/mh-phys.
