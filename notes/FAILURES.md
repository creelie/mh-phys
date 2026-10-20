# Riemann hypothesis: source check and engine loop (19 Oct 2026)

Verdict: RH is not closed. Nothing below closes it, and no step of the loop reaches Re s = 1/2.

## 1. The cited results, read from their primary sources

| Claim | Source (copy in `sources/`) | What it proves | Status |
|---|---|---|---|
| Quasi-RH by OpenAI | openai/math family 003, commit adc7f12: "The Quasi-Riemann Hypothesis: A Zero-Free Half-Plane Re s > 7/8" (30 Sep 2026); a shorter 11/12 version (5 Oct); "Uniform exclusion of Landau–Siegel zeros" (1 Oct) | zeta and every Dirichlet L-function have no zeros with Re s > 7/8 (Hecke L over Q(√-3) too); Siegel zeros excluded in the form (1-β) log q ≥ c | Unrefereed manuscripts. OpenAI's `lean/docs/003.md` says the 7/8 bound and the Siegel gap are formalised (comparator statements); not rebuilt here |
| "Highest limit of critical zeros" by Anthropic | Alpöge–Furman, arXiv:2608.13637v2, argument found by Claude; independent proof by Lamzouri, arXiv:2609.02882v2 | at least 2/3 of zeros (0.6725 with the Montgomery–Taylor window) are simple and on the critical line; at least 5/6 (0.8362) distinct | arXiv preprints, not yet refereed; the first is stated to be checked in Lean 4 |
| Higher constants in the zip | preprints 230435 and 230240 (Deep): 68.18% simple on the line, 84.09% distinct, called the ceiling of the bandwidth-one method | proportions only | Preprints.org, not peer reviewed |

Why they do not give RH, even together: the quasi-RH clears Re s > 7/8 only, so the strip 1/2 < Re s ≤ 7/8 is untouched. A proportion of zeros on the line, even 100%, still allows a density-zero set off the line (Bohr–Landau already showed in 1914 that almost all zeros lie within any ε of the line). Neither statement says anything about an individual zero at, say, 0.6 + iγ.

## 2. Loop (failure, rule change, retry)

**Round 1. Combine quasi-RH with the 67% proportion.**
Fails: the two statements control different things (a fixed half-plane, and a proportion), and their conjunction is consistent with a zero at 0.6 + iγ.
Rule: a route must bound every individual zero in 1/2 < σ ≤ 7/8, not a proportion.

**Round 2. Uniform Collar Bound route (zip: 222464, 222943, 223509, 224154, 224427, 224691, 224739, 226914).**
All claim E_2(X,H) = Σ_{H<h≤2H} C_h(X-h)^2 ≪ X^{2+ε}/H and deduce L(X) ≪ X^{1/2+ε}.
Fails: exact values (`notes/loop/ucb_check.py`, output `notes/loop/ucb_check.txt`, FFT autocorrelation of λ up to 4·10^6, collar identity L(X)^2 = X + 2ΣC_h checked exactly) give E_2 ≈ X·H, the square-root size of each C_h. At X = 4·10^6, H = X^0.8 the claimed bound is too small by a factor 8,498, and the ratio grows like H^2/X. The claimed bound is below the true size whenever H > X^{1/2}. Deep's own later preprints 218992 and 229091 already reject this route (fourth-moment floor; uniform |C_h| ≪ Y^θ gives only σ > (1+θ)/2, so even θ = 1/2 reaches only σ > 3/4).
Rule: reject any correlation bound smaller than the random-model size √X per shift; absolute-value bounds on C_h cannot reach the line.

**Round 3. Nyman–Beurling completeness (zip: 224525).**
Fails at Theorem `thm:mellin`: orthogonality to ρ_n for integer n only does not force F(1-s)ζ(s)/s = 0 on the critical line; that would need orthogonality against all real dilations. Completeness of {ρ_n} is equivalent to RH, so the step assumes what it proves. 201997 (Deep et al.) states correctly that its framework does not prove RH.
Rule: a completeness proof must not pass through vanishing of the Mellin transform on the line without a reason specific to the integer dilations.

**Round 4. Chowla route (zip: 210769, 209855).**
Fails: 210769 cites "the proved Chowla correlations"; fixed-shift Chowla is not proved (known: two-point logarithmic Chowla, Tao 2016; odd order, Tao–Teräväinen). Even full Chowla gives only σ > 3/4 by the barrier theorem of 229091. 209855 states its last averaged input as open.
Rule: correlation routes need cancellation across shifts, not shift-by-shift bounds.

**Round 5. Push the 7/8 half-plane down.**
The 7/8 comes from the sextic large sieve and the balanced row estimate of the OpenAI manuscript (§9). Lowering it to 1/2 would need square-root cancellation in that sieve, which is RH-strength for the Hecke family over Q(√-3). No reduction found. A density method (proportion → 1) cannot help either: the trace-moment obstruction in creelie/Riemann-H shows bounded even moments only give β ≤ 1/2 + 1/(bλ) for isolated zeros.
Rule kept open: a proof must exclude each zero in 1/2 < σ ≤ 7/8; no current input does.

## 3. What would count as progress here

- A zero-free half-plane Re s > θ with θ < 7/8 (new, still not RH).
- A proportion above 68.18% on the line (new, still not RH).
Neither is RH, and neither was obtained in this loop.

**Round 6. Push the proportion past the 68.18% ceiling (folder `ceiling/`, notes in `ceiling/NOTES.md`).**
Closed: 230435/230240's 68.18% is the method's ceiling, not a bound; the best window gives 67.25%. The ceiling LP rebuilt (period 8: 0.68600, period 12: 0.68375). Adding three-point (Rudnick–Sarnak range) data lifts it to 0.70880 at period 8 (converged) and about 0.7043–0.7057 at period 12, by local pricing (not certified). The cyclic cubic trace has frequency sum at most 2, so it sits in the Rudnick–Sarnak range.
Not closed: a cubic rank–trace inequality valid for off-line pairs of positive depth. Without it there is no new lower bound.
Rule: higher-order inputs only count once an inequality turns them into a bound on simple zeros for every configuration.

**Round 7. Make the cubic certificate safe against off-line pairs.**
The period-8 dual uses the cubic term with a negative sign, so an off-line pair, whose cubic contribution grows faster than its quadratic one, can break the inequality. Modelling off-line pairs as weighted configurations (`ceiling/lp_depth.py`) fails calibration: with pair data alone it drops to 0.584, below the proved 0.6725, so it is not a valid model of what the zeros allow.
Not closed. Rule: an inequality for higher traces needs an inertia statement for the off-line blocks, not a configuration bound.

**Round 8. Use only the three-point data that known zero-density estimates control.**
An off-line pair at depth y makes the triple term grow like e^{4π(a+b)y/N}, and zero-density estimates with constant c make such pairs a proportion e^{−2πcy}. So only triples with a + b < cN/2 can be used. With c = 1 (Jutila) the ceiling rises by 0.0002 at period 12; with c = 8/7 (Conrey, announced) by 0.0006; with c = 3/2 by 0.0013. The 2-point gain needs c ≈ 2.
Closed, negative. Rule: information beyond pair correlation counts only in the range where off-line zeros are provably rare.

**Round 9. Use the full stability defect of long blocks (scripts in `checks/` and `verify/`, notes in `notes/UNCAPPED.md`).**
Closed, positive. The cap that earlier versions put on the defect of each eigenvalue is not needed: zeros at mutual distance at least 11/20 give an overlap matrix with every eigenvalue below 1/(δ1 K0) = 1.97903 < 2 (sharp, by a Fejér majorant), so runs of well-separated zeros pay their full energy, and short gaps are paid by the pair they form (defect at least 0.719). With an eight-zero window inequality certified in interval arithmetic (10^10 branch-and-bound nodes), liminf N_0^s/N ≥ 0.6731133552, above both 0.6725007 (Alpöge–Furman) and 0.6730085 (zeta-simple-zeros). Periodic chain patterns, summed over all pairs, bound the method below 0.67357 for every weight (`checks/chain_all.py`).
Rule: before certifying a bigger window, compute the chain limit of the method over all pairs (a truncated pair sum is too low by about 5e-5); once the gain left is under 0.0005, information beyond the Gram matrix of simple zeros is needed.

**Round 10. Weight every pair and every gap of the window separately (scripts `checks/weights_lp.py`, `verify/bbp_ww.c`).**
Closed, positive. The window inequality only needs the pair weights at each distance to add up to at most 2 and the linear weights to add up to μ; the weights may depend on the place in the window. A cutting-plane linear programme (maximise the smallest local minimum of F over the weights, adding the new minimisers as cuts) moved the eight-zero value from 0.0046290 (uniform) to 0.0047731 at μ = 0.0025, and spreads the minimum over eight extremal configurations. Rounded down, the weights give F ≥ 19/4000, certified in 21,067,790,614 branch-and-bound nodes, so liminf N_0^s/N ≥ 0.6731985970, above 0.6731134 (uniform eight-zero window) and 0.6731524 (uniform nine-zero window, not finished). The restricted chain energy bounds every weighting of eight zeros by 0.6732159.
Rule: optimise the weights of a window by linear programming before lengthening the window; it gains more than one extra zero and costs only twice the nodes.
Rule: the limit for a window of R points is the chain energy with pairs at most R − 1 apart; compute it before choosing R and μ.

**Round 11. Weighted nine-zero window (scripts `checks/weights_lp.py`, `verify/bbp_ww.c`; outputs `checks/weights_lp9.txt`, `verify/bbpww9.log`).**
Closed, positive. Forty-five rounds of the cutting-plane programme at μ = 0.004 reach 0.6732640, within 10^-6 of the nine-point chain limit 0.6732646. Rounded down, the weights give a minimum of F near 0.0070747, and F ≥ 7/1000 is certified in 52,891,178,741 nodes (30 minutes on four threads), so liminf N_0^s/N ≥ (1000/993)(H − 19999/5000000) = 0.6732133975, above 0.6731986 (weighted eight-zero window). Uniform weights at the same μ give only 0.6731263.
Rule: the node count rises steeply as the certified threshold approaches the true minimum of F; certify first at a margin near 1 percent (half an hour at nine zeros) and tighten afterwards (7.03/1000 would give 0.6732340).
Rule: longer windows gain slowly (limit 0.6732908 for ten zeros, 0.6733252 for twelve, 0.67357 for all); the method is within 0.0004 of its ceiling.

**Round 12. Weighted ten-zero window (Anaxia 3.1.1; scripts `checks/weights_lp.py`, `verify/bbp_ww.c`; outputs `checks/weights_lp10.txt`, `verify/bbpww10.log`).**
Before the round the engine's main branch was read for new construction techniques; the new ones (gadget graphs, CEGAR forcers, reflexive covers, the tensor engine for abelian varieties; and, read again before the second-order checker, the barrier ledger, the atlas, the superstring and blow-up modules to pull request 33) do not act on zero windows, and the cutting-plane loop of `weights_lp.py` is already a CEGAR loop. Closed, positive. Forty-six rounds of the cutting-plane programme at μ = 0.004 reach 0.6732900, within 10^-6 of the ten-point chain limit 0.6732908. Rounded down, the weights give a minimum of F of 0.0071126. With the cell bound alone the target 22/3125 was out of reach in practice: one of 64 shards took 88,467,547,937 nodes and 57 minutes. A second-order Taylor bound on each box, from Arb tables of w and w′ at the grid points and of lower bounds for w″ on cells, cut that shard to 33,813,175 nodes and 8 seconds, and the whole search at the sharper target 889/125000 = 0.007112 to 1,424,333,714 nodes in 2 minutes 48 seconds. The target 0.0071125 fails at a single-cell box, so 0.007112 is the end for these weights at cells of width 1/4000. Hence liminf N_0^s/N ≥ (125000/124111)(H − 7999/2000000) = 0.6732896395, above 0.6732134 (nine zeros). Uniform weights at the same μ give only 0.6731519.
Rule: the search is heavy-tailed across shards (one shard of 64 ran ten times longer than its neighbours); run the shards in one resumable loop that skips logged shards, so that a stopped job loses at most one shard.
