# Beyond the 68.18% ceiling: zeros simple and on the critical line

Target: an unconditional proportion of zeros of ζ that are simple and on the critical line, above the bandwidth-one ceiling 0.6818287 of Alpöge–Furman (arXiv:2608.13637, §7.2).

## State of the art (checked 19 Oct 2026)

| Bound | Source | Status |
|---|---|---|
| 0.6725 simple on the line, 0.8362 distinct | Alpöge–Furman 2608.13637 (Claude); Lamzouri 2609.02882 | preprints, Lean for AF |
| 0.67300853 | creelie/zeta-simple-zeros (GPT-generated, 7 consecutive zeros) | unrefereed |
| 0.6818 (ceiling) | AF §7.2: an upper limit for every certificate that uses only pair correlation in support [-1,1] | proved limit, not a lower bound |

**Correction to preprint 230435 / 230240 (Deep).** These claim 68.18% simple on the line and 84.09% distinct, using a window with R(ψ) = 1.3182. No window does that: the minimum of R(ψ) over all windows on [-1/2,1/2] is the Montgomery–Taylor value 1.3274993 (checked by solving the discretised quadratic problem, n = 2000: 1.3274992), so 2 − R ≤ 0.6725. The 0.6818 is AF's ceiling (an upper limit for the method), not a bound that has been proved. The correct constants from that argument are 0.6725 and 0.8362.

## Step 1 (closed, computation): rebuild the ceiling

`lp.py` minimises the expected simple fraction over N-periodic configurations of N unit masses, with marks 1 (simple) and 2 (double, which is the same as an off-line pair of depth 0), subject to E|p_j|² = j for 0 < j < N (Montgomery's form factor). It uses column generation with multistart pricing.

| N | pair data only |
|---|---|
| 8 | 0.68600 |
| 12 | 0.68375 |
| 16 | 0.68332 (three-point: about 0.706, stopped before convergence) |
| 256 | 0.68183 (AF, exact LP) |

## Step 2 (closed, computation): add the three-point data

The extra constraint is E[p_a p_b conj p_{a+b}] = 0 for a + b < N. This is the CUE value: it follows from Schur orthogonality, because every partition of a + b ≤ N is admissible. The frequency vector (a, b, −(a+b))/N has Σ|ξ| = 2(a+b)/N < 2, which is the Rudnick–Sarnak range.

| N | pair only | pair + three-point |
|---|---|---|
| 8 | 0.68600 | **0.70880** |
| 12 | 0.68375 | **0.70497** (Lagrangian bound 0.70495) |

So the three-point data excludes the extremal configurations for the pair data, and raises the ceiling by about two points. Information inside the Rudnick–Sarnak range therefore goes past the 68.18% barrier. The lower values in the logs come from local pricing, so they are not certified.

## Step 3 (closed, observation; to be written carefully): the cubic trace is in range

AF §7.2(e) say that tr G̃^k at X ≍ T is available only for k = 1. For the cyclic product tr G̃³ = Σ K(γ1−γ2)K(γ2−γ3)K(γ3−γ1), whose window frequencies are u_i ∈ [−1/2, 1/2], the frequencies ξ_i = u_i − u_{i−1} satisfy Σ|ξ_i| = 2(max u − min u) ≤ 2. That is the same boundary as the Hilbert–Schmidt norm (|α| ≤ 1). For k = 4 the sum reaches 4, which is why AF's HL*(4) needs prime pairs. So the third moment sits in the Rudnick–Sarnak range, with the same edge issue as the pair case. This needs (i) the Rudnick–Sarnak evaluation for sums over all zeros, including complex ordinates, and (ii) the edge Σ|ξ| = 2, handled the way Montgomery's |α| = 1 is.

## Step 4 (not closed): a cubic rank–trace inequality

To turn Steps 2–3 into a theorem, one needs an inequality

  rank P1 ≥ a1 tr G̃ + a2 tr G̃² + a3 tr G̃³ − c·(s2 + p)

valid for every configuration, including off-line pairs of positive depth. AF's quadratic version works by Sylvester inertia. For a cubic, a (1,1) block with a large negative eigenvalue −ν contributes −a3 ν³, so a3 < 0 is unsafe without a bound on negative eigenvalues. A positive a3 fails for large positive eigenvalues. The Hilbert–Schmidt evaluation bounds Σν² ≤ (4/3)N but not max ν. This is the open step. Possible ways in:
- the dual certificate of the N = 8 three-point LP: check the sign it puts on the cubic constraint;
- capping eigenvalues using N(t, t+1) ≪ log t, together with the zero-density estimate for deep off-line zeros.

## Step 5 (closed, computation): the sign of the cubic term in the certificate

The dual of the period-8 three-point LP (`dual_N8_t1.npz`, value 0.708804, priced to 1e-7) reads

  s ≥ 1 − Σ_j c_j |p_j|² − Σ_{a,b} w_{ab} Re(p_a p_b conj p_{a+b}),

with every c_j ≥ 0 (except c_7 = −0.0003) and every w_{ab} > 0 (from 0.0024 to 0.0232, total 0.117). So the certificate uses the cubic term with a negative sign: this is the unsafe sign of Step 4. An off-line pair of depth η enters frequency j with weight 2cosh(jη), and its cubic contribution grows like cosh³, faster than the quadratic one. The 70.9% therefore survives only if deep pairs are controlled. `lp_depth.py` adds pairs of depth up to a cap H as configurations; Selberg's density theorem makes pairs of depth ≥ η (in units of 1/log T) a proportion ≪ e^{−cη}, which is what a cap would have to be paid for.

## Step 6 (closed, negative): off-line pairs of positive depth

`lp_depth.py` at period 8 with three-point data: depth ≤ 0.5 gives 0.7013 (still falling); depth ≤ 2 falls below 0.621 before the solver fails. Calibration with pair data only at depth ≤ 2 gives 0.584, which is below the proved 0.6725. So the configuration model with cosh weights admits configurations that the real zeros exclude (Alpöge–Furman handle off-line pairs through the signature of the Hermitian form, not through configurations). The depth runs are therefore not evidence either way. What stands is Step 5: the certificate needs the cubic term with the sign that deep off-line pairs exploit, so a proof needs an inertia argument for cubic forms, which Sylvester's law does not give.

## Next (open): Hölder in the safe direction

For a positive semidefinite part with eigenvalues λ ≥ 0, Hölder gives rank ≥ (Σλ²)³/(Σλ³)², which needs an upper bound on the cubic trace, the safe direction. An off-line pair whose block is traceless has eigenvalues ±ν, and its cubic trace cancels exactly; the cross terms between blocks are what remains to control. Test: evaluate this bound on the CUE values of tr G̃² and tr G̃³ for the Montgomery–Taylor window.

Result of the test (20 CUE samples, N = 400, Fejér kernel): Cauchy–Schwarz gives rank/N ≥ 0.893, Hölder gives only 0.758. The safe-direction cubic bound is weaker than the quadratic one, so it brings nothing.

## Step 7 (closed, calculation): why the cubic sign cannot be repaired, and how much of the cubic data is usable

A polynomial certificate q(λ) = aλ + bλ² + cλ³ for counting positive eigenvalues of an indefinite form needs q ≤ 1 for λ > 0 and q ≤ 0 for λ < 0. Letting λ → +∞ forces c ≤ 0, and letting λ → −∞ forces c ≥ 0, so c = 0. A cubic term can enter only when the negative eigenvalues, which come from off-line pairs, are bounded or rare.

How rare they are: a pair 1/2 ± δ + iγ, at depth y = δ log T/(2π) in mean-spacing units, makes p_j grow like e^{2πjy/N}. So the triple term p_a p_b conj p_{a+b} grows like e^{4π(a+b)y/N}. Zero-density estimates N(σ,T) ≪ T^{1−c(σ−1/2)} log T make such pairs a proportion ≪ e^{−2πcy}. The averaged triple term converges when (a+b)/N < c/2. Known values of c:
- Selberg: 1/4 ([Simonič, arXiv:1910.08274](https://arxiv.org/abs/1910.08274));
- Jutila: 1 − ε;
- Conrey: 8/7 − ε (announced).

So with Jutila's estimate, only triples with a + b < N/2 can be used. `lp.py N 1 0.5` measures the ceiling with only those triples.

| N | pair only | a+b < N/2 (Jutila, c = 1) | a+b < 4N/7 (Conrey, c = 8/7) | a+b < 3N/4 (c = 3/2) | all triples |
|---|---|---|---|---|---|
| 8 | 0.685998 | 0.686088 | | | 0.708804 |
| 12 | 0.683750 | 0.683932 | 0.684357 | 0.685026 | 0.704965 |

The usable triples add only 0.0001–0.0002. Almost all of the gain comes from triples with a + b ≥ N/2, and those need a zero-density constant c > 1, beyond Jutila's. Even c = 3/2 gives only +0.0013. The gain sits at the edge a + b ≈ N, where Σ|ξ| ≈ 2, and using it would need a zero-density constant close to 2. The three-point route is therefore closed with known density estimates.
