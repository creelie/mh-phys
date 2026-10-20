# More than 67.32 percent of the zeta zeros are simple and on the critical line

Deep Bhattacharjee (ORCID [0000-0003-0466-750X](https://orcid.org/0000-0003-0466-750X); itsdeep@live.com, d.bhattacharjee@erl-forschung.de)
and Ushashi Bhattacharya (ORCID [0009-0002-2254-3914](https://orcid.org/0009-0002-2254-3914)).

Archived on Zenodo: release 1.2.0 [doi:10.5281/zenodo.23484951](https://doi.org/10.5281/zenodo.23484951) (bound 0.6732133974); all versions [doi:10.5281/zenodo.23480714](https://doi.org/10.5281/zenodo.23480714); release 1.1.0 [doi:10.5281/zenodo.23483313](https://doi.org/10.5281/zenodo.23483313) (bound 0.6731985970); release 1.0.0 [doi:10.5281/zenodo.23480715](https://doi.org/10.5281/zenodo.23480715) (bound 0.6731133552).

This repository contains the paper *More Than 67.32 Percent of the Zeta Zeros
Are Simple and on the Critical Line* (`tex/main.tex`, compiled to `main.pdf`;
REVTeX 4-2) and the programs that check its one computational step.

## Result

Let N(T,2T) count the nontrivial zeros of ζ with T < γ ≤ 2T, and N₀ˢ(T,2T)
those that are simple and lie on the critical line. Then

    liminf N₀ˢ(T,2T)/N(T,2T) ≥ (125000/124111)(H − 7999/2000000) = 0.6732896395…,

where H = 3/2 − cot(1/√2)/√2 = 0.6725007036… is the Montgomery–Taylor
constant. The previous bounds were 0.6725007 (Alpöge–Furman, arXiv:2608.13637;
Lamzouri, arXiv:2609.02882), 0.6730085 (creelie/zeta-simple-zeros),
0.6731134 (release 1.0.0 of this repository, uniform weights on eight zeros),
0.6731986 (release 1.1.0, weighted eight-zero windows) and 0.6732134
(release 1.2.0, weighted nine-zero windows).

## Status

**Proved by hand in the paper:**

- the rank–trace framework of Alpöge and Furman, as used here (cited);
- the stability lemma and its pinching corollary;
- the eigenvalue lemma: points at mutual distance at least 11/20 have an
  overlap matrix with every eigenvalue at most 1/(δ₁K₀) = 1.97903 < 2, by a
  Fejér majorant; the bound is sharp;
- the block decomposition into runs, pairs and singletons, which lets the
  full defect of long blocks be used without a cap;
- that the window inequality may use any weights: each pair of the ten
  zeros and each gap gets its own weight, subject only to the pair weights at
  each distance adding up to at most 2 (Lemma 9, Proposition 10).

**Checked by computer (rigorous interval arithmetic):**

- the ten-zero window inequality F ≥ 889/125000 = 0.007112 for the weights of
  Table I (`verify/weights10.txt`, linear total 7999/2000000; Proposition 11),
  by branch and bound with two lower bounds on a box: one from a table of
  lower bounds for the kernel on cells, and a second-order Taylor bound from
  tables of w and w′ at grid points and of lower bounds for w″ on cells, all
  computed in 128-bit ball arithmetic: 1,424,333,714 nodes in 2 minutes 48
  seconds on four threads (`verify/bbpcf10.log`);
- the nine-zero inequality F ≥ 7/1000 for the weights of release 1.2.0
  (`verify/weights9.txt`, linear total 19999/5000000), which gives 0.6732134:
  52,891,178,741 nodes (`verify/bbpww9.log`);
- the eight-zero inequality F ≥ 19/4000 for the weights of release 1.1.0
  (`verify/weights8.txt`, linear total 12499/5000000), which gives
  0.6731986: 21,067,790,614 nodes (`verify/bbpww8.log`);
- the same inequality with uniform weights, F₈ ≥ 13/3125 for μ₀ = 1/3200,
  which gives the bound 0.6731133 of release 1.0.0: 10,039,560,518 nodes
  (`verify/bbp8_3200_final.log`).

**Chosen by computer (no proof rests on it):** the weights, by a cutting-plane
linear programme (`checks/weights_lp.py`, `weights_lp10.txt`); the
programme reaches 0.6732900 at ten zeros, and the restricted chain energy
shows that no weights for ten zeros give more than 0.6732908
(`checks/chain_window.py`, `chain_window9.txt`, `chain_window10.txt`;
0.6732646 for nine zeros, 0.6733072 for eleven, 0.6733252 for twelve). The
minimum of F for the rounded weights is 0.00711257
(`checks/weights_minima10.txt`), so the certified 0.007112 is within
6·10⁻⁷ of it; the target 0.0071125 fails at a single cell, so cells of width
1/4000 cannot go further with these weights.

**Bounded:** no window length or choice of weights lifts the method above
0.67357 (periodic gap patterns, `checks/chain_all.py`).

**Not claimed:** this is not a proof of the Riemann hypothesis, and it does not
bear on it. A bound of this kind is insensitive to a set of zeros of density
zero. It rests on the framework of an arXiv preprint (Alpöge–Furman) that is
not yet refereed.

## Layout

| Path | Contents |
|---|---|
| `tex/` | LaTeX source, figures, figure generator `tex/figures/src/make_figures.py` |
| `main.pdf`, `preprintSimpleZerosCriticalLine/` | compiled paper |
| `MasterSimpleZerosCriticalLine_masterzipfinal2.zip` | everything in one archive: this repository, the paper and the Hexagon submission package |
| `verify/table.py` | lower bounds for w = k² on cells of width 1/4000 (python-flint, Arb, 128 bits); `w4000.bin` (46,000 cells, sha256 `8ea9e020…afa42d5`) and `w4000_80k.bin` (80,000 cells, same first 46,000) |
| `verify/dtable.py` | second-order tables: lower bounds for w and enclosures of w′ at the grid points m/4000, lower bounds for w″ on cells (Arb, 128 bits); `d4000_80k.bin` (sha256 `8915be68…baca14c105`) |
| `verify/bbp_cf.c`, `verify/weights10.txt` | second-order weighted branch-and-bound verifier and the weights of Table I |
| `verify/bbp_ww.c` | weighted verifier with the cell bound only (`weights9.txt`: release 1.2.0; `weights8.txt`: release 1.1.0) |
| `verify/bbp.c`, `verify/bb.c` | parallel and serial verifiers for uniform weights (round-toward-minus-infinity arithmetic) |
| `verify/*.log` | logs of the certified runs (`bbpcf10.log`: the weights of Table I, ten zeros; `bbpww9.log`: nine zeros, release 1.2.0; `bbpww8.log`: eight zeros, release 1.1.0; uniform μ₀ = 1/3200 and 1/2000) |
| `checks/weights_lp.py`, `round_weights.py`, `weights_minima.py` | linear programme for the weights, rounding, local minima of the weighted F (outputs `weights_lp10.txt`, `ww_R10_mu0.004.npz`, `weights_minima10.txt`; the same for R = 8, 9) |
| `checks/chain_window.py`, `chain_window8.txt`, `chain_window9.txt`, `chain_window10.txt` | limit for a window of R points (pairs at most R − 1 apart) |
| `checks/eig_lemma.py` | Arb check of the majorant inequality and of the eigenvalue bound; lattice eigenvalues |
| `checks/localmin.py`, `lm*.txt` | local minima of F_R for R = 7, 8, 9 |
| `checks/chain_all.py`, `chain_all.txt` | upper limit of the method over every weight |
| `checks/tilt.py`, `window_H.py` | tilted windows |
| `ceiling/` | periodic linear programmes for the ceiling of pair and three-point data |
| `notes/` | working notes and the failure log of the search |

## Reproduce

```
python3 -m pip install python-flint==0.9.0 numpy scipy
python3 verify/table.py 4000 46000 w4000_check.bin   # prints the sha256 above
python3 verify/table.py 4000 80000 w80k_check.bin   # equals verify/w4000_80k.bin
python3 verify/dtable.py 4000 80000 4 d80k_check.bin   # equals verify/d4000_80k.bin
gcc -O2 -frounding-math -fopenmp -o bbp_cf verify/bbp_cf.c -lm
OMP_NUM_THREADS=4 ./bbp_cf 10 4000 verify/w4000_80k.bin verify/d4000_80k.bin verify/weights10.txt 889 125000   # ends with VERIFIED R=10
gcc -O2 -frounding-math -fopenmp -o bbp_ww verify/bbp_ww.c -lm
OMP_NUM_THREADS=4 ./bbp_ww 9 4000 verify/w4000_80k.bin verify/weights9.txt 7 1000 2000000   # release 1.2.0, ends with VERIFIED R=9
OMP_NUM_THREADS=4 ./bbp_ww 8 4000 verify/w4000_80k.bin verify/weights8.txt 19 4000 400000   # release 1.1.0
gcc -O2 -frounding-math -fopenmp -o bbp verify/bbp.c -lm
OMP_NUM_THREADS=4 ./bbp 8 4000 verify/w4000.bin 1 3200 13 3125   # uniform weights, ends with VERIFIED
python3 checks/eig_lemma.py
```

The ten-zero certificate takes under three minutes on four threads; the
nine-zero one, with the cell bound only, 30 minutes, the eight-zero weighted
one 26 minutes and the uniform one about 15 minutes on two cores. The
second-order bound is what makes the ten-zero run short: with the cell bound
alone and the weaker target 22/3125, one of 64 shards of the search took
88,467,547,937 nodes and 57 minutes, against 33,813,175 nodes and 8 seconds
with both bounds. The workflow `.github/workflows/verify.yml` rebuilds all
three tables and repeats the ten-zero, nine-zero and uniform certificates on
every push. Longer windows are certified by `.github/workflows/farm.yml`, which
splits the same search into shards across runners, two at a time by default (`bbp … target index count`);
every shard rebuilds the same list of work boxes and prints its checksum, and
the last job checks that all shards verified the same list.

## Credits

The route to the argument was found with the learn-from-failure protocol of
Anaxia, the BLLM engine (creelie/eBLLM, release 3.1.1,
[doi:10.5281/zenodo.23484395](https://doi.org/10.5281/zenodo.23484395)),
described in D. Bhattacharjee and U. Bhattacharya, *BLLM: Beyond the LLM with
Deep Human-Mind-AI*. Deep Bhattacharjee carried out the computations with AI
assistance from Claude (Anthropic).

## Licence

MIT.
