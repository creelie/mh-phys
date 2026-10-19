# More than 67.31 percent of the zeta zeros are simple and on the critical line

Deep Bhattacharjee (ORCID [0000-0003-0466-750X](https://orcid.org/0000-0003-0466-750X); itsdeep@live.com, d.bhattacharjee@erl-forschung.de)
and Ushashi Bhattacharya (ORCID [0009-0002-2254-3914](https://orcid.org/0009-0002-2254-3914)).

Archived on Zenodo: DOI to be added after the first release.

This repository contains the paper *More Than 67.31 Percent of the Zeta Zeros
Are Simple and on the Critical Line* (`tex/main.tex`, compiled to `main.pdf`;
REVTeX 4-2) and the programs that check its one computational step.

## Result

Let N(T,2T) count the nontrivial zeros of ζ with T < γ ≤ 2T, and N₀ˢ(T,2T)
those that are simple and lie on the critical line. Then

    liminf N₀ˢ(T,2T)/N(T,2T) ≥ (3125/3112)(H − 7/3200) = 0.6731133552…,

where H = 3/2 − cot(1/√2)/√2 = 0.6725007036… is the Montgomery–Taylor
constant. The previous bounds were 0.6725007 (Alpöge–Furman, arXiv:2608.13637;
Lamzouri, arXiv:2609.02882) and 0.6730085 (creelie/zeta-simple-zeros).

## Status

**Proved by hand in the paper:**

- the rank–trace framework of Alpöge and Furman, as used here (cited);
- the stability lemma and its pinching corollary;
- the eigenvalue lemma: points at mutual distance at least 11/20 have an
  overlap matrix with every eigenvalue at most 1/(δ₁K₀) = 1.97903 < 2, by a
  Fejér majorant; the bound is sharp;
- the block decomposition into runs, pairs and singletons, which lets the
  full defect of long blocks be used without a cap.

**Checked by computer (rigorous interval arithmetic):**

- the eight-zero window inequality F₈ ≥ 13/3125 with μ₀ = 1/3200
  (Proposition 11), by branch and bound over a table of lower bounds for the
  kernel computed in 128-bit ball arithmetic: 10,039,560,518 nodes, verified
  (`verify/bbp8_3200_final.log`).

**Bounded:** no window length lifts the method above 0.67357, for any weight (periodic gap
patterns, `checks/chain_all.py`).

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
| `verify/bbp.c`, `verify/bb.c` | parallel and serial branch-and-bound verifiers (round-toward-minus-infinity arithmetic) |
| `verify/*.log` | logs of the certified runs (μ₀ = 1/3200 and 1/2000) |
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
gcc -O2 -frounding-math -fopenmp -o bbp verify/bbp.c -lm
OMP_NUM_THREADS=4 ./bbp 8 4000 verify/w4000.bin 1 3200 13 3125   # ends with VERIFIED
python3 checks/eig_lemma.py
```

The certificate takes about 15 minutes on two cores. The workflow
`.github/workflows/verify.yml` rebuilds both tables and repeats it on every
push. Longer windows are certified by `.github/workflows/farm.yml`, which
splits the same search into shards across runners, two at a time by default (`bbp … target index count`);
every shard rebuilds the same list of work boxes and prints its checksum, and
the last job checks that all shards verified the same list.

## Credits

The route to the argument was found with the learn-from-failure protocol of
Anaxia, the BLLM engine (creelie/eBLLM, release 3.0.5,
[doi:10.5281/zenodo.23469389](https://doi.org/10.5281/zenodo.23469389)),
described in D. Bhattacharjee and U. Bhattacharya, *BLLM: Beyond the LLM with
Deep Human-Mind-AI*. Deep Bhattacharjee carried out the computations with AI
assistance from Claude (Anthropic).

## Licence

MIT.
