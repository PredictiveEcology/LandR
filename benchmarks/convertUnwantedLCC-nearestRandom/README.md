# `convertUnwantedLCC(method = )`

Supporting evidence for how `convertUnwantedLCC()` allocates replacement land-cover classes:
why the deterministic nearest-class rule briefly present in 1.2.0.9004 was **removed**, and
why both surviving methods weight by local abundance.

## Why

1.2.0.9004 replaced an iterative `spread2()` search — whose run time grew with the square of
the blob radius, and which could run for hours without finishing — with a deterministic
nearest-available allocation. That fixed the blow-up, but it also changed *what gets
imputed*.

The old search sampled among all valid cells within the radius at which it first found one,
so a class was picked in proportion to how much of it was nearby. Taking only the nearest
class gives no weight to abundance, and wherever two or more classes tie at the minimum
distance it must fall back on a tie-break. The 1.2.0.9004 tie-break took the lowest class
code, every time.

Ties are not rare: **35–41%** of unwanted pixels on real landscapes. And because the Canada
LCC codes run non-vegetated → non-forest vegetation → forest, "lowest code wins" is not a
neutral rule. It moved roughly **one in fourteen** unwanted pixels out of forest altogether.
That rule has therefore been removed; it is reconstructed in these scripts only as the
baseline the current methods are scored against.

Both surviving methods draw one of the pixel's available classes with probability
proportional to that class's abundance in the pixel's neighbourhood — the smallest window
reaching its nearest available class, i.e. the window at which `spread2()` would have
stopped. The window is rectangular because `spread2(directions = 8)`'s was too, and counts
come from a summed-area table, so a window 1500 cells across costs the same as one 3 cells
across. They differ only in where the draw comes from:

* **`"nearestWeighted"`** (default) keys it on the pixel's ground position — deterministic,
  no `set.seed()`, stable under `Cache()`, and because the key is the cell *centre* rather
  than the cell *index*, a grid-aligned crop reproduces its parent raster cell for cell.
* **`"nearestRandom"`** draws from the RNG, for when replicates should genuinely differ.

## The tie-break, measured

`05_bias_diagnostics.R`. "Share given the lowest tied class" is 100% by construction for the
removed rule; the old algorithm's value is the target. "Adjacency" is the share of
rook-adjacent unwanted-pixel pairs assigned the same class — a directional or patch artifact
would push it up.

| landscape | pixels with a tie | lowest-tied share — spiral (old) | lowest-code (REMOVED) | nearestWeighted | nearestRandom |
|---|---:|---:|---:|---:|---:|
| medium  | 34.9% | 45.0% | **100.0%** | 48.3% | 46.1% |
| bigblob | 40.6% | 52.7% | **100.0%** | 47.7% | 47.5% |

| landscape | adjacency — spiral | lowest-code | nearestWeighted | nearestRandom |
|---|---:|---:|---:|---:|
| medium  | 72.6% | 74.2% | 70.1% | 71.0% |
| bigblob | 69.3% | 73.6% | 65.8% | 66.3% |

The bias was never spatial — adjacency barely moves — it was entirely in the class-code
tie-break. Both current methods resolve ties at close to the old algorithm's rate.

![bias diagnostics](fig2_bias_diagnostics.png)

## What the removed rule cost, in cover type

Pooled over four real landscapes, weighted by unwanted pixels (`06_class_bias_by_cover_type.R`).
Values are % of all unwanted pixels; labels from LandR's own crosswalk in
`prepInputs_NTEMS_LCC_FAO()`.

| class | cover type | spiral (old) | lowest-code (REMOVED) | nearestWeighted | nearestRandom |
|---:|---|---:|---:|---:|---:|
|  40 | bryoids    |  0.10 |  0.16 (1.61×) |  0.06 (0.60×) |  0.09 (0.97×) |
|  50 | shrubs     |  9.63 | **16.25 (1.69×)** |  9.77 (1.01×) |  9.20 (0.96×) |
| 100 | herbs      |  1.60 |  2.20 (1.37×) |  1.46 (0.91×) |  1.57 (0.98×) |
| 210 | coniferous | 66.59 | 70.70 (1.06×) | 67.34 (1.01×) | 67.00 (1.01×) |
| 220 | broadleaf  | 15.12 | **8.75 (0.58×)** | 14.93 (0.99×) | 15.21 (1.01×) |
| 230 | mixedwood  |  6.97 | **1.96 (0.28×)** |  6.45 (0.92×) |  6.93 (0.99×) |

| cover group | spiral (old) | lowest-code (REMOVED) | nearestWeighted | nearestRandom |
|---|---:|---:|---:|---:|
| non-forest vegetation | 11.33 | **18.61 (1.64×)** | 11.29 (**1.00×**) | 10.86 (0.96×) |
| forest                | 88.68 | **81.41 (0.92×)** | 88.72 (**1.00×**) | 89.14 (1.01×) |

The removed rule inflated shrubs by 69% and thinned broadleaf by 42% and mixedwood by 72%.
That is not cosmetic for a succession model: a pixel imputed as shrubs or herbs carries no
tree cohorts at all, and broadleaf/mixedwood → coniferous shifts the deciduous fraction that
drives `partitionBiomass()` and the fire regime. Both current methods land on the old
algorithm's forest/non-forest split to within 0–4%.

![class bias](fig3_class_bias.png)

## Were there better deterministic rules?

`07_tiebreak_candidates.R` scores the alternatives that keep determinism, by total-variation
distance from the old algorithm's class mix (lower is closer; `floor` is the old algorithm
against itself under different seeds):

| landscape | floor (old vs old) | lowest-code | modal among tied | modal overall | **nearestWeighted** | nearestRandom |
|---|---:|---:|---:|---:|---:|---:|
| small   | 0.0486 | 0.2953 | 0.0472 | 0.0472 | 0.0669 | 0.0157 |
| medium  | 0.0072 | 0.1155 | 0.0352 | 0.0497 | 0.0148 | 0.0098 |
| large   | 0.0043 | 0.1008 | 0.0407 | 0.0532 | 0.0099 | 0.0042 |
| bigblob | 0.0207 | 0.1456 | 0.0640 | 0.0693 | 0.0143 | 0.0041 |

Breaking ties by whichever tied class is locally most abundant (`modal among tied`), or
taking the locally dominant class outright (`modal overall`), is a large improvement on
lowest-code but still 3–9× above the noise floor: winner-take-all over-concentrates the
dominant class. Keeping the full weighting and making only the *draw* deterministic is what
lands at the floor.

Note when reading that table: `nearestWeighted` is a single realization, so it carries one
draw's worth of noise, whereas the `nearestRandom` column is averaged over three seeds. The
fair comparison for `nearestWeighted` is the floor column, not `nearestRandom`.

## Real landscapes: composition and cost

`03_method_comparison.R`. Per-pixel agreement **cannot** discriminate here — the old
algorithm does not agree with itself (53.9–79.2%), so that column is saturated by its own
tie-break noise. Composition does.

| landscape | ncell | unwanted | spiral (s) | nearestWeighted (s) | nearestRandom (s) | TVD: floor | lowest-code | nearestWeighted | nearestRandom |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| small   |    43,681 |   254 | 0.04 | 0.12 | 0.11 | 0.0486 | 0.2953 | 0.0669 | 0.0157 |
| medium  |   341,056 | 2,889 | 0.11 | 0.40 | 0.39 | 0.0072 | 0.1155 | 0.0148 | 0.0098 |
| large   | 1,175,056 | 6,770 | 0.41 | 1.47 | 1.41 | 0.0043 | 0.1008 | 0.0099 | 0.0042 |
| bigblob |   250,000 | 1,145 | 0.07 | 0.27 | 0.27 | 0.0207 | 0.1456 | 0.0143 | 0.0041 |

On these landscapes the unwanted pixels are *scattered*, so the old loop is already fast —
its cost is driven by blob radius, not raster size.

## Cost: no return of the blow-up

`04_scaling_all_methods.R` — self-contained, no private data: a single unwanted blob of
increasing radius.

| blob radius (cells) | ncell | unwanted | spiral (old) | nearest-allocation | +weighted draw |
|---:|---:|---:|---:|---:|---:|
|  10 |     676 |    316 |   0.12 s | 0.04 s | 0.05 s |
|  20 |   2,704 |  1,264 |   1.23 s | 0.02 s | 0.04 s |
|  40 |  10,816 |  5,024 |  26.89 s | 0.02 s | 0.04 s |
|  80 |  43,264 | 20,108 | **DNF** (>120 s cap; 9,980 unresolved) | 0.06 s | 0.09 s |
| 160 | 173,056 | 80,452 | **DNF** (>120 s cap; 68,052 unresolved) | 0.18 s | 0.25 s |

The weighted draw adds a second set of distance transforms (to size the window in cells,
independent of CRS and resolution) and one summed-area table per candidate class — all
`O(ncell)` and all independent of blob depth. It stays flat exactly where the old
implementation diverges.

![scaling](fig1_scaling_all_methods.png)

## Tests

`tests/testthat/test-cohorts.R`:

* `"...replaces unwanted classes from the neighbourhood"` — return schema, every unwanted
  pixel resolved, no `classesToReplace` emitted, per-ecoregion availability respected; both
  methods.
* `"...'nearestWeighted' is deterministic and crop-stable"` — identical across repeated
  calls, unaffected by the RNG state, and a grid-aligned crop reproduces the parent raster
  cell for cell.
* `"...draws are weighted by local abundance"` — one unwanted pixel whose 8 neighbours are
  one 210 (also the strictly nearest) and seven 220; must return 210 about 1 time in 8,
  which neither a nearest-only nor a uniform rule would do.
* `"...shows no bias toward low class codes"` — guards the removed defect directly.
* `"...methods agree when only one class can be chosen"`.

## Reproducing

The four real landscapes are private LandWeb data; the scaling sweep is self-contained.

```sh
LANDR_SRC=. Rscript benchmarks/convertUnwantedLCC-nearestRandom/04_scaling_all_methods.R
LCC_BENCH_DIR=<dir with v2_input_*.tif> LANDR_SRC=. \
  Rscript benchmarks/convertUnwantedLCC-nearestRandom/03_method_comparison.R
LANDR_SRC=. Rscript benchmarks/convertUnwantedLCC-nearestRandom/06_class_bias_by_cover_type.R
LCC_BENCH_DIR=<dir with v2_input_*.tif> LANDR_SRC=. \
  Rscript benchmarks/convertUnwantedLCC-nearestRandom/05_bias_diagnostics.R
LCC_BENCH_DIR=<dir with v2_input_*.tif> LANDR_SRC=. \
  Rscript benchmarks/convertUnwantedLCC-nearestRandom/07_tiebreak_candidates.R
```
