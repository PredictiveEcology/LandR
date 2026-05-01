# LANDISDisp benchmarks: R vs Rcpp

Generated: 2026-05-01 10:57:02.422845
Host: A159568/x86_64 | R: R version 4.5.2 (2025-10-31)

All timings include the call into LANDISDisp() (spiral + ward-prob prep + the
loop). The first call per size is excluded (warm-up); reported numbers are the
median over the listed reps. Both implementations use the same RNG stream, so
the output rows match exactly (verified by seed-locked tests).

| size | pgm cells | rcv rows | src rows | reps | R median (s) | Cpp median (s) | speedup | output rows | identical |
|------|-----------|---------:|---------:|-----:|------------:|--------------:|--------:|-----------:|:---------:|
| tiny | 400 | 13 | 8 | 5 | 0.044 | 0.013 | 3.38x | 48 | yes |
| small | 2500 | 41 | 41 | 5 | 0.067 | 0.015 | 4.64x | 448 | yes |
| medium | 14400 | 155 | 147 | 3 | 0.149 | 0.024 | 6.26x | 2877 | yes |
| large | 62500 | 432 | 462 | 2 | 0.428 | 0.149 | 2.88x | 14912 | yes |
| xlarge | 250000 | 1137 | 1187 | 1 | 1.552 | 0.485 | 3.20x | 47439 | yes |
| xlarge_dense | 250000 | 3302 | 753 | 1 | 2.912 | 0.706 | 4.13x | 43938 | yes |
| xxlarge | 640000 | 5927 | 2130 | 1 | 6.476 | 1.029 | 6.29x | 130761 | yes |
| xxxlarge | 9000000 | 1348 | 29060 | 1 | 6.227 | 1.532 | 4.06x | 171576 | yes |
