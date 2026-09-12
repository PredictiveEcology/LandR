// Summed-area (integral image) window counts, for convertUnwantedLCC(method = "nearestRandom").
//
// For each unwanted pixel p and each candidate land-cover class c we need the number of
// cells of class c inside the rectangular window centred on p with half-widths (kx, ky).
// The half-widths vary per pixel -- each is set by that pixel's distance to its nearest
// available class -- so a fixed-radius terra::focal() pass cannot be used, and a naive
// per-pixel window scan would cost O(kx * ky) per pixel. That would reintroduce exactly
// the blow-up this code exists to avoid: in the real Western-Alberta-Upland study area the
// deepest unwanted pixel sits 739 cells from any valid class, i.e. a ~1479 x 1479 window.
//
// A summed-area table reduces each window count to 4 array lookups, independent of window
// size: one O(ncell) prefix-sum pass per candidate class, then O(1) per (pixel, class).
//
// The window is rectangular (Chebyshev in cell units) rather than circular by design. The
// former spread2()-based implementation used the default directions = 8, so "reachable
// within k iterations" was itself a square window; counting over a rectangle therefore
// reproduces that implementation's frequency weighting more faithfully than a disc would.

#include <Rcpp.h>
#include <algorithm>
#include <vector>

using namespace Rcpp;

// Counts of each candidate class within a per-pixel rectangular window.
//
// lccVals     land-cover value of every cell, in terra's row-major cell order
// candClasses candidate classes to count, in the column order of the returned matrix
// nrow, ncol  dimensions of the raster lccVals came from
// cells0      0-based cell indices of the pixels to count around
// kx, ky      per-pixel window half-widths, in cells (same length as cells0)
//
// returns     an integer matrix, length(cells0) rows by length(candClasses) columns
// [[Rcpp::export]]
IntegerMatrix windowCountsByClassCpp(
    IntegerVector lccVals,
    IntegerVector candClasses,
    int nrow,
    int ncol,
    IntegerVector cells0,
    IntegerVector kx,
    IntegerVector ky
) {
  const R_xlen_t nCell = lccVals.size();
  const R_xlen_t nCand = candClasses.size();
  const R_xlen_t nPix = cells0.size();

  if (nrow < 1 || ncol < 1) {
    stop("nrow and ncol must both be >= 1");
  }
  if (nCell != static_cast<R_xlen_t>(nrow) * static_cast<R_xlen_t>(ncol)) {
    stop("length(lccVals) must equal nrow * ncol");
  }
  if (kx.size() != nPix || ky.size() != nPix) {
    stop("kx and ky must be the same length as cells0");
  }
  // NA_INTEGER is INT_MIN, so an NA (or negative) half-width would overflow the window bounds
  for (R_xlen_t i = 0; i < nPix; i++) {
    if (kx[i] == NA_INTEGER || ky[i] == NA_INTEGER || kx[i] < 0 || ky[i] < 0) {
      stop("kx and ky must be non-negative and non-NA");
    }
    if (kx[i] > ncol || ky[i] > nrow) {
      stop("kx and ky must not exceed the raster dimensions");
    }
  }

  IntegerMatrix out(nPix, nCand);

  // S is (nrow + 1) x (ncol + 1) with a zero first row/column, so that
  // S[(r + 1) * W + (c + 1)] is the count over the cell block [0..r] x [0..c].
  const R_xlen_t W = static_cast<R_xlen_t>(ncol) + 1;
  std::vector<int> S(static_cast<R_xlen_t>(nrow + 1) * W, 0);

  for (R_xlen_t j = 0; j < nCand; j++) {
    const int cc = candClasses[j];

    std::fill(S.begin(), S.end(), 0);
    for (int r = 0; r < nrow; r++) {
      const R_xlen_t vOff = static_cast<R_xlen_t>(r) * ncol;
      const R_xlen_t sPrev = static_cast<R_xlen_t>(r) * W;
      const R_xlen_t sCur = sPrev + W;
      for (int c = 0; c < ncol; c++) {
        const int v = lccVals[vOff + c];
        const int hit = (v != NA_INTEGER && v == cc) ? 1 : 0;
        S[sCur + c + 1] = hit + S[sPrev + c + 1] + S[sCur + c] - S[sPrev + c];
      }
    }

    for (R_xlen_t i = 0; i < nPix; i++) {
      const R_xlen_t cell = cells0[i];
      if (cell < 0 || cell >= nCell) {
        stop("cells0 contains an out-of-range cell index");
      }
      const int r = static_cast<int>(cell / ncol);
      const int c = static_cast<int>(cell % ncol);

      const int r0 = std::max(0, r - ky[i]);
      const int c0 = std::max(0, c - kx[i]);
      const int r1 = std::min(nrow - 1, r + ky[i]);
      const int c1 = std::min(ncol - 1, c + kx[i]);

      out(i, j) = S[static_cast<R_xlen_t>(r1 + 1) * W + (c1 + 1)] -
        S[static_cast<R_xlen_t>(r0) * W + (c1 + 1)] -
        S[static_cast<R_xlen_t>(r1 + 1) * W + c0] +
        S[static_cast<R_xlen_t>(r0) * W + c0];
    }
  }

  return out;
}
