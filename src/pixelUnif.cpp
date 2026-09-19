// Deterministic, seed-free uniform draws keyed on a pixel's ground position, for
// convertUnwantedLCC(method = "nearestWeighted").
//
// "nearestWeighted" and "nearestRandom" allocate identically -- each unwanted pixel takes
// one of its available classes with probability proportional to that class's abundance in
// the pixel's neighbourhood. They differ only in where the uniform comes from:
// "nearestRandom" calls runif() and so varies with the RNG seed, while "nearestWeighted"
// derives it here from the pixel itself. That makes the result reproducible without
// set.seed(), and stable under Cache(), which does not key on RNG state.
//
// Keyed on COORDINATES, not on cell index. A cell index is only meaningful relative to one
// raster's extent, so indexing would give a different answer on a cropped subset than on
// the full raster -- exactly the workflow (develop on a small crop, then scale up) where
// you most want the two to agree. Cell-centre coordinates are a property of the ground, so
// a grid-aligned crop reproduces the parent raster's allocation cell for cell.
//
// The coordinates are quantized to the raster grid before hashing -- ix = round(x / resx),
// iy = round(y / resy) -- so the key is an exact integer pair rather than a double whose
// low bits might differ between a crop and its parent. The quantization is anchored at the
// CRS origin, not at the raster's own extent, so any two rasters sharing a CRS, resolution
// and grid alignment agree. Changing the resolution or reprojecting changes the pixels
// themselves, and is expected to change the result.
//
// The mix is splitmix64's finalizer, applied to each axis and then to the combination.
// Inputs here are near-sequential integers, so a weak hash (a single linear congruential
// step, or R's bitwShiftL, which overflows past 2^31) would leave visible structure across
// neighbouring cells and bias the allocation; adjacent cells must land far apart.

#include <Rcpp.h>
#include <cmath>
#include <cstdint>

using namespace Rcpp;

static inline uint64_t splitmix64(uint64_t x) {
  x += 0x9E3779B97F4A7C15ULL;
  x = (x ^ (x >> 30)) * 0xBF58476D1CE4E5B9ULL;
  x = (x ^ (x >> 27)) * 0x94D049BB133111EBULL;
  return x ^ (x >> 31);
}

// Uniform in [0, 1) for each cell centre, deterministic and independent of the RNG state.
//
// x, y        cell-centre coordinates, in the raster's CRS
// resx, resy  raster resolution, used to quantize those coordinates onto the grid
// [[Rcpp::export]]
NumericVector pixelUnifCpp(NumericVector x, NumericVector y, double resx, double resy) {
  const R_xlen_t n = x.size();
  if (y.size() != n) {
    stop("x and y must be the same length");
  }
  if (!(resx > 0) || !(resy > 0)) {
    stop("resx and resy must be positive");
  }

  NumericVector out(n);
  for (R_xlen_t i = 0; i < n; i++) {
    if (NumericVector::is_na(x[i]) || NumericVector::is_na(y[i])) {
      stop("x and y must not be NA");
    }
    // two's-complement cast keeps negative coordinates (west/south of the origin) exact
    const int64_t ix = static_cast<int64_t>(std::llround(x[i] / resx));
    const int64_t iy = static_cast<int64_t>(std::llround(y[i] / resy));
    const uint64_t h = splitmix64(
      splitmix64(static_cast<uint64_t>(ix)) ^ (static_cast<uint64_t>(iy) + 0x9E3779B97F4A7C15ULL)
    );
    // top 53 bits -> [0, 1), matching the precision of a double
    out[i] = static_cast<double>(h >> 11) / 9007199254740992.0;
  }
  return out;
}
