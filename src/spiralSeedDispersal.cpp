// LANDIS-II spiral seed dispersal — Rcpp port of the inner loop of
// spiralSeedDispersalR() in R/seedDispersalLANDIS.R.
//
// The R wrapper does the lightweight preprocessing (spiral, wardProb-by-(dist,
// species), rcvFull) and hands fully-prepared arrays to this function. This
// function performs the spiral walk and the per-iteration uniform draws.
// The RNG calls below use R::unif_rand() directly, drawing the same number
// of uniforms in the same order as the R reference (which calls
// SpaDES.tools::runifC(sumHasSp) once per iteration), so seeded outputs are
// bit-identical to the R implementation.
//
// Source-presence representation
// -----------------------------
// The R reference materialises a (ncell × numSp) integer matrix where each
// (cell, species) entry is the species code (or NA). At landscape scale this
// is the dominant memory and prep-time cost (e.g., 9M × 7 = 63M ints, plus a
// data.table cartesian join to populate it). We avoid it entirely by using a
// per-pixelGroup species bitmask:
//   - pgv[cell]  = pixelGroup ID at that cell (NA if masked)
//   - srcPgBitmask[pg] = uint64 bitmask, bit (sp-1) set iff (pg, sp) ∈ dtSrc
// "Cell c has species s" reduces to a 1-cycle lookup:
//   const int pg = pgv[c];
//   const bool has = ((srcPgBitmask[pg] >> (sp - 1)) & 1ULL) != 0;
// Memory drops from ncell × numSp ints to maxPg uint64s. The R-side wrapper
// no longer constructs the matrix; the bitmask is built here from dtSrc.

#include <Rcpp.h>
#include <cstdint>

using namespace Rcpp;

// Inputs (built in the R wrapper):
//
//   pixelIndex_in : 1-based pixel index of each receiver row (length = numRcv)
//   speciesCode_in: 1-based species code of each receiver row (length = numRcv)
//   rowOrig_in    : 1-based row of each receiver (length = numRcv)
//   colOrig_in    : 1-based col of each receiver (length = numRcv)
//   seeddist_max_perRow : per-row max seed dispersal distance (numRcv)
//   spiralRow     : row offsets in spiral order (length = numSpiral)
//   spiralCol     : col offsets in spiral order
//   spiralCurDist : per-spiral-step distance in raw units
//                   (= spiral[, "dists"] * cellSize)
//   pgmRows, pgmCols : pixelGroupMap dimensions (terra cell numbering: cell =
//                      (row-1)*pgmCols + col, all 1-based)
//   pgv           : IntegerVector of length ncell with pixelGroup ID at each
//                   cell (NA_INTEGER for masked cells)
//   srcPg         : pixelGroup IDs from dtSrc (length nSrc)
//   srcSpeciesCode: species codes from dtSrc (length nSrc, 1-based)
//   numSp         : total species count (must be ≤ 64 for the bitmask)
//   wardProbByDist: numUniqueDists * numSp matrix; wardProbByDist(distIdx, sp - 1)
//                   = ward prob (already raised to successionTimestep upstream
//                   if needed)
//   activeSpMaxDist : per-species seeddistance_maxMinCellSize
//                     (length numSp + 1, 1-indexed by species code)
//   activeSpMax     : per-species seeddistance_max (same indexing)
//   verbose         : 0/1 — does not change RNG path; only controls whether
//                    DistOfSuccess is recorded (R reference does the same)
//   successionTimestep : already applied to wardProbByDist when > 1, but kept
//                    here for parity in case caller hasn't applied it
//   wardAlreadyExp  : true if caller has already raised probabilities to the
//                    successionTimestep power (avoids re-doing it).
//
// Returns a list with:
//   Success        : LogicalVector length numRcv
//   DistOfSuccess  : NumericVector length numRcv (NA_real_ if no success)

// [[Rcpp::export]]
List spiralLoopCpp(IntegerVector pixelIndex_in,
                   IntegerVector speciesCode_in,
                   IntegerVector rowOrig_in,
                   IntegerVector colOrig_in,
                   IntegerVector seeddist_max_perRow,
                   IntegerVector spiralRow,
                   IntegerVector spiralCol,
                   NumericVector spiralCurDist,
                   int pgmRows,
                   int pgmCols,
                   IntegerVector pgv,
                   IntegerVector srcPg,
                   IntegerVector srcSpeciesCode,
                   int numSp,
                   NumericMatrix wardProbByDist,
                   NumericVector activeSpMaxDist,
                   NumericVector activeSpMax,
                   double cellSize,
                   int successionTimestep,
                   int verbose,
                   bool wardAlreadyExp,
                   bool debug = false) {

  const int numRcv = pixelIndex_in.size();
  const int numSpiral = spiralRow.size();
  const int ncell = pgv.size();

  if (numSp > 64) {
    stop("spiralLoopCpp: numSp must be <= 64 (bitmask width). Got %d.", numSp);
  }

  // Build per-pixelGroup species bitmask. maxPg is the largest pg encountered
  // in pgv (or in srcPg if larger, defensively).
  int maxPg = 0;
  for (int c = 0; c < ncell; ++c) {
    const int pg = pgv[c];
    if (pg != NA_INTEGER && pg > maxPg) maxPg = pg;
  }
  const int nSrcRows = srcPg.size();
  for (int s = 0; s < nSrcRows; ++s) {
    const int pg = srcPg[s];
    if (pg != NA_INTEGER && pg > maxPg) maxPg = pg;
  }
  std::vector<std::uint64_t> srcPgBitmask((size_t) maxPg + 1, 0ULL);
  for (int s = 0; s < nSrcRows; ++s) {
    const int pg = srcPg[s];
    const int sp = srcSpeciesCode[s];
    if (pg == NA_INTEGER || sp == NA_INTEGER) continue;
    if (pg < 0 || pg > maxPg) continue;
    if (sp < 1 || sp > numSp) continue;
    srcPgBitmask[(size_t) pg] |= (1ULL << (sp - 1));
  }

  // Output vectors — initialized to NA / FALSE
  LogicalVector Success(numRcv, false);
  NumericVector DistOfSuccess(numRcv, NA_REAL);

  // Active receiver tracking. Parallel arrays mirror the R code's
  // `activeFullIndex` / `rowOrig` / `colOrig` / `speciesCode`. We additionally
  // maintain `active` — an in-place compacted list of receiver indices still
  // in play. This avoids paying O(numRcv) per spiral step on dropped receivers
  // (which on landscape-scale runs is ~70% of the inner-loop cost). `dropped`
  // is the marker used during a single spiral step to flag freshly-succeeded
  // receivers; we compact at end of step when any drops happened.
  //
  // (We also tried an AoS struct{row,col,sp,fullIdx} layout to improve cache
  // locality. It made the inner loop slightly tighter but the corresponding
  // 4× larger compaction copies cost more than that saved on every workload
  // we measured. SoA + an int[] active list is the faster shape here.)
  std::vector<int>  rowAct(numRcv);
  std::vector<int>  colAct(numRcv);
  std::vector<int>  spAct(numRcv);
  std::vector<int>  fullIdx(numRcv); // index back into the input rows
  std::vector<char> dropped(numRcv, 0);
  std::vector<int>  active;
  active.reserve(numRcv);

  for (int j = 0; j < numRcv; ++j) {
    rowAct[j]  = rowOrig_in[j];
    colAct[j]  = colOrig_in[j];
    spAct[j]   = speciesCode_in[j];
    fullIdx[j] = j;
    active.push_back(j);
  }

  // Active species set: indexed by species code (1..numSp). Element true
  // means that species is still in play (some receiver could still receive
  // seed before exceeding maxMinCellSize).
  std::vector<char> activeSp(numSp + 1, 1);
  // Species code 0 is unused.
  // Mark species with NA / 0 max as dropped immediately. (R version uses NA
  // for "no species present" sentinel; here 0 is used to mark unused slots.)

  double prevCurDist = spiralCurDist[0];
  bool   newCurDist = true;
  int    uniqueDistCounter = 0;
  double lastWardMaxProb = 1.0;

  // Number of currently-live receivers (for early termination)
  int numLive = numRcv;

  GetRNGstate();

  // Workspace arrays reused across iterations to avoid reallocation
  std::vector<int>    hasSpRows;       // receiver-row indices (j) with hasSp
  std::vector<double> wardForRow;      // ward prob lookup for this iteration

  for (int i = 0; i < numSpiral; ++i) {
    const double curDist = spiralCurDist[i];
    const int    sRow    = spiralRow[i];
    const int    sCol    = spiralCol[i];

    // Detect distance step transitions
    if (i > 0) {
      if (curDist > prevCurDist) {
        newCurDist = true;
        uniqueDistCounter += 1;
        prevCurDist = curDist;

        // Drop entire species whose seeddistance_maxMinCellSize is now
        // exceeded by curDist
        bool anySpeciesDropped = false;
        for (int sp = 1; sp <= numSp; ++sp) {
          if (activeSp[sp] && curDist > activeSpMaxDist[sp]) {
            activeSp[sp] = 0;
            anySpeciesDropped = true;
          }
        }

        // Compact `active` to drop receivers whose species was just dropped
        // OR whose per-row pmax(cellSize, seeddistance_max) < curDist.
        // Single-pass in-place compaction preserves order.
        // Skipped if no species transition (per-row max equals species max
        // here, so no row drop is possible without a species drop).
        if (anySpeciesDropped) {
          size_t w = 0;
          for (size_t r = 0; r < active.size(); ++r) {
            const int j = active[r];
            const int sp = spAct[j];
            if (!activeSp[sp]) {
              dropped[j] = 1;
              continue;
            }
            double rowMax = (double) seeddist_max_perRow[fullIdx[j]];
            if (rowMax < cellSize) rowMax = cellSize;
            if (curDist > rowMax) {
              dropped[j] = 1;
              continue;
            }
            active[w++] = j;
          }
          numLive -= (int) (active.size() - w);
          active.resize(w);
        }
      } else {
        newCurDist = false;
      }
    }

    // Lookup ward probabilities for this distance step. Mirrors R reference
    // which recomputes wardProbActual unconditionally inside the newCurDist
    // branch BEFORE any sumHasSp / pre-screen short-circuit. Updating later
    // (after the early-out) would leave wardForRow stale across distance
    // boundaries where a "no-draws-pass" iteration is sandwiched between
    // valid iterations.
    if (newCurDist) {
      wardForRow.assign(numSp + 1, 0.0);
      for (int sp = 1; sp <= numSp; ++sp) {
        double p = wardProbByDist(uniqueDistCounter, sp - 1);
        if (!wardAlreadyExp && successionTimestep > 1) {
          p = 1.0 - std::pow(1.0 - p, (double) successionTimestep);
        }
        if (p > 1.0) p = 1.0;
        wardForRow[sp] = p;
      }
    }

    // Iteration row/col offsets and per-receiver source lookup. Iterate the
    // compacted active list (so dropped & species-dropped receivers cost zero
    // per step). Order is preserved relative to the original receiver index
    // — important for RNG-stream parity with the R reference.
    hasSpRows.clear();

    for (size_t idx = 0; idx < active.size(); ++idx) {
      const int j = active[idx];
      const int sp = spAct[j];

      const int newRow = rowAct[j] + sRow;
      const int newCol = colAct[j] + sCol;
      if (newRow < 1 || newRow > pgmRows ||
          newCol < 1 || newCol > pgmCols) {
        continue; // out of bounds → no source
      }
      // Cell index using terra's row-major numbering, 0-based
      const int cell0 = (newRow - 1) * pgmCols + (newCol - 1);
      const int pg = pgv[cell0];
      if (pg == NA_INTEGER) continue;          // masked cell
      if (pg < 0 || pg > maxPg) continue;      // pg outside dtSrc range
      const std::uint64_t mask = srcPgBitmask[(size_t) pg];
      if (((mask >> (sp - 1)) & 1ULL) == 0ULL) continue;

      hasSpRows.push_back(j);
    }

    const int sumHasSp = (int) hasSpRows.size();
    if (sumHasSp == 0) {
      if (numLive == 0) break;
      continue;
    }
    if (debug) {
      Rprintf("[cpp] i=%d curDist=%.10f n=%d lastMax=%.17f udc=%d\n",
              i + 1, curDist, sumHasSp, lastWardMaxProb, uniqueDistCounter);
    }

    // Draw sumHasSp uniforms — same count, same order as runifC(sumHasSp)
    // in the R reference.
    // Allocate locally to avoid mixing draws between iterations.
    // Doing this in a single tight loop matches the R RNG advance pattern.
    std::vector<double> ran(sumHasSp);
    for (int k = 0; k < sumHasSp; ++k) {
      ran[k] = ::unif_rand();
    }

    // whRanLTprevMaxProb <- which(ran <= lastWardMaxProb)
    // (early-out optimisation from the R version)
    // First check if anything passes the screen.
    bool anyPass = false;
    for (int k = 0; k < sumHasSp; ++k) {
      if (ran[k] <= lastWardMaxProb) { anyPass = true; break; }
    }
    if (!anyPass) {
      // No draw could possibly succeed under any current ward prob.
      if (numLive == 0) break;
      continue;
    }

    // Decide successes. Two paths in R:
    //   if (i == 1) oo <- seq.int(length(ran))   # all hasSpRows succeed
    //   else        oo <- ran[whRan...] < wardRes per species
    // The R code unconditionally consumes `sumHasSp` draws; we already did.
    double newMax = 0.0;
    bool   anySuccess = false;
    if (i == 0) {
      // Self pixel: every receiver with a source on its own pixel succeeds
      for (int k = 0; k < sumHasSp; ++k) {
        const int j = hasSpRows[k];
        const int rcvRow = fullIdx[j];
        Success[rcvRow] = true;
        if (verbose >= 1) DistOfSuccess[rcvRow] = curDist;
        dropped[j] = 1;
        --numLive;
        anySuccess = true;
      }
      lastWardMaxProb = 1.0; // first iteration, anything possible later
    } else {
      for (int k = 0; k < sumHasSp; ++k) {
        if (ran[k] > lastWardMaxProb) continue;
        const int j  = hasSpRows[k];
        const int sp = spAct[j];
        const double w = wardForRow[sp];
        if (w > newMax) newMax = w;
        if (ran[k] < w) {
          const int rcvRow = fullIdx[j];
          Success[rcvRow] = true;
          if (verbose >= 1) DistOfSuccess[rcvRow] = curDist;
          dropped[j] = 1;
          --numLive;
          anySuccess = true;
        }
      }
      // Update lastWardMaxProb to min(1, max(wardRes)) across the surviving
      // receivers from this iteration's pre-screened set, mirroring R.
      // R only updates lastWardMaxProb when length(whRanLTprevMaxProb) > 0,
      // which we already established with anyPass.
      if (newMax > 1.0) newMax = 1.0;
      lastWardMaxProb = newMax;
    }

    // If any receivers were marked dropped this step (success path), compact
    // the active list so subsequent steps don't iterate them. Single in-place
    // pass; preserves order of remaining indices.
    if (anySuccess) {
      size_t w = 0;
      for (size_t r = 0; r < active.size(); ++r) {
        if (!dropped[active[r]]) {
          active[w++] = active[r];
        }
      }
      active.resize(w);
    }

    if (numLive == 0) break;
  }

  PutRNGstate();

  return List::create(
    _["Success"]       = Success,
    _["DistOfSuccess"] = DistOfSuccess
  );
}
