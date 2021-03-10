#include <numeric>
// #include <math.h> /* round, floor, ceil, trunc */
#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
Rcpp::IntegerVector which2(Rcpp::LogicalVector x) {
  Rcpp::IntegerVector v = Rcpp::seq(0, x.size()-1);
  return v[x];
}

// [[Rcpp::export]]
Rcpp::IntegerVector rmElem(Rcpp::IntegerVector x, IntegerVector toRm) {
  for(IntegerVector::iterator toRmIt = toRm.begin(); toRmIt != toRm.end(); ++toRmIt) {
    for (IntegerVector::iterator xIt = x.begin(); xIt != x.end(); ) {
      if (*xIt == *toRmIt) {
        xIt = x.erase(xIt);
      } else {
        ++xIt;
      }
    }
  }
  return x;
}



//' Ward seed dispersal using Rcpp
//'
//' This uses a spiral pattern outwards from the \code{cellCoords} cells on
//' a raster with the dimensions \code{numCols}, \code{numCells}, \code{xmin},
//' \code{ymin}, and \code{cellSize}. For each cell in \code{cellCoords},
//' it evaluates whether there is a successful "dispersal" \code{to} that cell
//' for the species that can disperse there as identified by \code{rcvSpeciesByIndex}.
//' It will search outwards testing each and every cell in the spiral until
//' the maximum distance is reached as specified in the 3rd column (named or unnamed)
//' of \code{speciesTable}.
//'
//' @param cellCoords Matrix, 2 columns, of x-y coordinates of the Receive cells
//' @param speciesVectorsList A list, where each element is a vector of NA and the speciesCode
//'   value of the list. The length of each vector MUST be the number of cells in the
//'   raster whose \code{cellCoords} are provided.
//' @param rcvSpeciesByIndex A list of length \code{NROW(cellCoords)} where each element
//'   is the vector of speciesCodes that are capable of being received in the
//'   corresponding \code{cellCoords}
//' @param speciesTable A numeric matrix with species traits. Must have column 3 be
//'   \code{seeddistance_max}, column 2 be \code{seeddistance_eff}, and sorted in
//'   increasing order on the first column, speciesCode. The speciesCode values must
//'   be \code{seq(1, NROW(speciesTable))}. The names of these columns is not important,
//'   only the position in the matrix
//' @param numCols Integer, number of columns in the raster whose \code{cellCoords}
//'   were provided
//' @param numRows Integer, number of rows in the raster whose \code{cellCoords}
//'   were provided
//' @param numCells Integer, number of cells in the raster whose \code{cellCoords}
//'   were provided
//' @param cellSize Integer, the \code{res(ras)[1]} of the raster whose \code{cellCoords}
//'   were provided
//' @param xmin Integer, the \code{xmin(ras)} of the raster whose \code{cellCoords}
//'   were provided
//' @param ymin Integer, the \code{ymin(ras)} of the raster whose \code{cellCoords}
//'   were provided
//' @param k Numeric, parameter passed to Ward dispersal kernel
//' @param b Numeric, parameter passed to Ward dispersal kernel
//' @param successionTimestep Integer, Same as Biomass_core.
//' @param verbose Numeric, length 1. Currently \code{0} (no messaging), the fastest option,
//'   \code{1} (some messaging) and \code{2} or greater (more messaging) are active. Default is
//'   \code{getOption("LandR.verbose", TRUE)}.
//' @return A logical matrix with ncols = \code{length(speciesVectorsList)} and nrows =
//'   \code{NROW(cellCoords)}, indicating whether that cellCoords successfully
//'   received seeds from each species.
//' @author Eliot McIntire
//' @export
// [[Rcpp::export]]
LogicalMatrix spiralSeedDispersal( IntegerMatrix cellCoords,
                                   Rcpp::List speciesVectorsList, List rcvSpeciesByIndex,
                                   NumericMatrix speciesTable,
                                   int numCols, int numRows, int numCells, int cellSize,
                                   int xmin, int ymin,
                                   double k, double b, double successionTimestep,
                                   double verbose = 0.0)
{

  int nCellsRcv(cellCoords.nrow());
  int nSpeciesEntries(speciesTable.nrow());
  int x, y, dx, dy, spiralIndex;
  // int spiralIndexMax = maxSpiral; // not really used except for debugging, it can be shrunk
  bool underMaxDist = true;

  // max distances by species
  NumericVector maxDist(1);
  NumericVector effDist(1);
  NumericVector dis(1);
  NumericVector maxDistsSpV = speciesTable(_, 2);
  double overallMaxDist = max(maxDistsSpV);
  double overallMaxDistCorner = overallMaxDist * sqrt(2);
  NumericVector effDistsSpV = speciesTable(_, 1);
  NumericVector effDistNotNA;
  NumericVector maxDistNotNA;
  NumericVector numActiveCellsByRcvSp(nSpeciesEntries);
  NumericVector numActiveCellsByRcvSpDone(nSpeciesEntries);
  double maxOfMaxDists, newOverallMaxDistCorner;

  IntegerVector speciesPixelRcvPool;

  // coordinates and distances
  NumericVector xCoord;
  NumericVector yCoord;
  NumericVector dis1;
  int pixelSrc, pixelVal;
  // bool notNegative, notTooBig,
  bool inequ;

  // output
  LogicalMatrix seedsArrivedMat(nCellsRcv, nSpeciesEntries);

  // dispersal specific stuff
  double ran;
  bool alreadyReceived;

  // Spiral mechanism
  x = y = dx = 0;
  spiralIndex = -1;
  dy = -1;
  int t = overallMaxDist;

  // messaging
  int floorOverallMaxDist = floor(overallMaxDist / 10);
  int moduloVal;
  if (cellSize < floorOverallMaxDist) {
    moduloVal = floorOverallMaxDist;
  } else {
    moduloVal = cellSize;
  }
  int curModVal = moduloVal;
  int curMessage = 0;
  double possCurModVal;
  int possCurMessage, disInt;

  // Primary "spiral" loop, one pixel at a time. Only calculate for the "generic" pixel,
  // which is effectively an offset from the "central" pixel;
  // then add this offset to cellCoods matrix of initial cells.
  // This will create a square-ish shape, i.e., make a square then add a single
  // pixel width around entire square to make a new slightly bigger square.
  while(underMaxDist == true) { // && spiralIndex < spiralIndexMax) {
    Rcpp::checkUserInterrupt();
    spiralIndex += 1;
    NumericVector numActiveCellsByRcvSp(nSpeciesEntries); // need to rezero
    numActiveCellsByRcvSp = numActiveCellsByRcvSp + numActiveCellsByRcvSpDone;

    xCoord = x * cellSize + cellCoords(_, 0);
    yCoord = y * cellSize + cellCoords(_, 1);
    dis1[0] = sqrt(pow(cellCoords(0, 0) - xCoord[0], 2.0) +
      pow(cellCoords(0, 1) - yCoord[0], 2.0) );
    // if (verbose >= 1) {
    //   Rcpp::Rcout << " --- 000 --- " << std::endl;
    //   Rcpp::Rcout << "dis1 " << dis1[0] << std::endl;
    // }


    ////////////////////////////////////
    // messaging for progress
    if (verbose >= 1) {
      disInt = floor(dis1[0]/ sqrt(2));
      possCurModVal = disInt % moduloVal;
      possCurMessage = floor(disInt / moduloVal) * moduloVal;
      if (possCurModVal < curModVal && possCurMessage > curMessage)  {
        curMessage = possCurMessage;
        Rcpp::Rcout << "Dispersal distance completed: " << curMessage << " of " << overallMaxDist << std::endl;
      }
      curModVal = possCurModVal;
    }
    // End messaging for progress
    ////////////////////////////////////////////////////

    if (dis1[0] <= ( overallMaxDist + cellSize) ) { // make sure to omit the corners of the square due to circle

      // Loop around each of the original cells
      for (int cellRcvInd = 0; cellRcvInd < nCellsRcv; ++cellRcvInd) {
        IntegerVector speciesPixelRcvPool = rcvSpeciesByIndex[cellRcvInd];

        if (speciesPixelRcvPool.length() > 0) {

          bool onMap = (xCoord[cellRcvInd] < numCols * cellSize + xmin) && (xCoord[cellRcvInd] > xmin) &&
            (yCoord[cellRcvInd] < numRows * cellSize + ymin) && (yCoord[cellRcvInd] > ymin);
          if (onMap) {
            // if (verbose >= 1) {
            //   Rcpp::Rcout << "speciesPixelRcvPool2 " << speciesPixelRcvPool << std::endl;
            // }

            pixelSrc = (numRows - (yCoord[cellRcvInd] - ymin - cellSize/2)/cellSize - 1) * numCols +
              (xCoord[cellRcvInd] - xmin - cellSize/2)/cellSize + 1;
            // pixelSrc = (numCols - ((yCoord[cellRcvInd] - ymin + cellSize/2)/cellSize)) * numCols +
            //   (xCoord[cellRcvInd] - xmin + cellSize/2)/cellSize;

            for(IntegerVector::iterator speciesPixelRcv = speciesPixelRcvPool.begin();
                speciesPixelRcv != speciesPixelRcvPool.end(); ) {

              maxDist = maxDistsSpV[*speciesPixelRcv - 1];
              if (dis1[0] > ( maxDist[0] * sqrt(2) + cellSize) ) {
                // remove that species from the species pool for that Rcv cell
                speciesPixelRcv = speciesPixelRcvPool.erase(speciesPixelRcv);
                rcvSpeciesByIndex[cellRcvInd] = speciesPixelRcvPool;

              } else { // within the square
                if (dis1[0] <= ( maxDist[0] + cellSize  ) ) { // make sure to omit the corners of the square due to circle
                  alreadyReceived = seedsArrivedMat(cellRcvInd, *speciesPixelRcv - 1);
                  if (!alreadyReceived) {

                    IntegerVector speciesVector = speciesVectorsList[*speciesPixelRcv - 1];
                    pixelVal = speciesVector[pixelSrc - 1];
                    // pixelVal = speciesVector[pixelSrc];
                    numActiveCellsByRcvSp[*speciesPixelRcv - 1] += 1;

                    // if (verbose >= 1) {
                    //   Rcpp::Rcout << "-------- " << std::endl;
                    //   Rcpp::Rcout << "speciesPixelRcvPool " << speciesPixelRcvPool << std::endl;
                    //   Rcpp::Rcout << "speciesPixelRcv " << *speciesPixelRcv << std::endl;
                    //   Rcpp::Rcout << "cellRcvInd " << cellRcvInd << std::endl;
                    //   Rcpp::Rcout << "nCellsRcv " << nCellsRcv << std::endl;
                    //   Rcpp::Rcout << "numCols " << numCols << " cellSize " << cellSize << " xmin " << xmin << " numRows " << numRows << std::endl;
                    //   Rcpp::Rcout << " xCoord " << xCoord << " yCoord " << yCoord <<  std::endl;
                    //   Rcpp::Rcout << "pixelSrc " << pixelSrc << " maxDist " << maxDist << std::endl;
                    //   Rcpp::Rcout << "maxDistsSpV " << maxDistsSpV <<  std::endl;
                    //   Rcpp::Rcout << "dis1 " << dis1[0] << " numActiveCellsByRcvSp " << numActiveCellsByRcvSp << std::endl;
                    //   Rcpp::Rcout << "pixelVal " << pixelVal << std::endl;
                    // }

                    inequ = 0; // default
                    if (pixelVal >= 0) { // covers NA which is -2147483648
                      effDist = effDistsSpV[*speciesPixelRcv - 1];
                      dis = dis1[0];
                      // if (verbose >= 1) {
                      //   Rcpp::Rcout << "-------- " << std::endl;
                      //   Rcpp::Rcout << "dis " << dis << std::endl;
                      // }

                      if (is_true(all(dis <= ( maxDist + cellSize  ) ))) {
                        if (dis[0] == 0) {
                          inequ = 1;
                        } else {
                          // Hard coded Ward dispersal kernel -- could not figure out how to use eval(Ward, ...)
                          NumericVector dispersalProb =
                            ifelse(cellSize <= effDist,
                                   ifelse(dis <= effDist,
                                          exp((dis - cellSize) * log(1 - k)/effDist) - exp(dis * log(1 - k)/effDist),
                                          (1 - k) * exp((dis - cellSize - effDist) * log(b)/maxDist) - (1 - k) * exp((dis - effDist) * log(b)/maxDist)),
                                          ifelse(dis <= cellSize,
                                                 exp((dis - cellSize) * log(1 - k)/effDist) - (1 - k) * exp((dis - effDist) * log(b)/maxDist),
                                                 (1 - k) * exp((dis - cellSize - effDist) * log(b)/maxDist) - (1 - k) * exp((dis - effDist) * log(b)/maxDist)));
                          dispersalProb = 1 - pow(1 - dispersalProb, successionTimestep);
                          ran = R::runif(0, 1);
                          inequ = ran < dispersalProb[0];
                        }

                        // update the final matrix with a TRUE, if dispersal was successful
                        seedsArrivedMat(cellRcvInd, *speciesPixelRcv - 1) = inequ || seedsArrivedMat(cellRcvInd, *speciesPixelRcv - 1);

                        if (inequ) {
                          rcvSpeciesByIndex[cellRcvInd] = speciesPixelRcvPool;
                          numActiveCellsByRcvSp[*speciesPixelRcv - 1] = numActiveCellsByRcvSp[*speciesPixelRcv - 1] - 1;
                          // Remove this one as it is no longer needed
                          speciesPixelRcv = speciesPixelRcvPool.erase(speciesPixelRcv);

                        } else {
                          ++speciesPixelRcv;
                        }
                      } else {
                        ++speciesPixelRcv;
                      }
                    } else {
                      ++speciesPixelRcv;
                    }
                  } else {
                    ++speciesPixelRcv;
                  }
                } else {
                  ++speciesPixelRcv;
                }

              }
            }
          } else {
            // for(IntegerVector::iterator speciesPixelRcv2 = speciesPixelRcvPool.begin();
            //     speciesPixelRcv2 != speciesPixelRcvPool.end(); ) {
            //   numActiveCellsByRcvSp[*speciesPixelRcv2 - 1] == 1;
            // }
          }
        }
      }

      // Pre-empt further searching if all rcv cells have received a given species
      // If, say, all rcv cells that can get species 2 all have gotten species 2, then
      // can remove species 2 and lower the overallMaxDistCorner accordingly
      LogicalVector rcvSpDone = numActiveCellsByRcvSp == 0;
      if (is_true(any(rcvSpDone))) {
        numActiveCellsByRcvSpDone = -1 * rcvSpDone;
        maxDistsSpV[rcvSpDone] = 0;
        maxOfMaxDists = max(maxDistsSpV);
        newOverallMaxDistCorner = maxOfMaxDists * sqrt(2);
        if (newOverallMaxDistCorner < overallMaxDistCorner) {
          IntegerVector rcvSpCodesDone = which2(rcvSpDone) + 1;
          // if (verbose >= 1) {
          //   Rcpp::Rcout << "Species " << rcvSpCodesDone << " complete. New max dispersal distance: " << maxOfMaxDists << std::endl;
          // }
          overallMaxDistCorner = newOverallMaxDistCorner;
          overallMaxDist = maxOfMaxDists;
        }
      }
    } // end skip corners of distance

    if (dis1[0] > ( overallMaxDistCorner + cellSize) ) {
      underMaxDist = false;
    }


    if( (x == y) || ((x < 0) && (x == -y)) || ((x > 0) && (x == 1 - y))){
      t = dx;
      dx = -dy;
      dy = t;
    }
    // xKeep.push_back(x);
    // yKeep.push_back(y);
    // spiralIndexKeep.push_back(spiralIndex);
    // disKeep.push_back(dis1[0]);

    x += dx;
    y += dy;

  }


  // return List::create(
  //   _["x"] = xKeep,
  //   _["y"] = yKeep,
  //   _["dis"] = disKeep,
  //   _["spiralIndexKeep"] = spiralIndexKeep,
  //   _["dispersalProbKeep"] = dispersalProbKeep,
  //   _["seedsArrived"] = seedsArrivedMat
  // );
  return seedsArrivedMat;
  // return dis1;
}

