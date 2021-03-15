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
//' This uses a spiral pattern outwards from the \code{receiveCellCoords} cells on
//' a raster with the dimensions \code{numCols}, \code{numCells}, \code{xmin},
//' \code{ymin}, and \code{cellSize}. For each cell in \code{receiveCellCoords},
//' it evaluates whether there is a successful "dispersal" \code{to} that cell
//' for the species that can disperse there as identified by \code{rcvSpeciesByIndex}.
//' It will search outwards testing each and every cell in the spiral until
//' the maximum distance is reached as specified in the 3rd column (named or unnamed)
//' of \code{speciesTable}.
//'
//' @param receiveCellCoords Matrix, 2 columns, of x-y coordinates of the Receive cells
//' @param srcListVectorBySp A list, where each element is a vector of NA and the speciesCode
//'   value of the list. The length of each vector MUST be the number of cells in the
//'   raster whose \code{receiveCellCoords} are provided.
//' @param rcvSpeciesByIndex A list of length \code{NROW(receiveCellCoords)} where each element
//'   is the vector of speciesCodes that are capable of being received in the
//'   corresponding \code{receiveCellCoords}
//' @param speciesTable A numeric matrix with species traits. Must have column 3 be
//'   \code{seeddistance_max}, column 2 be \code{seeddistance_eff}, and sorted in
//'   increasing order on the first column, speciesCode. The speciesCode values must
//'   be \code{seq(1, NROW(speciesTable))}. The names of these columns is not important,
//'   only the position in the matrix
//' @param numRcvSpeciesVec A vector of length and order of NROW(speciesTable) indicating
//'   how many of each receiving species there are. This will be decremented. If it reaches
//'   0, it will allow shortcutting of that species.
//' @param numCols Integer, number of columns in the raster whose \code{receiveCellCoords}
//'   were provided
//' @param numRows Integer, number of rows in the raster whose \code{receiveCellCoords}
//'   were provided
//' @param numCells Integer, number of cells in the raster whose \code{receiveCellCoords}
//'   were provided
//' @param cellSize Integer, the \code{res(ras)[1]} of the raster whose \code{receiveCellCoords}
//'   were provided
//' @param xmin Integer, the \code{xmin(ras)} of the raster whose \code{receiveCellCoords}
//'   were provided
//' @param ymin Integer, the \code{ymin(ras)} of the raster whose \code{receiveCellCoords}
//'   were provided
//' @param k Numeric, parameter passed to Ward dispersal kernel
//' @param b Numeric, parameter passed to Ward dispersal kernel
//' @param successionTimestep Integer, Same as Biomass_core.
//' @param verbose Numeric, length 1. Currently \code{0} (no messaging), the fastest option,
//'   \code{1} (some messaging) and \code{2} or greater (more messaging) are active. Default is
//'   \code{getOption("LandR.verbose", TRUE)}.
//' @return A logical matrix with ncols = \code{length(srcListVectorBySp)} and nrows =
//'   \code{NROW(receiveCellCoords)}, indicating whether that receiveCellCoords successfully
//'   received seeds from each species.
//' @author Eliot McIntire
//' @export
// [[Rcpp::export]]
LogicalMatrix spiralSeedDispersal( IntegerMatrix receiveCellCoords,
                                   Rcpp::List srcListVectorBySp, List rcvSpeciesByIndex,
                                   NumericMatrix speciesTable,
                                   IntegerVector numRcvSpeciesVec,
                                   int numCols, int numRows, int numCells, int cellSize,
                                   int xmin, int ymin,
                                   double k, double b, double successionTimestep,
                                   double verbose = 0.0)
{

  int nCellsRcv(receiveCellCoords.nrow());
  int nSpeciesEntries(speciesTable.nrow());
  int x, y, dx, dy, spiralIndex;
  // int spiralIndexMax = maxSpiral; // not really used except for debugging, it can be shrunk
  bool underMaxDist = true;

  // max distances by species
  NumericVector maxDist(1);
  NumericVector maxDistMinCellSize(1);
  NumericVector effDist(1);
  NumericVector dis(1);
  NumericVector maxDistsSpV = speciesTable(_, 2) ;
  NumericVector maxDistsSpVMinCellSize = pmax( speciesTable(_, 2), cellSize * 1.0 ) ;
  double overallMaxDist = max(maxDistsSpV);
  double overallMaxDistMinCellSize = max(maxDistsSpVMinCellSize);
  overallMaxDist = std::max( overallMaxDist, cellSize * 1.0 );
  double overallMaxDistCorner = overallMaxDist * sqrt(2);
  NumericVector effDistsSpV = speciesTable(_, 1);
  NumericVector effDistNotNA;
  NumericVector maxDistNotNA;
  // NumericVector numActiveCellsByRcvSp(nSpeciesEntries);
  // NumericVector numActiveCellsByRcvSpDone(nSpeciesEntries);
  double maxOfMaxDists, newOverallMaxDistCorner;

  IntegerVector speciesPixelRcvPool;
  IntegerVector cellRcvPool = Rcpp::seq(0, receiveCellCoords.nrow()-1);

  // coordinates and distances
  NumericVector xCoord;
  NumericVector yCoord;
  NumericVector dis1;
  int pixelSrc, pixelVal;
  // bool notNegative, notTooBig,
  bool inequ;

  // output
  LogicalMatrix seedsArrivedMat(nCellsRcv, nSpeciesEntries);
  LogicalVector rcvSpDone (nSpeciesEntries, 0);
  LogicalVector notYetOverMaxDist (nSpeciesEntries, 1);

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
  int nCellsVisitedOnMap = 0;
  int nCellsVisited = 0;

  // Primary "spiral" loop, one pixel at a time. Only calculate for the "generic" pixel,
  // which is effectively an offset from the "central" pixel;
  // then add this offset to cellCoods matrix of initial cells.
  // This will create a square-ish shape, i.e., make a square then add a single
  // pixel width around entire square to make a new slightly bigger square.
  while( (underMaxDist == true) ) { // } && (nCellsVisited < 20) ) {
    Rcpp::checkUserInterrupt();
    spiralIndex += 1;
    // NumericVector numActiveCellsByRcvSp(nSpeciesEntries); // need to rezero
    // numActiveCellsByRcvSp = numActiveCellsByRcvSp + numActiveCellsByRcvSpDone;

    xCoord = x * cellSize + receiveCellCoords(_, 0);
    yCoord = y * cellSize + receiveCellCoords(_, 1);
    dis1[0] = sqrt(pow(receiveCellCoords(0, 0) - xCoord[0], 2.0) +
      pow(receiveCellCoords(0, 1) - yCoord[0], 2.0) );

    LogicalVector overMaxDists = (dis1[0]  > maxDistsSpVMinCellSize * sqrt(2)) * notYetOverMaxDist;

    if (is_true(any(overMaxDists))) {
      notYetOverMaxDist = ( 1 - overMaxDists ) * notYetOverMaxDist ;
      maxDistsSpV[overMaxDists] = 0;
      maxDistsSpVMinCellSize[overMaxDists] = 0;
    }

    // if (verbose >= 3) {
    //   Rcpp::Rcout << " ---------------------------------------- spiralIndex " << spiralIndex <<  " dis1 " << dis1[0] << " overallMaxDist: " << overallMaxDist << std::endl;
    //   Rcpp::Rcout << "overMaxDists: " <<  overMaxDists << " notYetOverMaxDist " << notYetOverMaxDist  << std::endl;
    // }


    ////////////////////////////////////
    // messaging for progress
    // if (verbose >= 3) {
    //   disInt = floor(dis1[0]/ sqrt(2));
    //   possCurModVal = disInt % moduloVal;
    //   possCurMessage = floor(disInt / moduloVal) * moduloVal;
    //   if (possCurModVal < curModVal && possCurMessage > curMessage)  {
    //     curMessage = possCurMessage;
    //     Rcpp::Rcout << "Dispersal distance completed: " << curMessage << " of " << overallMaxDist << std::endl;
    //   }
    //   curModVal = possCurModVal;
    // }
    // End messaging for progress
    ////////////////////////////////////////////////////

    if (verbose >= 3) {
      Rcpp::Rcout << "overallMaxDist: " << overallMaxDist << " dis1[0] " << dis1[0] << std::endl;
    }

    if (dis1[0] <= ( overallMaxDist ) ) { // make sure to omit the corners of the square due to circle
      if (verbose >= 4) {
        Rcpp::Rcout << "overallMaxDist: " << overallMaxDist << " -- " << dis1[0] << std::endl;
      }

      // Loop around each of the original cells

      for (IntegerVector::iterator cellRcvIt = cellRcvPool.begin();
           cellRcvIt != cellRcvPool.end(); ) {

        // for (int cellRcvInd = 0; cellRcvInd < nCellsRcv; ++cellRcvInd) {



        IntegerVector speciesPixelRcvPool = rcvSpeciesByIndex[*cellRcvIt];

        nCellsVisited += 1;
        if (verbose >= 3) {
          Rcpp::Rcout << "nCellsVisited " << nCellsVisited << std::endl;
          Rcpp::Rcout << "speciesPixelRcvPool " << speciesPixelRcvPool << " cellRcvPool " << cellRcvPool <<  std::endl;
          Rcpp::Rcout << "*cellRcvIt " << *cellRcvIt << "; nCellsRcv " << nCellsRcv << std::endl;
          Rcpp::Rcout << " xCoord[*cellRcvIt] " << xCoord[*cellRcvIt] << " yCoord[*cellRcvIt] " << yCoord[*cellRcvIt] <<  std::endl;
          Rcpp::Rcout << "****** dis1: " << dis1[0] << "  maxDistsSpV " << maxDistsSpV << ";; maxDistsSpVMinCellSize: " << maxDistsSpVMinCellSize << std::endl;
          Rcpp::Rcout << "pixelVal " << pixelVal << std::endl;
        }
        if (speciesPixelRcvPool.length() > 0) {

          bool onMap = (xCoord[*cellRcvIt] < numCols * cellSize + xmin) && (xCoord[*cellRcvIt] > xmin) &&
            (yCoord[*cellRcvIt] < numRows * cellSize + ymin) && (yCoord[*cellRcvIt] > ymin);
          if (onMap) {
            nCellsVisitedOnMap += 1;
            if (verbose >= 3) {
              Rcpp::Rcout << "#### nCellsVisitedOnMap " << nCellsVisitedOnMap  << std::endl;
            }
            pixelSrc = (numRows - (yCoord[*cellRcvIt] - ymin - cellSize/2)/cellSize - 1) * numCols +
              (xCoord[*cellRcvIt] - xmin - cellSize/2)/cellSize + 1;
            if (verbose >= 3) {
              Rcpp::Rcout << "pixelSrc " << pixelSrc << " maxDist " << maxDist << std::endl;
            }


            for (IntegerVector::iterator speciesPixelRcv = speciesPixelRcvPool.begin();
                 speciesPixelRcv != speciesPixelRcvPool.end(); ) {

              maxDist = maxDistsSpV[*speciesPixelRcv - 1];
              maxDistMinCellSize = maxDistsSpVMinCellSize[*speciesPixelRcv - 1];
              if (dis1[0] > ( maxDistMinCellSize[0] * sqrt(2) ) ) {
                // remove that species from the species pool for that Rcv cell
                speciesPixelRcv = speciesPixelRcvPool.erase(speciesPixelRcv);
                rcvSpeciesByIndex[*cellRcvIt] = speciesPixelRcvPool;

              } else { // within the square
                if (dis1[0] <= ( maxDistMinCellSize[0]  ) ) {
                  // make sure to omit the corners of the square due to circle
                  // if (dis1[0] <= std::max( maxDist[0], cellSize * 1.0  ) ) { // make sure to omit the corners of the square due to circle
                  alreadyReceived = seedsArrivedMat(*cellRcvIt, *speciesPixelRcv - 1);
                  if (!alreadyReceived) {

                    // pixelSrc = (numCols - ((yCoord[*cellRcvIt] - ymin + cellSize/2)/cellSize)) * numCols +
                    //   (xCoord[*cellRcvIt] - xmin + cellSize/2)/cellSize;

                    IntegerVector speciesVector = srcListVectorBySp[*speciesPixelRcv - 1];
                    pixelVal = speciesVector[pixelSrc - 1];
                    // pixelVal = speciesVector[pixelSrc];
                    // numActiveCellsByRcvSp[*speciesPixelRcv - 1] += 1;

                    inequ = 0; // default
                    if (pixelVal >= 0) { // covers NA which is -2147483648
                      effDist = effDistsSpV[*speciesPixelRcv - 1];
                      dis = dis1[0];

                      if (is_true(all(dis <= pmax( maxDist, cellSize * 1.0  ) ))) {
                        if (dis[0] == 0.0) {
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
                        seedsArrivedMat(*cellRcvIt, *speciesPixelRcv - 1) = inequ || seedsArrivedMat(*cellRcvIt, *speciesPixelRcv - 1);

                        if (inequ) {
                          if (verbose >= 4) {
                            Rcpp::Rcout << "  Success! " << " speciesPixelRcv " << *speciesPixelRcv  << std::endl;
                          }

                          numRcvSpeciesVec[*speciesPixelRcv - 1] = numRcvSpeciesVec[*speciesPixelRcv - 1] - 1;

                            if (verbose >= 4) {
                              Rcpp::Rcout << "&&&&&&&&&&&&& numRcvSpeciesVec " << numRcvSpeciesVec << std::endl;
                            }

                          // numActiveCellsByRcvSp[*speciesPixelRcv - 1] = numActiveCellsByRcvSp[*speciesPixelRcv - 1] - 1;
                          // Remove this one as it is no longer needed
                          if (verbose >= 3) {
                            Rcpp::Rcout << "speciesPixelRcvPool - BEFORE " << speciesPixelRcvPool << "   speciesPixelRcv " << *speciesPixelRcv << std::endl;
                          }
                          speciesPixelRcv = speciesPixelRcvPool.erase(speciesPixelRcv);
                          if (verbose >= 3) {
                            Rcpp::Rcout << "speciesPixelRcvPool - AFTER " << speciesPixelRcvPool  << "   speciesPixelRcv " << *speciesPixelRcv << std::endl;
                          }
                          rcvSpeciesByIndex[*cellRcvIt] = speciesPixelRcvPool;

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

        if (speciesPixelRcvPool.length() == 0L) {
          cellRcvIt = cellRcvPool.erase(cellRcvIt);
        } else {
          ++cellRcvIt;
        }


      }

      // Pre-empt further searching if all rcv cells have received a given species
      // If, say, all rcv cells that can get species 2 all have gotten species 2, then
      // can remove species 2 and lower the overallMaxDistCorner accordingly
      // LogicalVector rcvSpDone = numActiveCellsByRcvSp == 0;
      LogicalVector rcvSpDone = numRcvSpeciesVec == 0L;

      if (is_true(any(rcvSpDone))) {
        // numActiveCellsByRcvSpDone = -1 * rcvSpDone;
        maxDistsSpV[rcvSpDone] = 0;
        maxDistsSpVMinCellSize[rcvSpDone] = 0;
        numRcvSpeciesVec[rcvSpDone] = -1; // make it out of contention for future assessments
        maxOfMaxDists = max(maxDistsSpV);
        newOverallMaxDistCorner = maxOfMaxDists * sqrt(2);
        IntegerVector rcvSpCodesDone = which2(rcvSpDone) + 1;
        if (verbose >= 4) {
          Rcpp::Rcout << "Species " << rcvSpCodesDone << " complete. New max dispersal distance: " << maxOfMaxDists << std::endl;
        }
        overallMaxDistCorner = newOverallMaxDistCorner;
        overallMaxDist = maxOfMaxDists;
      }

    } // end skip corners of distance

    if (dis1[0] > ( overallMaxDistCorner ) ) {
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
  if (verbose >= 3) {
    Rcpp::Rcout << "nCellsVisited " << nCellsVisited << " nCellsVisitedOnMap: " << nCellsVisitedOnMap << std::endl;
  }
  return seedsArrivedMat;
  // return dis1;
}

