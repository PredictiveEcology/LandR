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
  std::vector<int> cellRcvPool2 = as< std::vector<int> >(cellRcvPool);

  IntegerVector cellsWithSpNext (receiveCellCoords.nrow(), -100); // (nCellsWithSp);
  int nCellsWithSp = 0;
  int nCellsWithSpPrev = receiveCellCoords.nrow();

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
  int cellRcvCounter = 0;
  int currentRing = 1;
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
  while( (underMaxDist == true) ) { //} && spiralIndex < 10 ) { // } && (nCellsVisited < 20) ) {
    if (verbose >= 2) {
      Rcpp::Rcout << " &&&&&&&&&&& cellRcvPool.size() " << cellRcvPool.size() << " nCellsWithSp " << nCellsWithSp << " cellRcvPool2.size() " << cellRcvPool2.size() << std::endl;
    }

    Rcpp::checkUserInterrupt();
    spiralIndex += 1;
    // NumericVector numActiveCellsByRcvSp(nSpeciesEntries); // need to rezero
    // numActiveCellsByRcvSp = numActiveCellsByRcvSp + numActiveCellsByRcvSpDone;

    xCoord = x * cellSize + receiveCellCoords(_, 0);
    yCoord = y * cellSize + receiveCellCoords(_, 1);
    dis1[0] = sqrt(pow(receiveCellCoords(0, 0) - xCoord[0], 2.0) +
      pow(receiveCellCoords(0, 1) - yCoord[0], 2.0) );

    LogicalVector overMaxDists = (dis1[0]  > maxDistsSpVMinCellSize * sqrt(2)) * notYetOverMaxDist;
    if (verbose >= 3) {
      Rcpp::Rcout << " ---------------------------------------- spiralIndex " << spiralIndex <<  " dis1 " << dis1[0] << " overallMaxDist: " << overallMaxDist << std::endl;
    }

    bool anyNewOverMax = is_true(any(overMaxDists));

    if (anyNewOverMax) {
      notYetOverMaxDist = ( 1 - overMaxDists ) * notYetOverMaxDist ;
      maxDistsSpV[overMaxDists] = 0;
      maxDistsSpVMinCellSize[overMaxDists] = 0;
      if (verbose >= 3) {
        Rcpp::Rcout << "overMaxDists: " <<  overMaxDists << " notYetOverMaxDist " << notYetOverMaxDist  << " maxDistsSpVMinCellSize " <<  maxDistsSpVMinCellSize << std::endl;
      }
    }

    // if (verbose >= 3) {
    //   Rcpp::Rcout << " ---------------------------------------- spiralIndex " << spiralIndex <<  " dis1 " << dis1[0] << " overallMaxDist: " << overallMaxDist << std::endl;
    //   Rcpp::Rcout << "overMaxDists: " <<  overMaxDists << " notYetOverMaxDist " << notYetOverMaxDist  << std::endl;
    // }


    if (dis1[0] <= ( overallMaxDist ) ) { // make sure to omit the corners of the square due to circle
      if (verbose >= 3) {
        Rcpp::Rcout << "overallMaxDist: " << overallMaxDist << " dis1[0] " << dis1[0] << " cellRcvPool.size() " << cellRcvPool.size() << std::endl;
      }

      // Loop around each of the original cells

      if (nCellsWithSpPrev == 0L ) {
        if (verbose >= 3) {
          Rcpp::Rcout << "LALALALALALALAL" << std::endl;
        }

        }
      for (int cellRcvInt = 0; cellRcvInt != nCellsWithSpPrev; ++cellRcvInt ) {
      // for (IntegerVector::iterator cellRcvIt = cellRcvPool.begin();
      //      cellRcvIt != nCellsWithSpPrev; ) {
        IntegerVector speciesPixelRcvPool = rcvSpeciesByIndex[cellRcvInt];
        if (verbose >= 4) {
          Rcpp::Rcout << "               cellRcvInt: " << cellRcvInt << " speciesPixelRcvPool: " << speciesPixelRcvPool << std::endl;
        }

        nCellsVisited += 1;
        if (speciesPixelRcvPool.length() > 0) {

          bool onMap = (xCoord[cellRcvInt] < numCols * cellSize + xmin) && (xCoord[cellRcvInt] > xmin) &&
            (yCoord[cellRcvInt] < numRows * cellSize + ymin) && (yCoord[cellRcvInt] > ymin);
          if (verbose >= 4) {
            Rcpp::Rcout << "onMap " << onMap << std::endl;
          }

          if (onMap) {
            nCellsVisitedOnMap += 1;
            if (verbose >= 4) {
              Rcpp::Rcout << "#### nCellsVisitedOnMap " << nCellsVisitedOnMap  << std::endl;
            }
            pixelSrc = (numRows - (yCoord[cellRcvInt] - ymin - cellSize/2)/cellSize - 1) * numCols +
              (xCoord[cellRcvInt] - xmin - cellSize/2)/cellSize + 1;
            if (verbose >= 4) {
              Rcpp::Rcout << "pixelSrc " << pixelSrc << std::endl;
            }


            for (IntegerVector::iterator speciesPixelRcv = speciesPixelRcvPool.begin();
                 speciesPixelRcv != speciesPixelRcvPool.end(); ) {

              maxDist = maxDistsSpV[*speciesPixelRcv - 1];
              maxDistMinCellSize = maxDistsSpVMinCellSize[*speciesPixelRcv - 1];
              if (verbose >= 4) {
                Rcpp::Rcout << "maxDist " << maxDist << " *speciesPixelRcv " << *speciesPixelRcv << std::endl;
              }

              if (dis1[0] > ( maxDistMinCellSize[0] * sqrt(2) ) ) {
                // remove that species from the species pool for that Rcv cell
                if (verbose >= 4) {
                  Rcpp::Rcout << "Removing cellRcvInt: "  << cellRcvInt << "  *speciesPixelRcv " << *speciesPixelRcv << " speciesPixelRcvPool " << speciesPixelRcvPool << std::endl;
                }
                speciesPixelRcv = speciesPixelRcvPool.erase(speciesPixelRcv);
                rcvSpeciesByIndex[cellRcvInt] = speciesPixelRcvPool;
                if (verbose >= 4) {
                  Rcpp::Rcout << "Removed cellRcvInt: "  << cellRcvInt << "  *speciesPixelRcv " << *speciesPixelRcv << " speciesPixelRcvPool " << speciesPixelRcvPool << std::endl;
                }

              } else { // within the square
                if (dis1[0] <= ( maxDistMinCellSize[0]  ) ) {
                  // make sure to omit the corners of the square due to circle
                  alreadyReceived = seedsArrivedMat(cellRcvInt, *speciesPixelRcv - 1);
                  if (!alreadyReceived) {

                    IntegerVector speciesVector = srcListVectorBySp[*speciesPixelRcv - 1];
                    pixelVal = speciesVector[pixelSrc - 1];
                    if (verbose >= 4) {
                      Rcpp::Rcout << "pixelVal: "  << pixelVal << std::endl;
                    }

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
                        seedsArrivedMat(cellRcvInt, *speciesPixelRcv - 1) = inequ || seedsArrivedMat(cellRcvInt, *speciesPixelRcv - 1);

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
                          if (verbose >= 4) {
                            Rcpp::Rcout << "speciesPixelRcvPool - BEFORE " << speciesPixelRcvPool << "   speciesPixelRcv " << *speciesPixelRcv << std::endl;
                          }
                          speciesPixelRcv = speciesPixelRcvPool.erase(speciesPixelRcv);
                          if (verbose >= 4) {
                            Rcpp::Rcout << "speciesPixelRcvPool - AFTER " << speciesPixelRcvPool  << "   speciesPixelRcv " << *speciesPixelRcv << std::endl;
                          }
                          rcvSpeciesByIndex[cellRcvInt] = speciesPixelRcvPool;

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


        if (speciesPixelRcvPool.length() > 0L) {
          cellsWithSpNext[nCellsWithSp] = cellRcvInt;
          nCellsWithSp += 1;
        }
        if (verbose >= 3) {
          if (cellRcvInt < 3) {
            Rcpp::Rcout << "- nCellsWithSp " << nCellsWithSp << " cellRcvInt " << cellRcvInt << " leaving cellRcvPool.size() " << cellRcvPool.size() << std::endl;
          }
        }

        // ++cellRcvIt;
      }
      if (verbose >= 3) {
        Rcpp::Rcout << "Done dis1[0] " << dis1[0] << std::endl;
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

    if ( (dis1[0] > ( overallMaxDistCorner ) ) ) {
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
    ////////////////////////////////////
    // messaging for progress
    if (verbose >= 3) {
      Rcpp::Rcout << "B LALALA " << dis1[0] << std::endl;
    }
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
    if (verbose >= 3) {
      Rcpp::Rcout << "C LALALA " << dis1[0] << std::endl;
    }
    // End messaging for progress
    ////////////////////////////////////////////////////
    if (dis1[0] > (currentRing * cellSize * sqrt(2)) && nCellsWithSp > 0) {
      currentRing += 1;
      // IntegerVector cellRcvPool(nCellsWithSp);
      IntegerVector cellRcvPoolInd = Rcpp::seq(0, nCellsWithSp - 1);
      if (verbose >= 2) {
        Rcpp::Rcout << "cellRcvPoolInd.length() " << cellRcvPoolInd.length() << std::endl;
      }

      IntegerVector cellsWithSpNextNotNA = cellsWithSpNext[cellRcvPoolInd];

      // if (verbose >= 2) {
      //   Rcpp::Rcout << "starting next ring " << currentRing << " cellRcvPool.size() " << cellRcvPool.size() << " cellsWithSpNext.length() " << cellsWithSpNext.length()  << " nCellsWithSp" << nCellsWithSp << std::endl;
      // }

      // std::vector<int> tmp = as< std::vector<int> >(cellsWithSpNext[cellsWithSpNextNotNA])
      std::vector<int> cellRcvPool2 = as< std::vector<int> >(cellsWithSpNext[cellsWithSpNextNotNA]);
      if (verbose >= 2) {
        Rcpp::Rcout << "cellRcvPool2.size() " << cellRcvPool2.size() << std::endl;
      }
      cellRcvPool2.resize(nCellsWithSp);
      std::vector<double> yyy(nCellsWithSp);
      IntegerVector cellRcvPool(yyy.begin(), yyy.end());

      if (verbose >= 2) {
        Rcpp::Rcout << "cellRcvPool2.size() " << cellRcvPool2.size() << std::endl;
      }

      // IntegerVector cellRcvPool = cellsWithSpNext[cellsWithSpNextNotNA];
      // cellRcvPool.resize(nCellsWithSp);
      cellRcvPool = cellsWithSpNext[cellsWithSpNextNotNA]; // move shorter dataset into cellRcvPool
      IntegerVector cellsWithSpNext(nCellsWithSp); // rezero cellsWithSpNext
      nCellsWithSpPrev = nCellsWithSp;
      if (verbose >= 2) {
        Rcpp::Rcout << "A starting next ring " << currentRing << " cellRcvPool.size() " << cellRcvPool.size() << " cellRcvPool.end() " << cellRcvPool.end()  << " nCellsWithSp" << nCellsWithSp << std::endl;
      }
    }
    if (verbose >= 3) {
      Rcpp::Rcout << "D LALALA " << dis1[0] << std::endl;
    }

      if (nCellsWithSpPrev < 1) {
        underMaxDist = false;
      }

    if (verbose >= 2) {
      Rcpp::Rcout << "A2 starting next ring " << currentRing << " cellRcvPool.size() " << cellRcvPool.size() << " cellRcvPool.end() " << cellRcvPool.end()  << " nCellsWithSp" << nCellsWithSp << std::endl;
    }
    nCellsWithSp = 0;

    x += dx;
    y += dy;
    if (verbose >= 2) {
      Rcpp::Rcout << "B starting next ring " << currentRing << " cellRcvPool.size() " << cellRcvPool.size() << " cellRcvPool.end() " << cellRcvPool.end()  << " nCellsWithSp" << nCellsWithSp << std::endl;
    }

  }
  if (verbose >= 2) {
    Rcpp::Rcout << "C starting next ring " << currentRing << " cellRcvPool.size() " << cellRcvPool.size() << " cellRcvPool.end() " << cellRcvPool.end()  << " nCellsWithSp" << nCellsWithSp << std::endl;
  }


  // return List::create(
  //   _["x"] = xKeep,
  //   _["y"] = yKeep,
  //   _["dis"] = disKeep,
  //   _["spiralIndexKeep"] = spiralIndexKeep,
  //   _["dispersalProbKeep"] = dispersalProbKeep,
  //   _["seedsArrived"] = seedsArrivedMat
  // );
  if (verbose >= 2) {
    Rcpp::Rcout << "nCellsVisited " << nCellsVisited << " nCellsVisitedOnMap: " << nCellsVisitedOnMap << std::endl;
  }
  return seedsArrivedMat;
  // return dis1;
}


// IntegerVector rcppSeq = Rcpp::seq(0, receiveCellCoords.nrow()-1);
// std::vector<int> cellRcvPool = as< std::vector<int> >(rcppSeq);

