// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

// which2
Rcpp::IntegerVector which2(Rcpp::LogicalVector x);
RcppExport SEXP _LandR_which2(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::LogicalVector >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(which2(x));
    return rcpp_result_gen;
END_RCPP
}
// rmElem
Rcpp::IntegerVector rmElem(Rcpp::IntegerVector x, IntegerVector toRm);
RcppExport SEXP _LandR_rmElem(SEXP xSEXP, SEXP toRmSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::IntegerVector >::type x(xSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type toRm(toRmSEXP);
    rcpp_result_gen = Rcpp::wrap(rmElem(x, toRm));
    return rcpp_result_gen;
END_RCPP
}
// spiralSeedDispersal
LogicalMatrix spiralSeedDispersal(IntegerMatrix cellCoords, Rcpp::List speciesVectorsList, List rcvSpeciesByIndex, NumericMatrix speciesTable, int numCols, int numCells, int cellSize, int xmin, int ymin, double k, double b, double successionTimestep, double verbose, int maxSpiral);
RcppExport SEXP _LandR_spiralSeedDispersal(SEXP cellCoordsSEXP, SEXP speciesVectorsListSEXP, SEXP rcvSpeciesByIndexSEXP, SEXP speciesTableSEXP, SEXP numColsSEXP, SEXP numCellsSEXP, SEXP cellSizeSEXP, SEXP xminSEXP, SEXP yminSEXP, SEXP kSEXP, SEXP bSEXP, SEXP successionTimestepSEXP, SEXP verboseSEXP, SEXP maxSpiralSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< IntegerMatrix >::type cellCoords(cellCoordsSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type speciesVectorsList(speciesVectorsListSEXP);
    Rcpp::traits::input_parameter< List >::type rcvSpeciesByIndex(rcvSpeciesByIndexSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type speciesTable(speciesTableSEXP);
    Rcpp::traits::input_parameter< int >::type numCols(numColsSEXP);
    Rcpp::traits::input_parameter< int >::type numCells(numCellsSEXP);
    Rcpp::traits::input_parameter< int >::type cellSize(cellSizeSEXP);
    Rcpp::traits::input_parameter< int >::type xmin(xminSEXP);
    Rcpp::traits::input_parameter< int >::type ymin(yminSEXP);
    Rcpp::traits::input_parameter< double >::type k(kSEXP);
    Rcpp::traits::input_parameter< double >::type b(bSEXP);
    Rcpp::traits::input_parameter< double >::type successionTimestep(successionTimestepSEXP);
    Rcpp::traits::input_parameter< double >::type verbose(verboseSEXP);
    Rcpp::traits::input_parameter< int >::type maxSpiral(maxSpiralSEXP);
    rcpp_result_gen = Rcpp::wrap(spiralSeedDispersal(cellCoords, speciesVectorsList, rcvSpeciesByIndex, speciesTable, numCols, numCells, cellSize, xmin, ymin, k, b, successionTimestep, verbose, maxSpiral));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_LandR_which2", (DL_FUNC) &_LandR_which2, 1},
    {"_LandR_rmElem", (DL_FUNC) &_LandR_rmElem, 2},
    {"_LandR_spiralSeedDispersal", (DL_FUNC) &_LandR_spiralSeedDispersal, 14},
    {NULL, NULL, 0}
};

RcppExport void R_init_LandR(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
