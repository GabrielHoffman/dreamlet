// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppEigen.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// rowSums_by_chunk_sparse
Rcpp::NumericMatrix rowSums_by_chunk_sparse(Eigen::MappedSparseMatrix<double>& data, Rcpp::List idxlst, bool verbose);
RcppExport SEXP _dreamlet_rowSums_by_chunk_sparse(SEXP dataSEXP, SEXP idxlstSEXP, SEXP verboseSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Eigen::MappedSparseMatrix<double>& >::type data(dataSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type idxlst(idxlstSEXP);
    Rcpp::traits::input_parameter< bool >::type verbose(verboseSEXP);
    rcpp_result_gen = Rcpp::wrap(rowSums_by_chunk_sparse(data, idxlst, verbose));
    return rcpp_result_gen;
END_RCPP
}
// rowSums_by_chunk
Rcpp::NumericMatrix rowSums_by_chunk(Rcpp::NumericMatrix& data, Rcpp::List idxlst, bool verbose);
RcppExport SEXP _dreamlet_rowSums_by_chunk(SEXP dataSEXP, SEXP idxlstSEXP, SEXP verboseSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix& >::type data(dataSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type idxlst(idxlstSEXP);
    Rcpp::traits::input_parameter< bool >::type verbose(verboseSEXP);
    rcpp_result_gen = Rcpp::wrap(rowSums_by_chunk(data, idxlst, verbose));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_dreamlet_rowSums_by_chunk_sparse", (DL_FUNC) &_dreamlet_rowSums_by_chunk_sparse, 3},
    {"_dreamlet_rowSums_by_chunk", (DL_FUNC) &_dreamlet_rowSums_by_chunk, 3},
    {NULL, NULL, 0}
};

RcppExport void R_init_dreamlet(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}