// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppEigen.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// cbind_list_of_sparseMatrix
Eigen::SparseMatrix<double> cbind_list_of_sparseMatrix(Rcpp::List resList, bool verbose);
RcppExport SEXP _dreamlet_cbind_list_of_sparseMatrix(SEXP resListSEXP, SEXP verboseSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::List >::type resList(resListSEXP);
    Rcpp::traits::input_parameter< bool >::type verbose(verboseSEXP);
    rcpp_result_gen = Rcpp::wrap(cbind_list_of_sparseMatrix(resList, verbose));
    return rcpp_result_gen;
END_RCPP
}
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
// aggregateByColnames
Eigen::SparseMatrix<double> aggregateByColnames(Rcpp::List resList, Rcpp::List idLst, Rcpp::StringVector grpUniq);
RcppExport SEXP _dreamlet_aggregateByColnames(SEXP resListSEXP, SEXP idLstSEXP, SEXP grpUniqSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::List >::type resList(resListSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type idLst(idLstSEXP);
    Rcpp::traits::input_parameter< Rcpp::StringVector >::type grpUniq(grpUniqSEXP);
    rcpp_result_gen = Rcpp::wrap(aggregateByColnames(resList, idLst, grpUniq));
    return rcpp_result_gen;
END_RCPP
}
// aggregateByColnames1
Eigen::SparseMatrix<double> aggregateByColnames1(Eigen::MappedSparseMatrix<double>& spM, Rcpp::StringVector& colNames, Rcpp::StringVector& grpUniq);
RcppExport SEXP _dreamlet_aggregateByColnames1(SEXP spMSEXP, SEXP colNamesSEXP, SEXP grpUniqSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Eigen::MappedSparseMatrix<double>& >::type spM(spMSEXP);
    Rcpp::traits::input_parameter< Rcpp::StringVector& >::type colNames(colNamesSEXP);
    Rcpp::traits::input_parameter< Rcpp::StringVector& >::type grpUniq(grpUniqSEXP);
    rcpp_result_gen = Rcpp::wrap(aggregateByColnames1(spM, colNames, grpUniq));
    return rcpp_result_gen;
END_RCPP
}
// sumSpMatList
Eigen::SparseMatrix<double> sumSpMatList(Rcpp::List resList, bool verbose);
RcppExport SEXP _dreamlet_sumSpMatList(SEXP resListSEXP, SEXP verboseSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::List >::type resList(resListSEXP);
    Rcpp::traits::input_parameter< bool >::type verbose(verboseSEXP);
    rcpp_result_gen = Rcpp::wrap(sumSpMatList(resList, verbose));
    return rcpp_result_gen;
END_RCPP
}
// colsum_beachmat_matrix
NumericMatrix colsum_beachmat_matrix(RObject mat, IntegerVector groupHsh, IntegerVector grpUnq);
RcppExport SEXP _dreamlet_colsum_beachmat_matrix(SEXP matSEXP, SEXP groupHshSEXP, SEXP grpUnqSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::traits::input_parameter< RObject >::type mat(matSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type groupHsh(groupHshSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type grpUnq(grpUnqSEXP);
    rcpp_result_gen = Rcpp::wrap(colsum_beachmat_matrix(mat, groupHsh, grpUnq));
    return rcpp_result_gen;
END_RCPP
}
// colsum_beachmat_sparseMatrix
NumericMatrix colsum_beachmat_sparseMatrix(RObject mat, IntegerVector groupHsh, IntegerVector grpUnq);
RcppExport SEXP _dreamlet_colsum_beachmat_sparseMatrix(SEXP matSEXP, SEXP groupHshSEXP, SEXP grpUnqSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::traits::input_parameter< RObject >::type mat(matSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type groupHsh(groupHshSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type grpUnq(grpUnqSEXP);
    rcpp_result_gen = Rcpp::wrap(colsum_beachmat_sparseMatrix(mat, groupHsh, grpUnq));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_dreamlet_cbind_list_of_sparseMatrix", (DL_FUNC) &_dreamlet_cbind_list_of_sparseMatrix, 2},
    {"_dreamlet_rowSums_by_chunk_sparse", (DL_FUNC) &_dreamlet_rowSums_by_chunk_sparse, 3},
    {"_dreamlet_rowSums_by_chunk", (DL_FUNC) &_dreamlet_rowSums_by_chunk, 3},
    {"_dreamlet_aggregateByColnames", (DL_FUNC) &_dreamlet_aggregateByColnames, 3},
    {"_dreamlet_aggregateByColnames1", (DL_FUNC) &_dreamlet_aggregateByColnames1, 3},
    {"_dreamlet_sumSpMatList", (DL_FUNC) &_dreamlet_sumSpMatList, 2},
    {"_dreamlet_colsum_beachmat_matrix", (DL_FUNC) &_dreamlet_colsum_beachmat_matrix, 3},
    {"_dreamlet_colsum_beachmat_sparseMatrix", (DL_FUNC) &_dreamlet_colsum_beachmat_sparseMatrix, 3},
    {NULL, NULL, 0}
};

RcppExport void R_init_dreamlet(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
