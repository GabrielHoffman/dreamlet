#include "beachmat3/beachmat.h"
#include <vector>
#include <algorithm>
#include <unordered_map>
#include <string>

#include  "Rcpp.h"

using namespace Rcpp;

// Not necessary in a package context:
// [[Rcpp::depends(beachmat)]]

// [[Rcpp::export(rng=false)]]
NumericMatrix colsum_beachmat_matrix(RObject mat, IntegerVector groupHsh, IntegerVector grpUnq) {

    // initialize beachmat interface
    auto ptr = beachmat::read_lin_block(mat);
    std::vector<double> workspace(ptr->get_nrow());

    // initialize output matrix
    NumericMatrix output(ptr->get_nrow(), grpUnq.size());

    // for each column i
    for (size_t i = 0; i < ptr->get_ncol(); ++i) {

        // pointer to column i
        auto colptr = ptr->get_col(i, workspace.data());
 
        // implicit version:
        /// initialize NumericVector with values, then assign to matrix
        NumericMatrix::Column col = output(_, groupHsh(i)-1);

        col = col + NumericVector::import(colptr, colptr + ptr->get_nrow());
    }

    return output;
}


// *manual* version with explicit loop
// assign values for each row r
// for(size_t r = 0; r < ptr->get_nrow(); ++r) {

//     // subtract 1 to convert from 1-based in R to 0-based in C++
//     output(r, groupHsh(i)-1) +=  colptr[r];
// }

// implicit version:
/// initialize NumericVector with values, then assign to matrix
// NumericVector tmp = NumericVector::import(colptr,colptr + ptr->get_nrow());

// Rcout << colptr[0] << ' ' << tmp(0) << std::endl;

// output(_, groupHsh(i)-1) = output(_, groupHsh(i)-1)+ NumericVector::import(colptr, colptr + ptr->get_nrow());


// [[Rcpp::export(rng=false)]]
NumericMatrix colsum_beachmat_sparseMatrix(RObject mat, IntegerVector groupHsh, IntegerVector grpUnq) {

    // initialize beachmat interface
    auto ptr = beachmat::read_lin_sparse_block(mat);
    std::vector<double> workspace_x(ptr->get_nrow());
    std::vector<int> workspace_i(ptr->get_nrow());

    // initialize output matrix
    NumericMatrix output(ptr->get_nrow(), grpUnq.size());

    // for each column i
    for (size_t i = 0; i < ptr->get_ncol(); ++i) {

        auto indices = ptr->get_col(i, workspace_x.data(), workspace_i.data());

        auto xptr = indices.x;
        auto iptr = indices.i; 
        auto nnzero = indices.n;

        // *manual* version with explicit loop
        // assign values for each row r
        for(size_t r = 0; r < nnzero; ++r) {

            // subtract 1 to convert from 1-based in R to 0-based in C++
            output(iptr[r], groupHsh(i)-1) +=  xptr[r];
        }
    }

    return output;
}







// output[i] = std::accumulate(colptr, colptr + ptr->get_nrow(), 0.0);

// std::vector<double> to_vect(ptr->get_nrow());
// std::copy(colptr, colptr + ptr->get_nrow(), to_vect);


// set row and column names of result
// only works with matrix, not sparseMatrix.  So just do in R
// Rcpp::CharacterVector myIntToStr(grpUnq.begin(), grpUnq.end());
// colnames(output) = myIntToStr;
// rownames(output) = rownames(mat);