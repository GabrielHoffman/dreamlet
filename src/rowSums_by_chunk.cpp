// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

// we only include RcppEigen.h which pulls Rcpp.h in for us
#include <RcppEigen.h>

// via the depends attribute we tell Rcpp to create hooks for
// RcppEigen so that the build process will know what to do
//
// [[Rcpp::depends(RcppEigen)]]

// [[Rcpp::depends(RcppProgress)]]
#include <progress.hpp>
#include <progress_bar.hpp>

// [[Rcpp::export]]
Rcpp::NumericMatrix rowSums_by_chunk_sparse(Eigen::MappedSparseMatrix<double> &data, Rcpp::List idxlst, bool verbose) { 

    // initialize NumericMatrix with zero values
    Rcpp::NumericMatrix result(data.rows(), idxlst.size()); 

    Progress progbar(idxlst.size(), verbose);

    // loop thru list
    for(int i=0; i < idxlst.size(); i++){

        Rcpp::NumericVector idx(idxlst[i]);

        // loop thru column indeces
        for(int j=0; j < idx.size(); j++){

            // loop thru genes (i.e. rows)
            for(int k=0; k< data.rows(); k++){

                // Note that these indecies are zero-based
                // Importantly, never create a subset of the matrix
                //  just extract one element at a time
                double value = (double) data.coeff(k,idx(j)-1);

                if( value != 0 ) result(k,i) += value;
            }
        }
        progbar.increment(); 
    }

    return result;
}


// [[Rcpp::export]]
Rcpp::NumericMatrix rowSums_by_chunk(Rcpp::NumericMatrix &data, Rcpp::List idxlst, bool verbose) { 

    // initialize NumericMatrix with zero values
    Rcpp::NumericMatrix result(data.rows(), idxlst.size()); 

    Progress progbar(idxlst.size(), verbose);

    // loop thru list
    for(int i=0; i < idxlst.size(); i++){

        Rcpp::NumericVector idx(idxlst[i]);

        // loop thru column indeces
        for(int j=0; j < idx.size(); j++){

            // loop thru genes (i.e. rows)
            for(int k=0; k< data.rows(); k++){

                // Note that these indecies are zero-based
                // Importantly, never create a subset of the matrix
                //  just extract one element at a time
                double value = (double) data(k,idx(j)-1);

                if( value != 0 ) result(k,i) += value;
            }
        }
        progbar.increment(); 
    }

    return result;
}