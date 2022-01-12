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

typedef Eigen::MappedSparseMatrix<double> MSpMat;
typedef MSpMat::InnerIterator InIterMat;
typedef Eigen::SparseMatrix<double> SpMat;

// [[Rcpp::export]]
Rcpp::NumericMatrix rowSums_by_chunk_sparse(Eigen::MappedSparseMatrix<double> &data, Rcpp::List idxlst, bool verbose) { 

    // initialize NumericMatrix with zero values
    Rcpp::NumericMatrix result(data.rows(), idxlst.size()); 
    Eigen::SparseMatrix<double> res(data.rows(), idxlst.size());

    Progress progbar(idxlst.size(), verbose);

    // loop thru list
    for(int i=0; i < idxlst.size(); i++){

        Rcpp::NumericVector idx(idxlst[i]);

        // loop thru column indeces
        for(int j=0; j < idx.size(); j++){

            // loop thru genes (i.e. rows)
            // Iteratator only accesses non-zero elements, and is much faster
            for (InIterMat g_(data, idx(j)-1); g_; ++g_){
                result(g_.index(),i) += g_.value();
            }

            // old version access every element
            // for(int k=0; k< data.rows(); k++){

            //     // Note that these indecies are zero-based
            //     // Importantly, never create a subset of the matrix
            //     //  just extract one element at a time
            //     double value = (double) data.coeff(k,idx(j)-1);

            //     if( value != 0 ) result(k,i) += value;
            // }
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

// // [[Rcpp::export]]
// Eigen::SparseMatrix<double> aggregateByColnames(Rcpp::List resList, Rcpp::List idLst, Rcpp::StringVector grpUniq) { 

//     MSpMat spM_tmp = resList(0);

//     int ngenes = spM_tmp.rows();

//     // initialize return value
//     SpMat spMatFinal(ngenes, grpUniq.size());
   
//     // for each group
//     #if defined(_OPENMP)
//     #pragma omp parallel for num_threads(omp_get_thread_num())
//     #endif
//     for(int i=0; i<grpUniq.size(); i++){
        
//         if(i % 100 == 0) Rcpp::Rcout << i << std::endl;                

//         Rcpp::NumericVector grpResult(ngenes);

//          // for each batch of columns
//         for(int j=0; j<resList.size(); j++){
//             MSpMat spM = resList(j);
//             Rcpp::StringVector colNames = idLst(j);

//             // for each column in spM
//             for(int h=0; h<colNames.size(); h++){
//                 if( colNames(h) == grpUniq(i)){
//                     // loop thru genes (i.e. rows)
//                     for (InIterMat g_(spM, h); g_; ++g_){
//                         grpResult(g_.index()) += g_.value();
//                     }
//                     break;
//                 }                    
//             }
//         }
//         double value;
//         #pragma omp critical
//         for(int k=0; k<grpResult.size(); k++){
//             value = grpResult(k);
//             if( value != 0) spMatFinal.insert(k,i) = value;
//         }
//     }

//     Rcpp::Rcout << "makeCompressed" << std::endl; 

//     spMatFinal.makeCompressed();
//     return spMatFinal;
// }

typedef Eigen::Triplet<double> T;

// [[Rcpp::export]]
Eigen::SparseMatrix<double> aggregateByColnames(Rcpp::List resList, Rcpp::List idLst, Rcpp::StringVector grpUniq) { 

    int ngenes;
    std::vector<T> tripletList;
   
    // for each batch of columns
    for(int j=0; j<resList.size(); j++){
        MSpMat spM = resList(j);
        Rcpp::StringVector colNames = idLst(j);
        if( j==0 ) ngenes = spM.rows();

        // for each group
        for(int i=0; i<grpUniq.size(); i++){
            double value;
            // for each column in spM
            for(int h=0; h<colNames.size(); h++){
                if( colNames(h) == grpUniq(i)){
                    // loop thru genes (i.e. rows)
                    for (InIterMat g_(spM, h); g_; ++g_){
                        value = g_.value();
                        if( value != 0){
                            // store i,j,value triplet for sparseMatrix
                            tripletList.push_back(T(g_.index(),i, value) );
                            //Rcpp::Rcout << g_.index() << " " << i << " " << value << std::endl;
                        }
                    }
                }                    
            }
        }
    }

    Rcpp::Rcout << "Add triplets" << std::endl;

    // populate sparse matrix for return
    // values are added for repeated entries
    SpMat spMatFinal(ngenes, grpUniq.size());
    spMatFinal.setFromTriplets(tripletList.begin(), tripletList.end());

    Rcpp::Rcout << "makeCompressed" << std::endl;

    // convert format for return to R
    spMatFinal.makeCompressed();

    return spMatFinal;
}





