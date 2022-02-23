// // -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

// // we only include RcppEigen.h which pulls Rcpp.h in for us
// #include <RcppEigen.h>

// // via the depends attribute we tell Rcpp to create hooks for
// // RcppEigen so that the build process will know what to do
// //
// // [[Rcpp::depends(RcppEigen)]]

// // [[Rcpp::depends(RcppProgress)]]
// #include <progress.hpp>
// #include <progress_bar.hpp>

// typedef Eigen::MappedSparseMatrix<double> MSpMat;
// typedef MSpMat::InnerIterator InIterMat;
// typedef Eigen::SparseMatrix<double> SpMat;

// typedef Eigen::Triplet<double> T;

// // [[Rcpp::export]]
// Eigen::SparseMatrix<double> cbind_list_of_sparseMatrix(Rcpp::List resList, bool verbose) { 

//     int ngenes;
//     std::vector<T> tripletList;
//     int ncols=0; 
//     double value;
//     int colidx = 0;

//     Progress progbar(resList.size(), verbose);

//     // for each batch of columns
//     for(int j=0; j<resList.size(); j++){
//         MSpMat spM = resList(j);
//         if( j==0 ) ngenes = spM.rows();
//         ncols += spM.cols();

//     	for(int i=0; i<spM.cols(); i++){
// 	        for (InIterMat g_(spM, i); g_; ++g_){
// 	            value = g_.value();
// 	            if( value != 0){
// 	                // store i,j,value triplet for sparseMatrix
// 	                tripletList.push_back(T(g_.index(),colidx+i, value) );
// 	            }
// 	        }
// 	    }
// 	    colidx += spM.cols();

//         progbar.increment(); 
//     }

//     // populate sparse matrix for return
//     // values are added for repeated entries
//     SpMat spMatFinal(ngenes, ncols);
//     spMatFinal.setFromTriplets(tripletList.begin(), tripletList.end());

//     // convert format for return to R
//     spMatFinal.makeCompressed();

//     return spMatFinal;
// }


