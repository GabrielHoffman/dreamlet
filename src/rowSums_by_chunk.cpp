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

#include <unordered_map>
#include <string>

typedef Eigen::MappedSparseMatrix<double> MSpMat;
typedef MSpMat::InnerIterator InIterMat;
typedef Eigen::SparseMatrix<double> SpMat;

// [[Rcpp::export]]
Rcpp::NumericMatrix rowSums_by_chunk_sparse(Eigen::MappedSparseMatrix<double> &data, Rcpp::List idxlst, bool verbose) { 

    // initialize NumericMatrix with zero values
    Rcpp::NumericMatrix result(data.rows(), idxlst.size()); 
    // Eigen::SparseMatrix<double> res(data.rows(), idxlst.size());

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

typedef Eigen::Triplet<double> T;

// [[Rcpp::export]]
Eigen::SparseMatrix<double> aggregateByColnames(Rcpp::List resList, Rcpp::List idLst, Rcpp::StringVector grpUniq) { 

    int ngenes, i;
    double value;
    std::vector<T> tripletList;
    int ncol = grpUniq.size();

    // create and populate hash table storing entries of grpUniq and index   
    std::unordered_map<std::string, int> grpUniqHash;

    for(int i=0; i<grpUniq.size(); i++){    
        grpUniqHash[Rcpp::as<std::string>(grpUniq(i))] = i;
    } 

    // for each batch of columns
    for(int j=0; j<resList.size(); j++){
        MSpMat spM = resList(j);
        Rcpp::StringVector colNames = idLst(j);
        if( j==0 ) ngenes = spM.rows();
           
        // for each column in spM
        for(int h=0; h<colNames.size(); h++){
            // look up location of colname in hash
            i = grpUniqHash[Rcpp::as<std::string>(colNames(h))];

            // loop thru genes (i.e. rows)
            for (InIterMat g_(spM, h); g_; ++g_){
                value = g_.value();
                if( value != 0){
                    // store i,j,value triplet for sparseMatrix
                    tripletList.push_back( T(g_.index(), i, value) );
                   // Rcpp::Rcout << g_.index() << " " << i << " " << value << std::endl;
                }
            }
        }
    }

    //Rcpp::Rcout << "Triplets: " << tripletList.size() << std::endl;
    //Rcpp::Rcout << \"Add triplets\" << std::endl;

    // populate sparse matrix for return
    // values are added for repeated entries
    SpMat spMatFinal(ngenes, ncol);
    spMatFinal.setFromTriplets(tripletList.begin(), tripletList.end());

    // clear here to free memory
    tripletList.clear();

    //Rcpp::Rcout << \"makeCompressed...\";

    // convert format for return to R
    spMatFinal.makeCompressed();
    // Rcpp::Rcout << \"done\" << std::endl;

    return spMatFinal;
}


// [[Rcpp::export]]
Eigen::SparseMatrix<double> aggregateByColnames1(Eigen::MappedSparseMatrix<double> &spM, Rcpp::StringVector &colNames, Rcpp::StringVector &grpUniq) { 

    int i;
    double value;
    std::vector<T> tripletList;
    int ncol = grpUniq.size();
    int ngenes = spM.rows();

    // create and populate hash table storing entries of grpUniq and index   
    std::unordered_map<std::string, int> grpUniqHash;

    for(int i=0; i<grpUniq.size(); i++){    
        grpUniqHash[Rcpp::as<std::string>(grpUniq(i))] = i;
    } 
       
    // for each column in spM
    for(int h=0; h<colNames.size(); h++){
        // look up location of colname in hash
        i = grpUniqHash[Rcpp::as<std::string>(colNames(h))];

        // loop thru genes (i.e. rows)
        for (InIterMat g_(spM, h); g_; ++g_){
            value = g_.value();
            if( value != 0){
                // store i,j,value triplet for sparseMatrix
                tripletList.push_back( T(g_.index(), i, value) );
            }
        }
    }

    // populate sparse matrix for return
    // values are added for repeated entries
    SpMat spMatFinal(ngenes, ncol);
    spMatFinal.setFromTriplets(tripletList.begin(), tripletList.end());

    // clear here to free memory
    tripletList.clear();

    // convert format for return to R
    spMatFinal.makeCompressed();

    return spMatFinal;
}



// sum list of sparse matrices
// [[Rcpp::export]]
Eigen::SparseMatrix<double> sumSpMatList(Rcpp::List resList) { 
    
    MSpMat tmp = resList(0);
    int nrow = tmp.rows();
    int ncol = tmp.cols();

    SpMat spMatFinal(nrow, ncol);

    for(int j=0; j<resList.size(); j++){
        MSpMat spM = resList(j);
        spMatFinal += spM;
    }

    return spMatFinal;
}



