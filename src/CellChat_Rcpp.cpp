#include <RcppEigen.h>
#include <Rcpp.h>

using namespace Rcpp;

typedef Eigen::Triplet<double> T;
// Adapted from swne (https://github.com/yanwu2014/swne)
//[[Rcpp::export]]
Eigen::SparseMatrix<double> ComputeSNN(Eigen::MatrixXd nn_ranked, double prune) {
  std::vector<T> tripletList;
  int k = nn_ranked.cols();
  tripletList.reserve(nn_ranked.rows() * nn_ranked.cols());
  for(int j=0; j<nn_ranked.cols(); ++j){
    for(int i=0; i<nn_ranked.rows(); ++i) {
      tripletList.push_back(T(i, nn_ranked(i, j) - 1, 1));
    }
  }
  Eigen::SparseMatrix<double> SNN(nn_ranked.rows(), nn_ranked.rows());
  SNN.setFromTriplets(tripletList.begin(), tripletList.end());
  SNN = SNN * (SNN.transpose());
  for (int i=0; i < SNN.outerSize(); ++i){
    for (Eigen::SparseMatrix<double>::InnerIterator it(SNN, i); it; ++it){
      it.valueRef() = it.value()/(k + (k - it.value()));
      if(it.value() < prune){
        it.valueRef() = 0;
      }
    }
  }
  SNN.prune(0.0); // actually remove pruned values
  return SNN;
}
