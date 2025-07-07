// [[Rcpp::depends(RcppThread, RcppEigen)]]
#include <RcppThread.h>
#include <RcppEigen.h>
using Eigen::MatrixXd;

// [[Rcpp::export]]
Rcpp::NumericVector grad_first_tree_gaussian(const MatrixXd& U,
                                             const Rcpp::NumericVector& theta)
{
  const int n  = U.rows();                      // block size
  const int d  = U.cols();
  const int L1 = d - 1;                         // #edges tree 1
  Rcpp::NumericVector g(L1);
  
  // Φ^{-1}(u) & Φ^{-1}(v) for the whole block
  MatrixXd Z = U.unaryExpr([](double u){return R::qnorm5(u,0,1,1,0);});
  
  // OpenMP loop over the L1 edges (i, j=i+1) of the first tree
#pragma omp parallel for schedule(static)
  for(int k = 0; k < L1; ++k){
    int i = 0;          // root node 1 is conditioned, leaves 2..d
    int j = k + 1;
    double theta_k = theta[k];
    double rho  = std::tanh(theta_k);
    double one_minus_r2 = 1. - rho*rho;
    double denom = one_minus_r2*one_minus_r2;
    double gk = 0.;
    
    for(int t = 0; t < n; ++t){
      double x = Z(t,i);
      double y = Z(t,j);
      gk += ( rho*(x*x + y*y) - (1.+rho*rho)*x*y ) / denom;
    }
    g[k] = (1.-rho*rho) * gk;     // chain rule dθ/dρ
  }
  return g;
}

