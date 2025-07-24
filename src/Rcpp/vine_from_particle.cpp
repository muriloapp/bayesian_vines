// [[Rcpp::plugins(cpp17)]]
// [[Rcpp::plugins(latomic)]]                       // ⬅ we’ll define this below
// [[Rcpp::depends(RcppEigen, RcppThread, wdm, rvinecopulib)]]

#include <Rcpp.h>
#include <vinecopulib.hpp>          // <-- 0.6 layout uses this flat header
using namespace Rcpp;
using namespace vinecopulib;

// [[Rcpp::export]]
SEXP vine_from_particle_cpp(const NumericVector &theta,
                            const IntegerMatrix &struct_mat) {
  const int d   = struct_mat.nrow();
  auto pcs      = Vinecop::make_pair_copula_store(d);
  std::size_t k = 0;
  for (auto &tree : pcs)
    for (auto &edge : tree)
      edge = Bicop(BicopFamily::gaussian, 0,
                   Eigen::VectorXd::Constant(1, std::tanh(theta[k++])));
  Vinecop vc(std::move(pcs), as<Eigen::MatrixXi>(struct_mat));
  return wrap(vc);
}


