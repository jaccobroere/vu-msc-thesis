#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
List cpp_FSPLASH_precalc(const sp_mat& Vhat_d, const sp_mat& Dtilde_inv, const vec& sigma_hat, int m) {
  
  // Perform operations
  sp_mat XD1 = Vhat_d * Dtilde_inv;
  sp_mat X1 = XD1.cols(0, m - 1);
  mat X2 = mat(XD1.cols(m, XD1.n_cols - 1));
  mat X2_plus = solve(trans(X2) * X2, trans(X2)); 
  
  // Calculate I - P directly, P = X2 %*% X2_plus
  mat IminP = eye(X2.n_rows, X2.n_rows) - X2 * X2_plus; 
  vec ytilde = IminP * sigma_hat;
  mat Xtilde = IminP * X1;
  
  // Return a List
  return List::create(Named("ytilde") = ytilde,
                      Named("Xtilde") = Xtilde,
                      Named("X2_plus") = X2_plus,
                      Named("X1") = X1);
}
