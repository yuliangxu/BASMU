#include <RcppArmadillo.h>
#include "BASMU_help_fun.h"
using namespace Rcpp;
// [[Rcpp::depends(RcppArmadillo)]]



// [[Rcpp::export]]
arma::uvec complement(arma::uword start, arma::uword end, arma::uword n) {
  arma::uvec y1 = arma::linspace<arma::uvec>(0, start-1, start);
  arma::uvec y2 = arma::linspace<arma::uvec>(end+1, n-1, n-1-end);
  arma::uvec y = arma::join_cols(y1,y2);
  return y;
}

double adjust_acceptance(double accept,double sgm,double target = 0.1){
  double y = 1. + 1000.*(accept-target)*(accept-target)*(accept-target);
  if (y < .9)
    y = .9;
  if (y > 1.1)
    y = 1.1;
  sgm *= y;
  return sgm;
}

//' @title Get the null space of the design matrix G
//' @name get_H_mat
//' @useDynLib BASMU, .registration=TRUE
//' @export
// [[Rcpp::export(rng = false)]]
List get_H_mat(const::arma::mat G){
  arma::mat G_null = null(G);
  arma::mat H_inv = join_vert(G_null.t(),G);
  arma::mat H = inv(H_inv);
  return List::create(Named("H") = H,
                      Named("H_inv") = H_inv,
                      Named("G_null_t") = G_null.t());
}
// [[Rcpp::export(rng = false)]]
arma::mat hyperplane_MVN_multiple(const::arma::mat G,
                                  const::List H_mat,
                                  const::arma::vec sigma2_vec,
                                  const::arma::mat mu_mat){
  // prepare H mat
  arma::mat H = H_mat["H"];
  arma::mat H_inv = H_mat["H_inv"];
  arma::mat G_null_t = H_mat["G_null_t"];
  arma::uword q = G.n_rows;
  arma::uword n = G.n_cols;
  arma::mat H1 = H.cols(0,n-q-1);
  arma::mat H2 = H.cols(n-q,n-1);
  
  // get mean and variance for z1
  arma::mat Lambda11 = H1.t()* H1;
  arma::mat Lambda12 = H1.t()* H2;
  arma::mat Lambda11_inv = inv(Lambda11);
  arma::mat Lambda11_inv_sqrt = sqrtmat_sympd(Lambda11_inv);
  
  arma::mat mu_z1 = ( G_null_t + Lambda11_inv * Lambda12 * G * H_inv ) * mu_mat; // n by m
  arma::uword m = sigma2_vec.n_elem;
  arma::mat z1 =   Lambda11_inv_sqrt *  arma::randn(n-q,m);
  z1.each_row() %= sqrt(sigma2_vec.t());
  z1 += mu_z1;
  arma::mat x = H1*z1;
  return x.t();
}

double square(double y){
  return y*y;
}


// [[Rcpp::export]]
arma::vec rcpp_generate_beta(int n, double alpha, double beta) {
  arma::vec x = arma::randg(n, arma::distr_param(alpha,1.0));
  arma::vec y = arma::randg(n, arma::distr_param(beta,1.0));
  return x/(x+y);
}

arma::uvec arma_setdiff(arma::uvec x, arma::uvec y){
  
  x = arma::unique(x);
  y = arma::unique(y);
  
  for (arma::uword j = 0; j < y.n_elem; j++) {
    arma::uvec q1 = arma::find(x == y[j]);
    if (!q1.empty()) {
      x.shed_row(q1(0));
    }
  }
  
  return x;
}

//' @title Map basis coefficients in matrix form to the image space
//' @name Low_to_high
//' @useDynLib BASMU, .registration=TRUE
//' @export
// [[Rcpp::export]]
arma::mat Low_to_high(arma::mat& Low_mat, int p, Rcpp::List& Phi_Q,
                      Rcpp::List& region_idx, Rcpp::List& L_idx){
  int num_region = region_idx.size();
  int n = Low_mat.n_cols;
  arma::mat High_mat(p,n);
  for(arma::uword r=0; r<num_region; r++){
    arma::uvec p_idx = region_idx[r];
    arma::uvec L_range = L_idx[r];
    arma::mat Q = Phi_Q[r];
    High_mat.rows(p_idx) = Q*Low_mat.rows(L_range);
  }
  return High_mat;
}


//' @title Map basis coefficients to the image space
//' @name Low_to_high_vec
//' @useDynLib BASMU, .registration=TRUE
//' @export
// [[Rcpp::export]]
arma::colvec Low_to_high_vec(const arma::colvec& Low_vec, int p,
                             const Rcpp::List& Phi_Q,
                             const Rcpp::List& region_idx, 
                             const Rcpp::List& L_idx){
  int num_region = region_idx.size();
  arma::colvec High_vec(p,1);
  for(arma::uword r=0; r<num_region; r++){
    arma::uvec p_idx = region_idx[r];
    arma::uvec L_range = L_idx[r];
    arma::mat Q = Phi_Q[r];
    High_vec(p_idx) = Q*Low_vec(L_range);
  }
  return High_vec;
}
// [[Rcpp::export]]
arma::vec High_to_low_vec(arma::vec& High_vec, int L, Rcpp::List& Phi_Q,
                          Rcpp::List& region_idx, Rcpp::List& L_idx){
  int num_region = region_idx.size();
  arma::colvec Low_vec(L,1);
  for(arma::uword r=0; r<num_region; r++){
    arma::uvec p_idx = region_idx[r];
    arma::uvec L_range = L_idx[r];
    arma::mat Q = Phi_Q[r];
    Low_vec(L_range) = Q.t()*High_vec(p_idx);
  }
  return Low_vec;
  
}
// [[Rcpp::export]]
arma::mat High_to_low(const arma::mat& High_mat, int L, Rcpp::List& Phi_Q,
                      Rcpp::List& region_idx, Rcpp::List& L_idx){
  int num_region = region_idx.size();
  int n = High_mat.n_cols;
  arma::mat Low_mat = arma::zeros(L,n);
  for(arma::uword r=0; r<num_region; r++){
    arma::uvec p_idx = region_idx[r];
    arma::uvec L_range = L_idx[r];
    arma::mat Q = Phi_Q[r];
    Low_mat.rows(L_range) = Q.t()*High_mat.rows(p_idx);
  }
  return Low_mat;
  
}

