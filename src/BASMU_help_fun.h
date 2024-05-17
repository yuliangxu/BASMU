#ifndef BASMU_HELP_FUN_H
#define BASMU_HELP_FUN_H

#include <RcppArmadillo.h>

arma::uvec complement(arma::uword start, arma::uword end, arma::uword n);
double adjust_acceptance(double accept,double sgm,double target);
Rcpp::List get_H_mat(const arma::mat G);
arma::mat hyperplane_MVN_multiple(const arma::mat G,
                                  const Rcpp::List H_mat,
                                  const arma::vec sigma2_vec,
                                  const arma::mat mu_mat);

double square(double y);
arma::vec rcpp_generate_beta(int n, double alpha, double beta);

arma::uvec arma_setdiff(arma::uvec x, arma::uvec y);
arma::mat Low_to_high(arma::mat& Low_mat, int p, Rcpp::List& Phi_Q,
                      Rcpp::List& region_idx, Rcpp::List& L_idx);

arma::colvec Low_to_high_vec(const arma::colvec& Low_vec, int p,
                            const Rcpp::List& Phi_Q,
                            const Rcpp::List& region_idx, 
                            const Rcpp::List& L_idx);

arma::vec High_to_low_vec(arma::vec& High_vec, int L, Rcpp::List& Phi_Q,
                          Rcpp::List& region_idx, Rcpp::List& L_idx);

arma::mat High_to_low(const arma::mat& High_mat, int L, Rcpp::List& Phi_Q,
                      Rcpp::List& region_idx, Rcpp::List& L_idx);

#endif