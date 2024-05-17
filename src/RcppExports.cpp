// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// BASMU_outcome
Rcpp::List BASMU_outcome(arma::colvec& Y, arma::mat& M, arma::colvec& X, arma::mat& C, arma::uvec L_all, arma::mat& eta, arma::uvec& delta_sp, arma::uword num_region, Rcpp::List& region_idx, int n_mcmc, Rcpp::List& K, int stop_burnin, arma::colvec target_accept_vec, double a, double b, Rcpp::List& init, arma::colvec step, int interval_step, int interval_thin, int start_joint, double lambda, bool residual_update);
RcppExport SEXP _BASMU_BASMU_outcome(SEXP YSEXP, SEXP MSEXP, SEXP XSEXP, SEXP CSEXP, SEXP L_allSEXP, SEXP etaSEXP, SEXP delta_spSEXP, SEXP num_regionSEXP, SEXP region_idxSEXP, SEXP n_mcmcSEXP, SEXP KSEXP, SEXP stop_burninSEXP, SEXP target_accept_vecSEXP, SEXP aSEXP, SEXP bSEXP, SEXP initSEXP, SEXP stepSEXP, SEXP interval_stepSEXP, SEXP interval_thinSEXP, SEXP start_jointSEXP, SEXP lambdaSEXP, SEXP residual_updateSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::traits::input_parameter< arma::colvec& >::type Y(YSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type M(MSEXP);
    Rcpp::traits::input_parameter< arma::colvec& >::type X(XSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type C(CSEXP);
    Rcpp::traits::input_parameter< arma::uvec >::type L_all(L_allSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type eta(etaSEXP);
    Rcpp::traits::input_parameter< arma::uvec& >::type delta_sp(delta_spSEXP);
    Rcpp::traits::input_parameter< arma::uword >::type num_region(num_regionSEXP);
    Rcpp::traits::input_parameter< Rcpp::List& >::type region_idx(region_idxSEXP);
    Rcpp::traits::input_parameter< int >::type n_mcmc(n_mcmcSEXP);
    Rcpp::traits::input_parameter< Rcpp::List& >::type K(KSEXP);
    Rcpp::traits::input_parameter< int >::type stop_burnin(stop_burninSEXP);
    Rcpp::traits::input_parameter< arma::colvec >::type target_accept_vec(target_accept_vecSEXP);
    Rcpp::traits::input_parameter< double >::type a(aSEXP);
    Rcpp::traits::input_parameter< double >::type b(bSEXP);
    Rcpp::traits::input_parameter< Rcpp::List& >::type init(initSEXP);
    Rcpp::traits::input_parameter< arma::colvec >::type step(stepSEXP);
    Rcpp::traits::input_parameter< int >::type interval_step(interval_stepSEXP);
    Rcpp::traits::input_parameter< int >::type interval_thin(interval_thinSEXP);
    Rcpp::traits::input_parameter< int >::type start_joint(start_jointSEXP);
    Rcpp::traits::input_parameter< double >::type lambda(lambdaSEXP);
    Rcpp::traits::input_parameter< bool >::type residual_update(residual_updateSEXP);
    rcpp_result_gen = Rcpp::wrap(BASMU_outcome(Y, M, X, C, L_all, eta, delta_sp, num_region, region_idx, n_mcmc, K, stop_burnin, target_accept_vec, a, b, init, step, interval_step, interval_thin, start_joint, lambda, residual_update));
    return rcpp_result_gen;
END_RCPP
}
// BASMU_mediator
List BASMU_mediator(arma::mat& M, arma::colvec& X, arma::mat& C, arma::uvec L_all, arma::uword num_region, Rcpp::List& region_idx, int n_mcmc, Rcpp::List& K, int stop_burnin, double lambda, arma::colvec& target_accept_vec, Rcpp::List& init, int interval, arma::vec step, double a, double b, int interval_eta, int thinning, int begin_save_eta, int begin_update_eta, double eta_percentL, bool only_update_eta, bool update_eta_sgld);
RcppExport SEXP _BASMU_BASMU_mediator(SEXP MSEXP, SEXP XSEXP, SEXP CSEXP, SEXP L_allSEXP, SEXP num_regionSEXP, SEXP region_idxSEXP, SEXP n_mcmcSEXP, SEXP KSEXP, SEXP stop_burninSEXP, SEXP lambdaSEXP, SEXP target_accept_vecSEXP, SEXP initSEXP, SEXP intervalSEXP, SEXP stepSEXP, SEXP aSEXP, SEXP bSEXP, SEXP interval_etaSEXP, SEXP thinningSEXP, SEXP begin_save_etaSEXP, SEXP begin_update_etaSEXP, SEXP eta_percentLSEXP, SEXP only_update_etaSEXP, SEXP update_eta_sgldSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::traits::input_parameter< arma::mat& >::type M(MSEXP);
    Rcpp::traits::input_parameter< arma::colvec& >::type X(XSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type C(CSEXP);
    Rcpp::traits::input_parameter< arma::uvec >::type L_all(L_allSEXP);
    Rcpp::traits::input_parameter< arma::uword >::type num_region(num_regionSEXP);
    Rcpp::traits::input_parameter< Rcpp::List& >::type region_idx(region_idxSEXP);
    Rcpp::traits::input_parameter< int >::type n_mcmc(n_mcmcSEXP);
    Rcpp::traits::input_parameter< Rcpp::List& >::type K(KSEXP);
    Rcpp::traits::input_parameter< int >::type stop_burnin(stop_burninSEXP);
    Rcpp::traits::input_parameter< double >::type lambda(lambdaSEXP);
    Rcpp::traits::input_parameter< arma::colvec& >::type target_accept_vec(target_accept_vecSEXP);
    Rcpp::traits::input_parameter< Rcpp::List& >::type init(initSEXP);
    Rcpp::traits::input_parameter< int >::type interval(intervalSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type step(stepSEXP);
    Rcpp::traits::input_parameter< double >::type a(aSEXP);
    Rcpp::traits::input_parameter< double >::type b(bSEXP);
    Rcpp::traits::input_parameter< int >::type interval_eta(interval_etaSEXP);
    Rcpp::traits::input_parameter< int >::type thinning(thinningSEXP);
    Rcpp::traits::input_parameter< int >::type begin_save_eta(begin_save_etaSEXP);
    Rcpp::traits::input_parameter< int >::type begin_update_eta(begin_update_etaSEXP);
    Rcpp::traits::input_parameter< double >::type eta_percentL(eta_percentLSEXP);
    Rcpp::traits::input_parameter< bool >::type only_update_eta(only_update_etaSEXP);
    Rcpp::traits::input_parameter< bool >::type update_eta_sgld(update_eta_sgldSEXP);
    rcpp_result_gen = Rcpp::wrap(BASMU_mediator(M, X, C, L_all, num_region, region_idx, n_mcmc, K, stop_burnin, lambda, target_accept_vec, init, interval, step, a, b, interval_eta, thinning, begin_save_eta, begin_update_eta, eta_percentL, only_update_eta, update_eta_sgld));
    return rcpp_result_gen;
END_RCPP
}
// complement
arma::uvec complement(arma::uword start, arma::uword end, arma::uword n);
RcppExport SEXP _BASMU_complement(SEXP startSEXP, SEXP endSEXP, SEXP nSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::uword >::type start(startSEXP);
    Rcpp::traits::input_parameter< arma::uword >::type end(endSEXP);
    Rcpp::traits::input_parameter< arma::uword >::type n(nSEXP);
    rcpp_result_gen = Rcpp::wrap(complement(start, end, n));
    return rcpp_result_gen;
END_RCPP
}
// get_H_mat
List get_H_mat(const ::arma::mat G);
RcppExport SEXP _BASMU_get_H_mat(SEXP GSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::traits::input_parameter< const ::arma::mat >::type G(GSEXP);
    rcpp_result_gen = Rcpp::wrap(get_H_mat(G));
    return rcpp_result_gen;
END_RCPP
}
// hyperplane_MVN_multiple
arma::mat hyperplane_MVN_multiple(const ::arma::mat G, const ::List H_mat, const ::arma::vec sigma2_vec, const ::arma::mat mu_mat);
RcppExport SEXP _BASMU_hyperplane_MVN_multiple(SEXP GSEXP, SEXP H_matSEXP, SEXP sigma2_vecSEXP, SEXP mu_matSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::traits::input_parameter< const ::arma::mat >::type G(GSEXP);
    Rcpp::traits::input_parameter< const ::List >::type H_mat(H_matSEXP);
    Rcpp::traits::input_parameter< const ::arma::vec >::type sigma2_vec(sigma2_vecSEXP);
    Rcpp::traits::input_parameter< const ::arma::mat >::type mu_mat(mu_matSEXP);
    rcpp_result_gen = Rcpp::wrap(hyperplane_MVN_multiple(G, H_mat, sigma2_vec, mu_mat));
    return rcpp_result_gen;
END_RCPP
}
// rcpp_generate_beta
arma::vec rcpp_generate_beta(int n, double alpha, double beta);
RcppExport SEXP _BASMU_rcpp_generate_beta(SEXP nSEXP, SEXP alphaSEXP, SEXP betaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< double >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< double >::type beta(betaSEXP);
    rcpp_result_gen = Rcpp::wrap(rcpp_generate_beta(n, alpha, beta));
    return rcpp_result_gen;
END_RCPP
}
// Low_to_high
arma::mat Low_to_high(arma::mat& Low_mat, int p, Rcpp::List& Phi_Q, Rcpp::List& region_idx, Rcpp::List& L_idx);
RcppExport SEXP _BASMU_Low_to_high(SEXP Low_matSEXP, SEXP pSEXP, SEXP Phi_QSEXP, SEXP region_idxSEXP, SEXP L_idxSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat& >::type Low_mat(Low_matSEXP);
    Rcpp::traits::input_parameter< int >::type p(pSEXP);
    Rcpp::traits::input_parameter< Rcpp::List& >::type Phi_Q(Phi_QSEXP);
    Rcpp::traits::input_parameter< Rcpp::List& >::type region_idx(region_idxSEXP);
    Rcpp::traits::input_parameter< Rcpp::List& >::type L_idx(L_idxSEXP);
    rcpp_result_gen = Rcpp::wrap(Low_to_high(Low_mat, p, Phi_Q, region_idx, L_idx));
    return rcpp_result_gen;
END_RCPP
}
// Low_to_high_vec
arma::colvec Low_to_high_vec(const arma::colvec& Low_vec, int p, const Rcpp::List& Phi_Q, const Rcpp::List& region_idx, const Rcpp::List& L_idx);
RcppExport SEXP _BASMU_Low_to_high_vec(SEXP Low_vecSEXP, SEXP pSEXP, SEXP Phi_QSEXP, SEXP region_idxSEXP, SEXP L_idxSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::colvec& >::type Low_vec(Low_vecSEXP);
    Rcpp::traits::input_parameter< int >::type p(pSEXP);
    Rcpp::traits::input_parameter< const Rcpp::List& >::type Phi_Q(Phi_QSEXP);
    Rcpp::traits::input_parameter< const Rcpp::List& >::type region_idx(region_idxSEXP);
    Rcpp::traits::input_parameter< const Rcpp::List& >::type L_idx(L_idxSEXP);
    rcpp_result_gen = Rcpp::wrap(Low_to_high_vec(Low_vec, p, Phi_Q, region_idx, L_idx));
    return rcpp_result_gen;
END_RCPP
}
// High_to_low_vec
arma::vec High_to_low_vec(arma::vec& High_vec, int L, Rcpp::List& Phi_Q, Rcpp::List& region_idx, Rcpp::List& L_idx);
RcppExport SEXP _BASMU_High_to_low_vec(SEXP High_vecSEXP, SEXP LSEXP, SEXP Phi_QSEXP, SEXP region_idxSEXP, SEXP L_idxSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec& >::type High_vec(High_vecSEXP);
    Rcpp::traits::input_parameter< int >::type L(LSEXP);
    Rcpp::traits::input_parameter< Rcpp::List& >::type Phi_Q(Phi_QSEXP);
    Rcpp::traits::input_parameter< Rcpp::List& >::type region_idx(region_idxSEXP);
    Rcpp::traits::input_parameter< Rcpp::List& >::type L_idx(L_idxSEXP);
    rcpp_result_gen = Rcpp::wrap(High_to_low_vec(High_vec, L, Phi_Q, region_idx, L_idx));
    return rcpp_result_gen;
END_RCPP
}
// High_to_low
arma::mat High_to_low(const arma::mat& High_mat, int L, Rcpp::List& Phi_Q, Rcpp::List& region_idx, Rcpp::List& L_idx);
RcppExport SEXP _BASMU_High_to_low(SEXP High_matSEXP, SEXP LSEXP, SEXP Phi_QSEXP, SEXP region_idxSEXP, SEXP L_idxSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type High_mat(High_matSEXP);
    Rcpp::traits::input_parameter< int >::type L(LSEXP);
    Rcpp::traits::input_parameter< Rcpp::List& >::type Phi_Q(Phi_QSEXP);
    Rcpp::traits::input_parameter< Rcpp::List& >::type region_idx(region_idxSEXP);
    Rcpp::traits::input_parameter< Rcpp::List& >::type L_idx(L_idxSEXP);
    rcpp_result_gen = Rcpp::wrap(High_to_low(High_mat, L, Phi_Q, region_idx, L_idx));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_BASMU_BASMU_outcome", (DL_FUNC) &_BASMU_BASMU_outcome, 22},
    {"_BASMU_BASMU_mediator", (DL_FUNC) &_BASMU_BASMU_mediator, 23},
    {"_BASMU_complement", (DL_FUNC) &_BASMU_complement, 3},
    {"_BASMU_get_H_mat", (DL_FUNC) &_BASMU_get_H_mat, 1},
    {"_BASMU_hyperplane_MVN_multiple", (DL_FUNC) &_BASMU_hyperplane_MVN_multiple, 4},
    {"_BASMU_rcpp_generate_beta", (DL_FUNC) &_BASMU_rcpp_generate_beta, 3},
    {"_BASMU_Low_to_high", (DL_FUNC) &_BASMU_Low_to_high, 5},
    {"_BASMU_Low_to_high_vec", (DL_FUNC) &_BASMU_Low_to_high_vec, 5},
    {"_BASMU_High_to_low_vec", (DL_FUNC) &_BASMU_High_to_low_vec, 5},
    {"_BASMU_High_to_low", (DL_FUNC) &_BASMU_High_to_low, 5},
    {NULL, NULL, 0}
};

RcppExport void R_init_BASMU(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
