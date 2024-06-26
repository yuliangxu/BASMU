# Generated by using Rcpp::compileAttributes() -> do not edit by hand
# Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#' @title Scalar on image regression with unmeasured confounder
#' @name BASMU_outcome
#'
#' @description 
#' {A basis decomposition is used. The main coefficient beta follows STGP prior.
#' Kernel matrices need to be pre-speficified}
#'
#' 
#' @param Y The scalar outcome variable, n by 1
#' @param M The image predictor, p by n
#' @param X The scalar exposure variable, n by 1
#' @param C The q confounders, n by q
#' @param L_all A vector of length num_region, each element is an integer to indicate the number of basis in each region
#' @param eta A matrix of size p by n, the latent confounder
#' @param delta_sp A p-by-1 indicator vector to specify the non-zero elements in nu, initial value for delta_nu
#' @param num_region An integer, the total number of regions
#' @param region_idx A list object of length num_region, each element is a vector of
#' the indices of each voxel in that region. Note that this index starts from 0.
#' @param n_mcmc An integer to indicate the total number of MCMC iterations
#' @param K A list object of length num_region, the r-th element is a p_r by L_r matrix for the basis function
#' @param stop_burnin An integer to indicate from which iteration to stop burnin period.
#' Note that during burinin, the step size in MALA is adjusted every interval_step iterations.
#' @param start_joint An integer to indicate from which iteration to start updating parameters other than beta. Default = 0.
#' @param lambda A numeric variable to indicate the thresholding parameter lambda in STGP prior. Default = 0.
#' @param target_accept_vec A vector of length num_region. Each element is a numeric variable in (0,1).
#' This allows the user to define different target acceptance rate for each region in the MALA algorithm,
#' and the step size will be adjusted to meet the target acceptance rate.
#' @param a A numeric variable for the Inverse-Gamma(a,b), priors for \eqn{\sigma^2_Y,\sigma^2_\beta}
#' @param b A numeric variable for the Inverse-Gamma(a,b), priors for \eqn{\sigma^2_Y,\sigma^2_\beta}
#' @param init A list object that contains the following element
#' \itemize{
#'   \item theta_beta A vector of length L. Initial value for theta_beta
#'   \item gamma A numeric scalar, initial value for gamma
#'   \item zetay A vector of length q, intial value for zetay
#'   \item sigma_Y A numeric scalar, initial value for sigma_Y
#'   \item sigma_beta A numeric scalar, initial value for sigma_beta
#'   \item nu A vector of length p, initial value for nu
#'   \item cb A numeric scalar, initial value for cy (=0)
#'   \item D A vector of length L. Eigenvalues for all regions in the basis
#' }
#' @param step A numeric vector of length num_region, the initial step size for each region
#' @param interval_step An integer to denote how often to update the step size
#' @param interval_thin An integer to denote how often to save the MCMC samples for theta_alpha
#' @import Rcpp
#' @useDynLib BASMU, .registration=TRUE
#' @export
#' @return A List object with the following component
#' \itemize{
#' \item theta_alpha_mcmc
#' \item logll_mcmc
#' \item track_step
#' \item accept_block
#' \item emp_accept
#' \item gs A list object with the following component
#'   \itemize{
#'     \item theta_zetam_mcmc
#'     \item sigma_M2_inv_mcmc
#'     \item sigma_alpha2_inv_mcmc
#'     \item sigma_eta2_inv_mcmc
#'   }
#' \item Timer
#' }
#' 
BASMU_outcome <- function(Y, M, X, C, L_all, eta, delta_sp, num_region, region_idx, n_mcmc, K, stop_burnin, target_accept_vec, a, b, init, step, interval_step, interval_thin, start_joint = 0L, lambda = 0, residual_update = FALSE) {
    .Call(`_BASMU_BASMU_outcome`, Y, M, X, C, L_all, eta, delta_sp, num_region, region_idx, n_mcmc, K, stop_burnin, target_accept_vec, a, b, init, step, interval_step, interval_thin, start_joint, lambda, residual_update)
}

#' @title Image on scalar regression
#' @name BASMU_mediator
#'
#' @description 
#' {A basis decomposition is used. The main coefficient alpha follows STGP prior.
#' Kernel matrices need to be pre-speficified}
#'
#' @param M The image predictor, p by n
#' @param X The scalar exposure variable, n by 1
#' @param C The q confounders, n by q
#' @param L_all A vector of length num_region, each element is an integer to indicate the number of basis in each region
#' @param num_region An integer, the total number of regions
#' @param region_idx A list object of length num_region, each element is a vector of
#' the indices of each voxel in that region. Note that this index starts from 0.
#' @param n_mcmc An integer to indicate the total number of MCMC iterations
#' @param K A list object of length num_region, the r-th element is a p_r by L_r matrix for the basis function
#' @param stop_burnin An integer to indicate from which iteration to stop burnin period.
#' Note that during burinin, the step size in MALA is adjusted every interval_step iterations.
#' @param lambda A numeric variable to indicate the thresholding parameter lambda in STGP prior
#' @param target_accept_vec A vector of length num_region. Each element is a numeric variable in (0,1).
#' This allows the user to define different target acceptance rate for each region in the MALA algorithm,
#' and the step size will be adjusted to meet the target acceptance rate.
#' @param a A numeric variable for the Inverse-Gamma(a,b), priors for \eqn{\sigma^2_Y,\sigma^2_\beta}
#' @param b A numeric variable for the Inverse-Gamma(a,b), priors for \eqn{\sigma^2_Y,\sigma^2_\beta}
#' @param init A list object that contains the following element
#' \itemize{
#'   \item theta_alpha A vector of length L. Initial value for theta_beta
#'   \item theta_zetam A L by q matrix. Initial value for theta_zetam
#'   \item theta_eta A matrix (L by n). Initial value for theta_eta
#'   \item D A vector of length L. Eigenvalues for all regions in the basis
#'   \item sigma_M A numeric scalar, initial value for gamma
#'   \item zetay A vector of length q, intial value for zetay
#'   \item sigma_alpha A numeric scalar, intial value for sigma_alpha
#'   \item sigma_eta A numeric scalar, initial value for sigma_eta
#' }
#' @param step A numeric vector of length num_region, the initial step size for each region
#' @param interval An integer to denote how often to update the step size
#' @param interval_eta An integer to denote how often to update theta_eta
#' @param thinning An integer to indicate how often to save the MCMC samples for theta_alpha
#' @param display_progress True for displaying progress bar
#' @import Rcpp
#' @useDynLib BASMU, .registration=TRUE
#' @export
#' @return A List object with the following component
#' \itemize{
#' \item theta_alpha_mcmc
#' \item logll_mcmc
#' \item track_step
#' \item accept_block
#' \item emp_accept
#' \item gs A list object with the following component
#'   \itemize{
#'     \item theta_zetam_mcmc
#'     \item sigma_M2_inv_mcmc
#'     \item sigma_alpha2_inv_mcmc
#'     \item sigma_eta2_inv_mcmc
#'   }
#' \item Timer
#' }
#' 
BASMU_mediator <- function(M, X, C, L_all, num_region, region_idx, n_mcmc, K, stop_burnin, lambda, target_accept_vec, init, interval, step, a = 1, b = 1, interval_eta = 10L, thinning = 10L, begin_save_eta = 1000L, begin_update_eta = 0L, eta_percentL = 1, only_update_eta = FALSE, update_eta_sgld = FALSE) {
    .Call(`_BASMU_BASMU_mediator`, M, X, C, L_all, num_region, region_idx, n_mcmc, K, stop_burnin, lambda, target_accept_vec, init, interval, step, a, b, interval_eta, thinning, begin_save_eta, begin_update_eta, eta_percentL, only_update_eta, update_eta_sgld)
}

complement <- function(start, end, n) {
    .Call(`_BASMU_complement`, start, end, n)
}

#' @title Get the null space of the design matrix G
#' @name get_H_mat
#' @useDynLib BASMU, .registration=TRUE
#' @export
get_H_mat <- function(G) {
    .Call(`_BASMU_get_H_mat`, G)
}

hyperplane_MVN_multiple <- function(G, H_mat, sigma2_vec, mu_mat) {
    .Call(`_BASMU_hyperplane_MVN_multiple`, G, H_mat, sigma2_vec, mu_mat)
}

rcpp_generate_beta <- function(n, alpha, beta) {
    .Call(`_BASMU_rcpp_generate_beta`, n, alpha, beta)
}

#' @title Map basis coefficients in matrix form to the image space
#' @name Low_to_high
#' @useDynLib BASMU, .registration=TRUE
#' @export
Low_to_high <- function(Low_mat, p, Phi_Q, region_idx, L_idx) {
    .Call(`_BASMU_Low_to_high`, Low_mat, p, Phi_Q, region_idx, L_idx)
}

#' @title Map basis coefficients to the image space
#' @name Low_to_high_vec
#' @useDynLib BASMU, .registration=TRUE
#' @export
Low_to_high_vec <- function(Low_vec, p, Phi_Q, region_idx, L_idx) {
    .Call(`_BASMU_Low_to_high_vec`, Low_vec, p, Phi_Q, region_idx, L_idx)
}

High_to_low_vec <- function(High_vec, L, Phi_Q, region_idx, L_idx) {
    .Call(`_BASMU_High_to_low_vec`, High_vec, L, Phi_Q, region_idx, L_idx)
}

High_to_low <- function(High_mat, L, Phi_Q, region_idx, L_idx) {
    .Call(`_BASMU_High_to_low`, High_mat, L, Phi_Q, region_idx, L_idx)
}

