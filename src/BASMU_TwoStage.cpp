#include <RcppArmadillo.h>
#include "BASMU_help_fun.h"
#include <Rcpp/Benchmark/Timer.h>

using namespace Rcpp;
using namespace arma;


//' @title Scalar on image regression with unmeasured confounder
//' @name BASMU_outcome
//'
//' @description 
//' {A basis decomposition is used. The main coefficient beta follows STGP prior.
//' Kernel matrices need to be pre-speficified}
//'
//' 
//' @param Y The scalar outcome variable, n by 1
//' @param M The image predictor, p by n
//' @param X The scalar exposure variable, n by 1
//' @param C The q confounders, n by q
//' @param L_all A vector of length num_region, each element is an integer to indicate the number of basis in each region
//' @param eta A matrix of size p by n, the latent confounder
//' @param delta_sp A p-by-1 indicator vector to specify the non-zero elements in nu, initial value for delta_nu
//' @param num_region An integer, the total number of regions
//' @param region_idx A list object of length num_region, each element is a vector of
//' the indices of each voxel in that region. Note that this index starts from 0.
//' @param n_mcmc An integer to indicate the total number of MCMC iterations
//' @param K A list object of length num_region, the r-th element is a p_r by L_r matrix for the basis function
//' @param stop_burnin An integer to indicate from which iteration to stop burnin period.
//' Note that during burinin, the step size in MALA is adjusted every interval_step iterations.
//' @param start_joint An integer to indicate from which iteration to start updating parameters other than beta. Default = 0.
//' @param lambda A numeric variable to indicate the thresholding parameter lambda in STGP prior. Default = 0.
//' @param target_accept_vec A vector of length num_region. Each element is a numeric variable in (0,1).
//' This allows the user to define different target acceptance rate for each region in the MALA algorithm,
//' and the step size will be adjusted to meet the target acceptance rate.
//' @param a A numeric variable for the Inverse-Gamma(a,b), priors for \eqn{\sigma^2_Y,\sigma^2_\beta}
//' @param b A numeric variable for the Inverse-Gamma(a,b), priors for \eqn{\sigma^2_Y,\sigma^2_\beta}
//' @param init A list object that contains the following element
//' \itemize{
//'   \item theta_beta A vector of length L. Initial value for theta_beta
//'   \item gamma A numeric scalar, initial value for gamma
//'   \item zetay A vector of length q, intial value for zetay
//'   \item sigma_Y A numeric scalar, initial value for sigma_Y
//'   \item sigma_beta A numeric scalar, initial value for sigma_beta
//'   \item nu A vector of length p, initial value for nu
//'   \item cb A numeric scalar, initial value for cy (=0)
//'   \item D A vector of length L. Eigenvalues for all regions in the basis
//' }
//' @param step A numeric vector of length num_region, the initial step size for each region
//' @param interval_step An integer to denote how often to update the step size
//' @param interval_thin An integer to denote how often to save the MCMC samples for theta_alpha
//' @import Rcpp
//' @useDynLib BASMU, .registration=TRUE
//' @export
//' @return A List object with the following component
//' \itemize{
//' \item theta_alpha_mcmc
//' \item logll_mcmc
//' \item track_step
//' \item accept_block
//' \item emp_accept
//' \item gs A list object with the following component
//'   \itemize{
//'     \item theta_zetam_mcmc
//'     \item sigma_M2_inv_mcmc
//'     \item sigma_alpha2_inv_mcmc
//'     \item sigma_eta2_inv_mcmc
//'   }
//' \item Timer
//' }
//' 
// [[Rcpp::export(rng = false)]]
Rcpp::List BASMU_outcome(arma::colvec& Y, arma::mat& M,
                     arma::colvec& X, arma::mat& C,arma::uvec L_all,
                     arma::mat& eta, 
                     arma::uvec& delta_sp,
                     arma::uword num_region, Rcpp::List& region_idx,
                     int n_mcmc, Rcpp::List& K, int stop_burnin,
                     arma::colvec target_accept_vec,
                     double a, double b,
                     Rcpp::List& init,
                     arma::colvec step,
                     int interval_step,
                     int interval_thin,
                     int start_joint =0,
                     double lambda =0 ,
                     bool residual_update = false){
  Rcpp::Timer timer;
  timer.step("start of precomputation");
  arma::uword p = M.n_rows;
  arma::uword n = M.n_cols;
  arma::uvec L_cumsum = cumsum(L_all);
  arma::uword L_max = sum(L_all);
  
  arma::mat eta_t = eta.t();
  arma::uvec delta_nu = delta_sp;
  

  
  // input
  arma::colvec theta_beta = init["theta_beta"];
  arma::colvec D = init["D"];
  arma::colvec D_sqrt = sqrt(D);
  double gamma = init["gamma"], cy = init["cb"];
  arma::colvec zetay = init["zetay"];
  double sigma_Y = init["sigma_Y"], sigma_Y2 = sigma_Y*sigma_Y;
  double sigma_beta = init["sigma_beta"], sigma_beta2 = sigma_beta*sigma_beta;
  double sigma_nu2 = sigma_beta2;
  arma::colvec step_all;
  if(step.n_elem==1){
    step_all = step(0)*arma::ones(num_region);
  }else{
    step_all = step;
  }
  arma::colvec step_all_nu = step_all;
  // hyper parameter for inverse-Gamma
  double sigma_gamma2 = 1.0, sigma_zeta_y2 = 1.0,sigma_cy2=1.0;
  if(C.n_cols !=zetay.n_rows){
    Rcout<<"Error: dimensions of C and zetay don't match!"<<
      "dim of C = "<<size(C)<<"; dim of zetay = "<<size(zetay)<<std::endl;
    return Rcpp::List::create(Rcpp::Named("ERROR")=1);
  }
  
  arma::mat M_t = M.t(); arma::mat C_t = C.t();
  
  Rcpp::List Q_t(num_region) ;
  arma::colvec beta = arma::zeros(p,1);
  
  arma::colvec nu = init["nu"];
  
  
  for(int l=0; l<num_region; l++){
    arma::uvec idx = region_idx[l];
    arma::mat Q = K[l];
    Q_t[l] = Q.t();
    arma::uvec L_idx;
    if(l==0){
      L_idx = arma::linspace<arma::uvec>(0,L_cumsum(l)-1,L_all(l));
    }else{
      L_idx = arma::linspace<arma::uvec>(L_cumsum(l-1),L_cumsum(l)-1,L_all(l));
    }
    beta(idx) = Q*theta_beta(L_idx);
  }
  arma::colvec Y_star = Y - cy - gamma*X -  C*zetay ;
  
  
  arma::mat svd_m_U, svd_m_V;
  arma::vec svd_m_D;
  
  
  // svd on C_t
  arma::mat svd_C_U, svd_C_V;
  arma::vec svd_C_D;
  arma::svd_econ(svd_C_U,svd_C_D,svd_C_V,C,"both","std");

  
  //return
  arma::mat theta_beta_mcmc_thin = arma::zeros(L_max,n_mcmc/interval_thin);
  arma::mat nu_mcmc_thin = arma::zeros(p,n_mcmc/interval_thin);
  
  arma::colvec logll_mcmc_Y = arma::zeros(n_mcmc);
  arma::mat track_step = arma::zeros(num_region,n_mcmc);
  arma::mat emp_accept = arma::zeros( n_mcmc/interval_step,num_region);
  arma::mat accept_block = arma::zeros(n_mcmc,num_region);
  arma::cube time_seg = arma::zeros(n_mcmc,num_region,5);

  // gs returns
  arma::colvec gamma_mcmc = arma::zeros(n_mcmc,1);
  arma::mat zetay_mcmc = arma::zeros(zetay.n_elem,n_mcmc);
  arma::colvec cy_mcmc = arma::zeros(n_mcmc,1);
  arma::colvec sigma_beta2_mcmc = arma::zeros(n_mcmc,1);
  arma::colvec sigma_nu2_mcmc = arma::zeros(n_mcmc,1);
  arma::colvec sigma_Y2_mcmc = arma::zeros(n_mcmc,1);
  
  // check if region_idx starts from 0!
  arma::uword min_region_idx=1;
  for(arma::uword l =0; l<num_region;l++){
    arma::uvec region_idx0 = region_idx[l];
    arma::uword min_l = min(region_idx0);
    if(min_l<min_region_idx){min_region_idx=min_l;}
  }
  
  if(min_region_idx>0){
    Rcout<<"Error: region_idx does not start from 0!"<<
      "min(region_idx[0]) = "<<min_region_idx<<std::endl;
    return Rcpp::List::create(Rcpp::Named("ERROR")=1,
                              Rcpp::Named("region_idx")=region_idx);
  }
  
  
  arma::uword all_iter=0;
  arma::uword num_block = num_region;
  arma::uvec delta_in_block, delta_in_block_nu;
  timer.step("start of iteration");
  for(int iter=0; iter<n_mcmc; iter++){
    
    if(iter==stop_burnin){
      timer.step("stop of burnin");        // record the starting point
    }
    Y_star = Y- cy - gamma*X -  C*zetay;
    
    
    if(!residual_update){
      // start block update
      double log_target_density;
      double log_target_density_new;
      arma::colvec tb_temp;
      arma::colvec temp_new;
      arma::colvec temp;
      double temp_2norm2;
      
      double log_target_density_nu;
      double log_target_density_new_nu;
      arma::colvec temp_nu;
      arma::colvec tb_temp_nu;
      arma::colvec temp_new_nu;
      double temp_2norm2_nu;
      
      for(arma::uword m=0; m < num_region; m++){

        // ------------- update theta_beta ------------- //
        // ------------- time_seg0 ------------- //
        clock_t t0;t0 = clock();
        
        arma::uvec delta = find(abs(beta)>lambda);
        arma::uvec idx = region_idx[m];
        arma::uvec delta_Q = find(abs(beta(idx))>lambda);
        arma::mat Q = K[m];
        arma::uvec L_idx;
        if(m==0){
          L_idx = arma::linspace<arma::uvec>(0,L_cumsum(m)-1,L_all(m));
        }else{
          L_idx = arma::linspace<arma::uvec>(L_cumsum(m-1),L_cumsum(m)-1,L_all(m));
        }
        
        t0 = clock() - t0;
        double sys_t0 = ((double)t0)/CLOCKS_PER_SEC;
        time_seg(iter,m,0) = sys_t0;
        // ------------- time_seg1 ------------- //
        clock_t t1;t1 = clock();
        
        arma::mat K_block_t = Q_t[m];
        delta_in_block = intersect(idx,delta);
        arma::colvec temp_idx = M_t.cols(delta_in_block)*(beta.rows(delta_in_block)-sign(beta.rows(delta_in_block))*lambda)/1.00;
        
        
        t1 = clock() - t1;
        double sys_t1 = ((double)t1)/CLOCKS_PER_SEC;
        time_seg(iter,m,1) = sys_t1;
        // ------------- time_seg2 ------------- //
        clock_t t2;t2 = clock();
        if(m==0){
          temp = Y_star-M_t.cols(delta)*(beta.rows(delta)-sign(beta.rows(delta))*lambda)/1.00;
          temp -= eta_t * (nu%delta_nu);
          temp_2norm2 = dot(temp,temp);
        }
        
        arma::colvec grad_f = -theta_beta(L_idx)/D(L_idx)/sigma_beta2+
          K_block_t.cols(delta_Q)*M.rows(delta_in_block)*temp/1.00/sigma_Y2;
        double step = step_all(m);
        arma::colvec theta_beta_diff = step*grad_f+sqrt(2*step)*arma::randn(size(grad_f));
        arma::colvec theta_beta_new_block = theta_beta(L_idx)+theta_beta_diff;
        // grad_f_all(L_idx) = grad_f;
        t2 = clock() - t2;
        double sys_t2 = ((double)t2)/CLOCKS_PER_SEC;
        time_seg(iter,m,2) = sys_t2;
        
        
        
        
        
        // MH step
        // ------------- time_seg3 ------------- //
        clock_t t3;t3 = clock();
        tb_temp = theta_beta(L_idx)/D_sqrt(L_idx);
        log_target_density = -0.5*dot(tb_temp,tb_temp)/sigma_beta2-
          0.5*temp_2norm2/sigma_Y2;
        
        arma::colvec beta_new = beta;
        beta_new(idx) += Q*theta_beta_diff;
        arma::uvec delta_new = arma::find(abs(beta_new)>lambda);
        arma::uvec delta_Q_new = arma::find(abs(beta_new(idx))>lambda);
        arma::uvec delta_in_block_new = intersect(idx,delta_new);
        
        arma::colvec x1 = arma::zeros(n,1); arma::colvec x2 = arma::zeros(n,1);
        bool b1 = delta_in_block.n_elem>0;
        bool b2 = delta_in_block_new.n_elem>0;
        if(b1){
          x1 = M_t.cols(delta_in_block)*(beta.rows(delta_in_block)-sign(beta.rows(delta_in_block))*lambda)/1.00;
        }
        if(b2){
          x2 = M_t.cols(delta_in_block_new)*(beta_new.rows(delta_in_block_new)-sign(beta_new.rows(delta_in_block_new))*lambda)/1.00;
        }
        
        if(!b1 && b2){
          temp_new = temp-x2;
        }
        if(b1 && !b2){
          temp_new = temp+x1;
        }
        if(b1 && b2){
          temp_new = temp+x1-x2;
        }
        if(!b1 && !b2){
          temp_new = temp;
        }
        
        t3 = clock() - t3;
        double sys_t3 = ((double)t3)/CLOCKS_PER_SEC;
        time_seg(iter,m,3) = sys_t3;
        // ------------- time_seg4 ------------- //
        clock_t t4;t4 = clock();
        tb_temp = theta_beta_new_block/D_sqrt(L_idx);
        log_target_density_new = -0.5*dot(tb_temp,tb_temp)/sigma_beta2-
          0.5*dot(temp_new,temp_new)/sigma_Y2;
        arma::colvec grad_f_new = -theta_beta_new_block/D(L_idx)/sigma_beta2+
          K_block_t.cols(delta_Q_new)*M.rows(delta_in_block_new)*temp_new/1.00/sigma_Y2; // L x 1
        
        double log_q = -1/4/step * pow(norm(-theta_beta_diff-step*grad_f_new,2),2);//accu, change this part, %
        double log_q_new = -1/4/step * pow(norm(theta_beta_diff-step*grad_f,2),2);
        double rho = log_target_density_new + log_q - log_target_density - log_q_new;
        
        if(log(arma::randu())<=rho){
          theta_beta(L_idx) = theta_beta_new_block;
          beta = beta_new;
          temp = temp_new;
          log_target_density=log_target_density_new;
          accept_block(iter,m) = 1;
        }
        
        t4 = clock() - t4;
        double sys_t4 = ((double)t4)/CLOCKS_PER_SEC;
        time_seg(iter,m,4) = sys_t4;
        
        
        
        
        
        
      }// end of region update
      
      
      if( (iter%interval_step==0) & (iter>0)  ){
        arma::uvec u = arma::linspace<arma::uvec>(iter-interval_step,iter-1,interval_step);
        emp_accept.row(iter/interval_step-1) = mean(accept_block.rows(u),0);
        if(iter<stop_burnin){
          
          for(arma::uword l = 0; l<num_block; l++){
            step_all(l)  = adjust_acceptance(emp_accept(iter/interval_step-1,l),step_all(l),target_accept_vec(l));
            if(step_all(l)>1){step_all(l)=1;}
          }
          
        }
        
      }
      arma::uword back =  ceil( (n_mcmc - stop_burnin)/(n_mcmc/100.0) );
      if(iter==stop_burnin & stop_burnin>back){
        arma::uvec u = arma::linspace<arma::uvec>(iter-back,iter-1,back);
        step_all = exp(mean(log(track_step.cols(u)),1));
      }
      if( (iter%interval_step==0) & (iter>0)  ){
        arma::uvec u = arma::linspace<arma::uvec>(iter-interval_step,iter-1,interval_step);
        emp_accept.row(iter/interval_step-1) = mean(accept_block.rows(u),0);
      }
      
      track_step.col(iter) = step_all;
      all_iter = all_iter+1;
      
    }
    
    
    
    // -------------- Update GS for all other parameters --------------
    if(iter>start_joint){
      arma::uvec delta = find(abs(beta)>lambda);
      
      // arma::colvec M_beta_term = M_t.cols(delta)*(beta.rows(delta)-sign(beta.rows(delta))*lambda)/1.00;
      arma::colvec M_beta_term = M_t*beta;
      
      
      // ------------------------- update theta_nu using fixed delta_sp ------------------------- //
      
      // 1. nu
      
      arma::vec Y_res =  Y - cy - M_beta_term - C * zetay - gamma*X;
      
      // split into 1 and 1 set
      arma::uvec delta_nu1 = find(delta_nu != 0);
      arma::uvec delta_nu0 = find(delta_nu == 0);

      // update delta_nu1 indices
      // svd on eta_t with delta_nu1
      if(delta_nu1.n_elem>0){

        arma::mat eta_t1 = eta_t.cols(delta_nu1);
        arma::svd_econ(svd_m_U,svd_m_D,svd_m_V,eta_t1,"both","std");
        
        arma::vec svd_m_Y = svd_m_U.t()*Y_res;
        double tau2 = sigma_nu2/sigma_Y2;

        if(delta_nu1.n_elem>Y.n_elem){
          arma::mat svd_m_D_Vt = arma::diagmat(svd_m_D)*svd_m_V.t();
          arma::vec alpha1 = arma::randn(delta_nu1.n_elem) * sqrt(sigma_nu2);
          arma::vec alpha2 = arma::randn(n) * sqrt(sigma_Y2);
          arma::vec nu_new = alpha1 + tau2 * svd_m_V *
            diagmat(svd_m_D/(1+tau2*svd_m_D%svd_m_D)) *
            (svd_m_Y - svd_m_D_Vt*alpha1-alpha2);
          nu(delta_nu1) = nu_new;
        }else{
          arma::vec nu_var = 1/(svd_m_D%svd_m_D/sigma_Y2 + 1/sigma_nu2);
          arma::vec nu_mean = nu_var % (svd_m_D%svd_m_Y/sigma_Y2);
          arma::vec nu_new = arma::randn(delta_nu1.n_elem) % sqrt(nu_var) + nu_mean;
          nu_new = svd_m_V * nu_new;
          nu(delta_nu1) = nu_new;

          // return Rcpp::List::create(Rcpp::Named("ERROR")=1,
          //                         Rcpp::Named("eta_t1") = eta_t1,
          //                         Rcpp::Named("svd_m_D")=svd_m_D,
          //                         Rcpp::Named("svd_m_V")=svd_m_V,
          //                         Rcpp::Named("svd_m_U")=svd_m_U,
          //                         Rcpp::Named("delta_nu1")=delta_nu1,
          //                         Rcpp::Named("delta_nu0")=delta_nu0,
          //                         Rcpp::Named("sigma_nu2")=sigma_nu2,
          //                         Rcpp::Named("sigma_Y2")=sigma_Y2,
          //                         Rcpp::Named("Y")=Y,
          //                         Rcpp::Named("M_beta_term")=M_beta_term,
          //                         Rcpp::Named("beta")=beta,
          //                         Rcpp::Named("Y_res")=Y_res,
          //                         Rcpp::Named("nu_mean") = svd_m_V*nu_mean,
          //                         Rcpp::Named("nu_new") = nu_new,
          //                         Rcpp::Named("nu") = nu);
        }
        
      }
      

      // update on delta_nu0 indices
      if(delta_nu0.n_elem>0){
        arma::vec nu0 = arma::randn(delta_nu0.n_elem) * sqrt(sigma_nu2);
        nu(delta_nu0) = nu0;
      }

      if(nu.has_nan() ){
        arma::mat eta_t1 = eta_t.cols(delta_nu1);
        Rcout<<"iter="<<iter<<"; nu has nan! or testing"<<std::endl;
        return Rcpp::List::create(Rcpp::Named("ERROR")=1,
                                  Rcpp::Named("eta_t1") = eta_t1,
                                  Rcpp::Named("svd_m_D")=svd_m_D,
                                  Rcpp::Named("svd_m_V")=svd_m_V,
                                  Rcpp::Named("svd_m_U")=svd_m_U,
                                  Rcpp::Named("delta_nu1")=delta_nu1,
                                  Rcpp::Named("delta_nu0")=delta_nu0,
                                  Rcpp::Named("Y")=Y,
                                  Rcpp::Named("beta") = beta,
                                  Rcpp::Named("gamma")=gamma,
                                  Rcpp::Named("zetay")=zetay,
                                  Rcpp::Named("nu")=nu,
                                  Rcpp::Named("delta_nu")=delta_nu,
                                  Rcpp::Named("nu_mcmc_thin") = nu_mcmc_thin);
      }


      
      
      
      // update delta_nu
      arma::colvec delta_prior = arma::zeros(p,1);
      delta_prior += 0.5;
      
      arma::colvec logL_all = arma::zeros(p);
      arma::vec res_nu = Y_res - eta_t*(nu%delta_nu);
      for(arma::uword j=0;j< p;j++){
        arma::uword idx = j;

        arma::vec res1, res0;
        if(delta_nu(idx)==1){
          res1 = res_nu ;
          res0 = res_nu + eta_t.col(idx)*nu(idx);
        }else{
          res1 = res_nu - eta_t.col(idx)*nu(idx);
          res0 = res_nu ;
        }
         
         
        double logL1 = -0.5/sigma_Y2*dot(res1,res1);
        double logL0 = -0.5/sigma_Y2*dot(res0,res0);
        double p1 = exp(logL1 - logL0); p1*=delta_prior(idx)/(1-delta_prior(idx));
        double p0 = 1/(p1+1); p1=1-p0;
        
        
        
        if(arma::randu()<p1){
          delta_nu(idx) = 1;
          res_nu = res1;
        }else{
          delta_nu(idx) = 0;
          res_nu = res0;
        }
        logL_all(idx) = logL1;
        
        
      }

      // // update delta_prior, use beta(1,1) as a prior
      // delta_prior = rbeta(p,1,1);

      
      
      
      // // 
      arma::colvec eta_nu_term = eta_t*(nu%delta_nu);
      
      // 2. gamma
      double post_sigma_gamma2 = 1/(arma::sum(arma::dot(X,X))/sigma_Y2 + 1/sigma_gamma2);
      double temp2 = dot(X, Y - cy - M_beta_term - C * zetay - eta_nu_term);
      double mu_gamma = post_sigma_gamma2 * temp2/sigma_Y2;
      gamma = arma::randn() * sqrt(post_sigma_gamma2) + mu_gamma;
      gamma_mcmc(iter) = gamma;
      
      arma::colvec temp_cy = Y - M_beta_term - C * zetay - gamma*X - eta_nu_term;
      
      if(!residual_update){
        // // 3. zeta_y, use svd update
        arma::vec Y_res_zetay = Y - cy - M_beta_term - gamma*X - eta_nu_term;
        arma::vec Y_res_zetay_U = svd_C_U.t()*Y_res_zetay;
        arma::vec pos_Var = 1/(svd_C_D%svd_C_D/sigma_Y2 + 1/sigma_zeta_y2);
        arma::vec pos_mean = pos_Var % (svd_C_D%Y_res_zetay_U/sigma_Y2);
        arma::vec zetay = pos_mean + arma::randn(C.n_cols)%sqrt(pos_Var);
        zetay = svd_C_V * zetay;
        zetay_mcmc.col(iter) = zetay;
        
        // update sigma_zeta_y2
        double a_zeta_y = a + zetay.n_elem/2;
        double b_zeta_y = b + dot(zetay,zetay)/2;
        sigma_zeta_y2 = 1/arma::randg( arma::distr_param(a_zeta_y,1/b_zeta_y) );
        
        
        // //   5. sigma_beta
        double sigma_beta_a = a + theta_beta.n_elem/2;
        double sigma_beta_b = b + dot(theta_beta,theta_beta/D)/2;
        sigma_beta2 = 1/arma::randg( arma::distr_param(sigma_beta_a,1/sigma_beta_b) );
        sigma_beta2_mcmc(iter) = sigma_beta2;
        
      }
      // // sigma_nu
      double sigma_nu_a = a + nu.n_elem/2;
      double sigma_nu_b = b + dot(nu, nu)/2;
      sigma_nu2 = 1/arma::randg( arma::distr_param(sigma_nu_a,1/sigma_nu_b) );
      sigma_nu2_mcmc(iter) = sigma_nu2;
      
      // // 
      // // //   6. sigma_Y
      arma::colvec res = Y - M_beta_term - C * zetay - gamma*X - eta_nu_term - cy;
      double sigma_Y_b = b + dot(res,res)/2;
      sigma_Y2 = 1/arma::randg( arma::distr_param(a + Y.n_elem/2,1/sigma_Y_b) );
      sigma_Y2_mcmc(iter) = sigma_Y2;
      
      arma::uword five_percent = floor(n_mcmc/20);
      if(iter%five_percent ==0){
        Rcout<<"iter="<<iter<<"; sigma_Y2 = "<<sigma_Y2<<"; sigma_nu2="<<sigma_nu2<<std::endl;
      }
      // Rcout<<"iter="<<iter<<"; sigma_Y2 = "<<sigma_Y2<<"; sigma_nu2="<<sigma_nu2<<std::endl;
    }
    
    // --------------------- summarize return --------------------- //
    
    // arma::uvec delta_Y = find(abs(beta)>lambda);
    Y_star = Y- cy - gamma*X -  C*zetay;
    arma::colvec temp_Y = Y_star - M_t*beta;
    
    temp_Y -= eta_t*(nu%delta_nu);
    logll_mcmc_Y(iter) = -dot(temp_Y,temp_Y)/2/sigma_Y2;
    logll_mcmc_Y(iter) += -0.5*n*log(sigma_Y2);
    
    if( (iter%interval_thin==0) & (iter>0) ){
      theta_beta_mcmc_thin.col(iter/interval_thin-1)=theta_beta;
      nu_mcmc_thin.col(iter/interval_thin-1)=nu%delta_nu;
    }

    
  }
  timer.step("end of iterations");
  List gs = Rcpp::List::create(Rcpp::Named("gamma_mcmc")=gamma_mcmc,
                               Rcpp::Named("zetay_mcmc")= zetay_mcmc,
                               Rcpp::Named("cy_mcmc")= cy_mcmc,
                               Rcpp::Named("sigma_beta2_mcmc")= sigma_beta2_mcmc,
                               Rcpp::Named("sigma_nu2_mcmc")= sigma_nu2_mcmc,
                               Rcpp::Named("sigma_Y2_mcmc")= sigma_Y2_mcmc);
  
  return Rcpp::List::create(
    Rcpp::Named("theta_beta_mcmc_thin") = theta_beta_mcmc_thin,
    Rcpp::Named("nu_mcmc_thin") = nu_mcmc_thin,
    Rcpp::Named("logll_mcmc_Y") = logll_mcmc_Y,
    Rcpp::Named("track_step") = track_step,
    Rcpp::Named("accept_blcok") = accept_block,
    Rcpp::Named("emp_accept") = emp_accept,
    Rcpp::Named("gs") = gs,
    Rcpp::Named("time_seg") = time_seg,
    Rcpp::Named("Timer")=timer
  );
}



//' @title Image on scalar regression
//' @name BASMU_mediator
//'
//' @description 
//' {A basis decomposition is used. The main coefficient alpha follows STGP prior.
//' Kernel matrices need to be pre-speficified}
//'
//' @param M The image predictor, p by n
//' @param X The scalar exposure variable, n by 1
//' @param C The q confounders, n by q
//' @param L_all A vector of length num_region, each element is an integer to indicate the number of basis in each region
//' @param num_region An integer, the total number of regions
//' @param region_idx A list object of length num_region, each element is a vector of
//' the indices of each voxel in that region. Note that this index starts from 0.
//' @param n_mcmc An integer to indicate the total number of MCMC iterations
//' @param K A list object of length num_region, the r-th element is a p_r by L_r matrix for the basis function
//' @param stop_burnin An integer to indicate from which iteration to stop burnin period.
//' Note that during burinin, the step size in MALA is adjusted every interval_step iterations.
//' @param lambda A numeric variable to indicate the thresholding parameter lambda in STGP prior
//' @param target_accept_vec A vector of length num_region. Each element is a numeric variable in (0,1).
//' This allows the user to define different target acceptance rate for each region in the MALA algorithm,
//' and the step size will be adjusted to meet the target acceptance rate.
//' @param a A numeric variable for the Inverse-Gamma(a,b), priors for \eqn{\sigma^2_Y,\sigma^2_\beta}
//' @param b A numeric variable for the Inverse-Gamma(a,b), priors for \eqn{\sigma^2_Y,\sigma^2_\beta}
//' @param init A list object that contains the following element
//' \itemize{
//'   \item theta_alpha A vector of length L. Initial value for theta_beta
//'   \item theta_zetam A L by q matrix. Initial value for theta_zetam
//'   \item theta_eta A matrix (L by n). Initial value for theta_eta
//'   \item D A vector of length L. Eigenvalues for all regions in the basis
//'   \item sigma_M A numeric scalar, initial value for gamma
//'   \item zetay A vector of length q, intial value for zetay
//'   \item sigma_alpha A numeric scalar, intial value for sigma_alpha
//'   \item sigma_eta A numeric scalar, initial value for sigma_eta
//' }
//' @param step A numeric vector of length num_region, the initial step size for each region
//' @param interval An integer to denote how often to update the step size
//' @param interval_eta An integer to denote how often to update theta_eta
//' @param thinning An integer to indicate how often to save the MCMC samples for theta_alpha
//' @param display_progress True for displaying progress bar
//' @import Rcpp
//' @useDynLib BASMU, .registration=TRUE
//' @export
//' @return A List object with the following component
//' \itemize{
//' \item theta_alpha_mcmc
//' \item logll_mcmc
//' \item track_step
//' \item accept_block
//' \item emp_accept
//' \item gs A list object with the following component
//'   \itemize{
//'     \item theta_zetam_mcmc
//'     \item sigma_M2_inv_mcmc
//'     \item sigma_alpha2_inv_mcmc
//'     \item sigma_eta2_inv_mcmc
//'   }
//' \item Timer
//' }
//' 
 // [[Rcpp::export(rng = false)]]
 List BASMU_mediator(arma::mat& M,
                     arma::colvec& X, arma::mat& C,arma::uvec L_all,
                     arma::uword num_region, Rcpp::List& region_idx,
                     int n_mcmc, Rcpp::List& K, int stop_burnin,
                     double lambda, arma::colvec& target_accept_vec,
                     Rcpp::List& init,
                     int interval,arma::vec   step,
                     double a=1, double b=1,
                     int interval_eta = 10,
                     int thinning = 10,
                     int begin_save_eta = 1000,
                     int begin_update_eta = 0,
                     double eta_percentL = 1,
                     bool only_update_eta = false,
                     bool update_eta_sgld = false){
   Rcpp::Timer timer;
   // set_seed(1);
   timer.step("start of precomputation");
   arma::uword p = M.n_rows;
   arma::uword n = M.n_cols;
   
   
   
   double a_eta_step = 1; double b_eta_step = 1; double gamma_eta_step = 0.55;
   if(update_eta_sgld){
     Rcout<<"use SGLD to update eta, load a_eta_step, b_eta_step, gamma_eta_step from init"<<std::endl;
     a_eta_step = init["a_eta_step"];
     b_eta_step = init["b_eta_step"];
     gamma_eta_step = init["gamma_eta_step"];
   }
   
   arma::uvec L_cumsum = cumsum(L_all);
   arma::uword L_max = L_cumsum(num_region-1);
   arma::uvec eta_Lidx_update = arma::linspace<arma::uvec>( 0, floor(eta_percentL*L_max)-1,
                                                            floor(eta_percentL*L_max) );    
   // input
   arma::vec   theta_alpha = init["theta_alpha"];
   arma::mat theta_zetam =  init["theta_zetam"]; //L by q
   arma::vec   D = init["D"];
   arma::uword q = theta_zetam.n_cols;
   
   double sigma_M = init["sigma_M"], sigma_M2 = sigma_M*sigma_M, sigma_M2_inv = 1/sigma_M2;
   double sigma_alpha = init["sigma_alpha"], sigma_alpha2 = sigma_alpha*sigma_alpha, sigma_alpha2_inv = 1/sigma_alpha2;
   double sigma_eta = init["sigma_eta"], sigma_eta2_inv = 1/sigma_eta/sigma_eta;
   arma::mat theta_eta = init["theta_eta"];double sigma_zetam2_inv = 100;
   
   // prepare for constrained eta update
   arma::mat G = (join_horiz(X,C)).t();
   Rcpp::List H_mat = get_H_mat(G);
   
   arma::vec  step_all;
   if(step.n_elem==1){
     step_all = step(0)*arma::ones(num_region);
   }else{
     step_all = step;
   }
   
   arma::mat C_t = C.t();
   Rcpp::List M_star_pre_eta(num_region);
   arma::vec  alpha = arma::zeros(p,1);
   arma::mat zetam = arma::zeros(p,q);
   arma::vec   X2_sum_allsample_q = arma::zeros(q,1);
   arma::mat XcXq_sumsq = arma::zeros(q-1,q); // arma::sum_i C[-j]*C_j
   arma::vec  XXq_sumsq = arma::zeros(q,1); // arma::sum_i C*C_j
   for(int j=0; j<q; j++){
     X2_sum_allsample_q(j) += arma::sum(C_t.row(j) %C_t.row(j));
     // XqYstar_term_allsample.col(j) += Y_star_b * trans(X_q.row(j));
     XXq_sumsq(j) += accu(X.t() %C_t.row(j));
     arma::uvec c_j = complement(j, j, q);
     XcXq_sumsq.col(j) += C_t.rows(c_j) * trans(C_t.row(j));
   }
   
   for(int l=0; l<num_region; l++){
     arma::uvec idx = region_idx[l];
     arma::mat Q = K[l];
     // M_star_pre_eta[l] = Q_t*(M.rows(idx) - arma::ones(idx.n_elem,1)*zetam.t()*C);
     arma::uvec L_idx;
     if(l==0){
       L_idx = arma::linspace<arma::uvec>(0,L_cumsum(l)-1,L_all(l));
     }else{
       L_idx = arma::linspace<arma::uvec>(L_cumsum(l-1),L_cumsum(l)-1,L_all(l));
     }
     
     alpha(idx) = Q*theta_alpha(L_idx);
     zetam.rows(idx) = Q*theta_zetam.rows(L_idx);
     
     
   }
   
   
   //return
   int total_mcmc = n_mcmc/thinning;
   arma::mat  theta_alpha_mcmc = arma::zeros(L_max,total_mcmc);
   arma::vec   logll_mcmc = arma::zeros(total_mcmc);
   arma::mat track_step = arma::zeros(num_region,n_mcmc);
   arma::mat emp_accept = arma::zeros( n_mcmc/interval,num_region);
   arma::vec   accept = arma::zeros(n_mcmc*num_region);
   arma::mat accept_block = arma::zeros(n_mcmc,num_region);
   arma::uword all_iter=0;
   arma::uword num_block = num_region;
   arma::mat theta_eta_res(L_max,n);
   arma::mat M_reg_res(L_max,n);
   
   // GS: initialize mcmc sequences
   arma::cube theta_zetam_mcmc = arma::zeros(L_max,q,total_mcmc );
   arma::vec  sigma_M2_inv_mcmc = arma::zeros(total_mcmc,1);
   arma::vec  sigma_alpha2_inv_mcmc = arma::zeros(total_mcmc,1);
   arma::mat theta_eta_temp(L_max,n);
   arma::mat theta_eta_mean(L_max,n); int eta_counter = 0;
   arma::vec  sigma_eta2_inv_mcmc = arma::zeros(total_mcmc,1);
   timer.step("start of iteration");
   for(int iter=0; iter<n_mcmc; iter++){
     double logll_M = 0;
     // if(only_update_eta){
     //   prog.increment();
     // }
     
     if(!only_update_eta){
       if(iter==stop_burnin){
         // start the timer
         timer.step("stop of burnin");        // record the starting point
       }
       arma::vec   alpha_new = alpha;
       
       // Rcout<<"11111"<<std::endl;
       // check if region_idx starts from 0!
       // start block update
       for(arma::uword m=0; m < num_region; m++){
         // Rcout<<"iter="<<iter<<"; region="<<m<<std::endl;
         arma::uvec delta = find(abs(alpha)>lambda);
         arma::uvec idx = region_idx[m];
         arma::uvec delta_Q = find(abs(alpha(idx))>lambda);
         arma::mat Q = K[m];
         arma::uvec L_idx;
         if(m==0){
           L_idx = arma::linspace<arma::uvec>(0,L_cumsum(m)-1,L_all(m));
         }else{
           L_idx = arma::linspace<arma::uvec>(L_cumsum(m-1),L_cumsum(m)-1,L_all(m));
         }
         arma::uvec delta_in_block = intersect(idx,delta);
         
         arma::mat K_block_t = Q.t();
         // arma::mat M_star = K_block_t*(M.rows(idx) - ones(idx.n_elem,1)*zetam.t()*C_t) - theta_eta.rows(L_idx);
         arma::mat M_star = K_block_t*M.rows(idx) - theta_zetam.rows(L_idx) * C_t - theta_eta.rows(L_idx);
         arma::mat temp = M_star - K_block_t.cols(delta_Q)*(alpha.rows(delta_in_block)-sign(alpha.rows(delta_in_block))*lambda)*X.t();
         arma::mat temp_X = temp;
         temp_X.each_row()%= X.t();
         arma::vec  temp_sum = arma::sum(temp_X,1)*sigma_M2_inv;
         arma::vec  grad_f = -theta_alpha(L_idx)/D(L_idx)*sigma_alpha2_inv+
           K_block_t.cols(delta_Q)*Q.rows(delta_Q) *temp_sum; // L x 1, use smooth function in grad
         
         double step = step_all(m);
         arma::vec   theta_alpha_diff = step*grad_f+sqrt(2*step)*arma::randn(size(grad_f));
         arma::vec  theta_alpha_new_block = theta_alpha(L_idx) + theta_alpha_diff;
         
         
         // MH step
         double log_target_density = -0.5*square(norm(theta_alpha(L_idx)/sqrt(D(L_idx)),2))*sigma_alpha2_inv-
           0.5*square(arma::norm(temp,"fro"))*sigma_M2_inv;
         
         
         alpha_new = alpha;// change this
         alpha_new(idx) += Q*theta_alpha_diff;
         // arma::uvec delta_new = find(abs(alpha_new)>lambda);// find for one region
         arma::uvec delta_Q_new = find(abs(alpha_new(idx))>lambda);
         arma::uvec delta_in_block_new = idx(delta_Q_new);
         // arma::uvec delta_in_block_new = intersect(idx,delta_new);
         arma::vec   alpha_new_region = alpha_new.rows(delta_in_block_new);
         arma::mat temp_new = M_star - K_block_t.cols(delta_Q_new)*(alpha_new_region-sign(alpha_new_region)*lambda)*X.t();
         
         arma::mat temp_X_new = temp_new;
         temp_X_new.each_row()%= X.t();
         arma::vec  temp_sum_new = arma::sum(temp_X_new,1)*sigma_M2_inv;
         arma::vec  grad_f_new = -theta_alpha_new_block/D(L_idx)*sigma_alpha2_inv+
           K_block_t.cols(delta_Q_new)*Q.rows(delta_Q_new) *temp_sum_new; // L x 1, use smooth function in grad
         // Rcout<<"test 6-5"<<std::endl;
         double log_target_density_new = -0.5*square(norm( theta_alpha_new_block/sqrt(D(L_idx)),2))*sigma_alpha2_inv-
           0.5*square(arma::norm(temp_new,"fro"))*sigma_M2_inv;
         //         // Rcout<<"test 7"<<std::endl;
         double log_q = -1/4/step * square(norm(-theta_alpha_diff-step*grad_f_new,2));
         double log_q_new = -1/4/step * square(norm(theta_alpha_diff-step*grad_f,2));
         double rho = log_target_density_new + log_q - log_target_density - log_q_new;
         
         if(log(arma::randu())<=rho){
           theta_alpha(L_idx) = theta_alpha_new_block;
           alpha = alpha_new;
           accept(all_iter) = 1;
           accept_block(iter,m) = 1;
           temp = temp_new;
         }
         
         
         
         
         
       }// true end of block update
       
       if( (iter%interval==0) & (iter>0)  ){
         arma::uvec u = arma::linspace<arma::uvec>(iter-interval,iter-1,interval);
         emp_accept.row(iter/interval-1) = mean(accept_block.rows(u),0);
         if(iter<stop_burnin){
           arma::vec  sigma_t = sqrt(2*step_all);
           for(arma::uword l = 0; l<num_block; l++){
             sigma_t(l) = adjust_acceptance(emp_accept(iter/interval-1,l),sigma_t(l),
                     target_accept_vec(l));
             step_all(l) = sigma_t(l)*sigma_t(l)/2;
             if(step_all(l)>1){step_all(l)=1;}
           }
           
         }
       }
       //         // when stop burnin, choose the average of last few steps
       arma::uword back =  (n_mcmc - stop_burnin)/(n_mcmc/10) ;
       if(iter==stop_burnin & stop_burnin > back){
         
         arma::uvec u = arma::linspace<arma::uvec>(iter-back,iter-1,back);
         step_all = exp(mean(log(track_step.cols(u)),1));
       }
       
       if( (iter%interval==0) & (iter>0)  ){
         arma::uvec u = arma::linspace<arma::uvec>(iter-interval,iter-1,interval);
         emp_accept.row(iter/interval-1) = mean(accept_block.rows(u),0);
       }
       
       track_step.col(iter) = step_all;
       all_iter = all_iter+1;
       
     }
     
     
     
     //     // -------------- Update all other parameters using GS --------------
     arma::mat Mstar_alpha_term = arma::zeros(size(theta_eta));
     for(arma::uword m=0; m<num_region; m++){
       arma::uvec delta = find(abs(alpha)>lambda);
       arma::uvec idx = region_idx[m];
       arma::uvec delta_Q = find(abs(alpha(idx))>lambda);
       arma::uvec delta_in_block = intersect(idx,delta);
       arma::mat Q = K[m];
       arma::uvec L_idx;
       if(m==0){
         L_idx = arma::linspace<arma::uvec>(0,L_cumsum(m)-1,L_all(m));
       }else{
         L_idx = arma::linspace<arma::uvec>(L_cumsum(m-1),L_cumsum(m)-1,L_all(m));
       }
       arma::mat K_block_t = Q.t();
       Mstar_alpha_term.rows(L_idx) = K_block_t.cols(delta_Q)*(M.rows(delta_in_block)-( alpha.rows(delta_in_block)-sign(alpha.rows(delta_in_block))*lambda)*X.t() );
     }
     
     if(!only_update_eta){
       // // 2. zetam
       arma::mat zetam_res = Mstar_alpha_term - theta_eta;// L by n
       arma::mat mean_zetam = arma::mat(L_max,q);
       
       for(int j =0; j<q; j++){
         arma::vec  Sigma_zetam_j = 1/(1/D*sigma_zetam2_inv + sigma_M2_inv*X2_sum_allsample_q(j));
         arma::uvec c_j = complement(j, j, q);
         // change the following line when theta_eta needs to be updated
         arma::vec  mean_zetam_j = zetam_res*C.col(j) - theta_zetam.cols(c_j) * XcXq_sumsq.col(j);
         mean_zetam_j %= Sigma_zetam_j*sigma_M2_inv;
         mean_zetam.col(j) = mean_zetam_j;
         theta_zetam.col(j) = arma::randn(L_max,1)%sqrt(Sigma_zetam_j) +  mean_zetam_j;
       }
       // arma::mat Sigma_zetam_inv = sigma_zetam2_inv + sigma_M2_inv*q2_sum*C_t*C;
       // arma::mat Sigma_zetam = inv_sympd(Sigma_zetam_inv );
       // arma::mat zeta_res = Mstar_alpha_term - theta_eta;
       // arma::rowvec temp_zeta = q_vec.t()*zeta_res;
       // arma::rowvec mu_zetam = temp_zeta*C *Sigma_zetam* sigma_M2_inv;
       // zetam = arma::mvnrnd( mu_zetam.t(), Sigma_zetam);
       
     }
     
     
     // Rcout<<"1"<<std::endl;
     // 3. theta_eta
     
     if(iter%interval_eta ==0 & iter>begin_update_eta){
       theta_eta_res =  Mstar_alpha_term - theta_zetam*C_t;
       arma::colvec eta_Sigma_vec = 1/(sigma_M2_inv + sigma_eta2_inv/D);
       theta_eta_res.each_col() %= eta_Sigma_vec*sigma_M2_inv;
       if(update_eta_sgld){
         double step_eta = a_eta_step*pow((b_eta_step+iter),gamma_eta_step);
         
         arma::mat grad_theta_eta = theta_eta - theta_eta_res;
         grad_theta_eta.each_col() %= 1/eta_Sigma_vec;
         // double step_eta = step_allregion(r)/eta_step;
         arma::mat theta_eta_inc = -step_eta/2*grad_theta_eta + randn(size(grad_theta_eta))*sqrt(step_eta);
         theta_eta.rows(eta_Lidx_update) += theta_eta_inc.rows(eta_Lidx_update);
         
       }else{
         
         theta_eta = hyperplane_MVN_multiple(G,H_mat,eta_Sigma_vec,theta_eta_res.t());
         
       }
       
       // arma::mat theta_eta_new = arma::randn(L_max,n);
       // theta_eta_new.each_col() %= arma::sqrt(eta_Sigma_vec);
       // theta_eta_new += theta_eta_res;
       // theta_eta.rows(eta_cutoff_range) = theta_eta_new.rows(eta_cutoff_range);
       if(iter>= begin_save_eta){
         theta_eta_mean += theta_eta;
         eta_counter += 1;
       }
     }
     
     // Rcout<<"2"<<std::endl;
     // 4. sigma_alpha
     sigma_alpha2_inv = arma::randg( arma::distr_param(a + L_max*0.5, 1/(b + dot(theta_alpha,theta_alpha/D)/2)) );
     
     // 5. sigma_M
     M_reg_res = Mstar_alpha_term - theta_zetam*C_t - theta_eta;
     double M_norm = norm(M_reg_res,"fro");
     sigma_M2_inv = arma::randg( arma::distr_param(a + n*L_max/2,1/(b + M_norm*M_norm/2) ) );
     
     
     //   6. sigma_zetam -> half-cauchy
     // sigma_zetam2_inv = arma::randg( arma::distr_param(a + zetam.n_elem/2,1/(b + dot(zetam,zetam)/2) ) );
     // sigma_zetam2_inv_mcmc(iter) = sigma_zetam2_inv;
     // Rcout<<"3"<<std::endl;
     //   7. sigma_eta -> giving too large variance
     if(iter%interval_eta ==0 & iter>begin_update_eta){
       theta_eta_temp = theta_eta;
       theta_eta_temp.each_col() %= 1/sqrt(D);
       double eta_norm = norm(theta_eta_temp,"fro");
       sigma_eta2_inv = arma::randg( arma::distr_param(a + 0.5*n*L_max,
                                                       1/(b + eta_norm*eta_norm/2)) );
       
     }
     
     //
     // Rcout<<"4"<<std::endl;
     //     // --------------------- arma::summarize return --------------------- //
     logll_M = -0.5*square(arma::norm(M_reg_res,"fro"))*sigma_M2_inv;
     if(iter%thinning == 0){
       int iter_mcmc = iter/thinning;
       theta_alpha_mcmc.col(iter_mcmc) = theta_alpha;
       logll_mcmc(iter_mcmc) = logll_M + 0.5*n*L_max*log(sigma_M2_inv);
       theta_zetam_mcmc.slice(iter_mcmc) = theta_zetam;
       sigma_alpha2_inv_mcmc(iter_mcmc) = sigma_alpha2_inv;
       sigma_M2_inv_mcmc(iter_mcmc) = sigma_M2_inv;
       sigma_eta2_inv_mcmc(iter_mcmc) = sigma_eta2_inv;
       
       
     }
     
     int int_show = n_mcmc/20;
     if(iter%int_show==0){
       Rcout<<"iter="<<iter<<"; sigma_eta2="<<1/sigma_eta2_inv<<"; sigma_Y2="<<1/sigma_M2_inv<<std::endl;
     }
     
     
   }
   
   theta_eta_mean /= eta_counter;
   
   
   
   List gs = Rcpp::List::create(Rcpp::Named("theta_eta") = theta_eta,
                                Rcpp::Named("theta_zetam_mcmc")= theta_zetam_mcmc,
                                Rcpp::Named("sigma_M2_inv_mcmc")= sigma_M2_inv_mcmc ,
                                Rcpp::Named("sigma_alpha2_inv_mcmc")= sigma_alpha2_inv_mcmc,
                                Rcpp::Named("sigma_eta2_inv_mcmc")=  sigma_eta2_inv_mcmc);
   timer.step("end of iterations");
   return Rcpp::List::create(Rcpp::Named("theta_alpha_mcmc")=theta_alpha_mcmc,
                             Rcpp::Named("theta_eta_mean")=theta_eta_mean,
                             Rcpp::Named("logll_mcmc")=logll_mcmc,
                             Rcpp::Named("track_step")=track_step,
                             Rcpp::Named("accept_blcok")=accept_block,
                             Rcpp::Named("emp_accept")=emp_accept,
                             Rcpp::Named("gs")=gs,
                             Rcpp::Named("Timer")=timer
   );
   
 }

