library(ggplot2)
library(viridis)
library(BayesGPfit)
library(RSpectra)


rowcol_idx = function(idx_i, mat){
  n_row = nrow(mat)
  n_col = ncol(mat)
  col_idx = floor(idx_i/n_row)
  row_idx = idx_i - col_idx*n_row
  col_idx = col_idx + 1
  
  cbind(row_idx , col_idx)
}


plot_joint = function(joint_bayes,datsim){
  par(mfrow = c(3,3))
  plot(joint_bayes$mcmc$loglik_Y[-1],main = "loglik_Y")
  plot(joint_bayes$mcmc$theta_beta[1,-1],main = "theta_beta[1,-1]")
  plot(1/joint_bayes$mcmc$inv_sigmasq_Y,main ="sigmasq_Y")
  
  plot(joint_bayes$mcmc$loglik_M[-1],main = "loglik_M[-1]")
  plot(joint_bayes$mcmc$theta_alpha[1,-1],main="theta_alpha[1,-1]")
  plot(1/joint_bayes$mcmc$inv_sigmasq_M[-1],main = "sigmasq_M[-1]")
  
  if(is.na(any(joint_bayes$mcmc$theta_eta_mean))){
    # could plot the initial value
  }else{
    plot(c(joint_bayes$mcmc$theta_eta_mean), c(datsim$theta_eta),main = "theta_eta v truth");abline(0,1,col="red")
  }
  
  plot(joint_bayes$mcmc$sigma_eta[-1],main = "sigma_eta[-1]")
  plot(apply(joint_bayes$mcmc$theta_nu,1,mean),datsim$true_params$theta_nu,main = "theta_nu v truth");abline(0,1,col="red")
  par(mfrow = c(1,1))
}


matern_kernel = function(x,y,nu,l=1){
  d = sqrt(sum((x-y)^2))/l
  y = 2^(1-nu)/gamma(nu)*(sqrt(2*nu)*d)^nu*besselK(sqrt(2*nu)*d,nu)
  return(y)
}

#' generate_matern_basis2
#' 
#' @export
generate_matern_basis2 = function(grids, region_idx_list, L_vec,scale = 2,nu = 1/5,
                                  show_progress = FALSE){
  if(nu=="vec"){
    nu_vec = region_idx_list["nu_vec"]
  }
  num_block = length(region_idx_list)
  Phi_D = vector("list",num_block)
  Phi_Q = vector("list",num_block)
  Lt = NULL; pt = NULL
  for(i in 1:num_block){
    if(show_progress){
      print(paste("Computing basis for block ",i))
    }
    p_i = length(region_idx_list[[i]])
    kernel_mat = matrix(NA,nrow = p_i, ncol=p_i)
    for(l in 1:p_i){
      if(nu=="vec"){
        kernel_mat[l,] = apply(grids[region_idx_list[[i]],],1,matern_kernel,y=grids[region_idx_list[[i]],][l,],nu = nu_vec[i],l=scale)
      }else{
        kernel_mat[l,] = apply(grids[region_idx_list[[i]],],1,matern_kernel,y=grids[region_idx_list[[i]],][l,],nu = nu,l=scale)
      }
    }
    diag(kernel_mat) = 1
    K = eigs_sym(kernel_mat,L_vec[i])
    K_QR = qr(K$vectors)
    Phi_Q[[i]] = qr.Q(K_QR )
    Phi_D[[i]] = K$values
    Lt = c(Lt, length(Phi_D[[i]]))
    pt = c(pt, dim(Phi_Q[[i]])[1])
  }
  return(list(Phi_D = Phi_D,
              region_idx_block = region_idx_list,
              Phi_Q = Phi_Q,L_all = Lt,p_length=pt))
}

#' @title plot_img
#' @param img A vector of input image
#' @param grids_df A data frame to indicate the position of pixels
#' @import ggplot2
#' @import viridis
#' @export
plot_img = function(img, grids_df,title="img",col_bar = NULL, legend = T){
  if(length(img)!=dim(grids_df)[1]){
    print("dimension of img and grids do not match!")
  }else{
    g = ggplot(grids_df, aes(x=x1,y=x2)) +
      geom_tile(aes(fill = img)) +
      scale_fill_viridis_c(limits = col_bar, oob = scales::squish)+
      ggtitle(title)+
      theme(plot.title = element_text(size=20),legend.text=element_text(size=10))
    if(!legend){
      g = g + theme(legend.position="none")
    }
  }
  
  g
  
}
Soft_threshold = function(x,lambda){
  return( (x-sign(x)*lambda)*(abs(x)>lambda))
}

#' @title generate image mediation data with unobserved confounders
#' @export
Unconfounded_mediation_data_constrained_eta = function(beta, alpha,n,Q,D,
                                                       lambda,nu = NULL,
                                                       lambda_nu=0,
                                                       q=2,
                                           sigma_M = 0.1, 
                                           sigma_Y = 0.1, 
                                           nu_scale = 1,
                                           sigma_eta = 0.1,
                                           gamma = 2){
  
  p = length(beta)
  X = rnorm(n)
  C = matrix(rnorm(n*q),q,n)*0.1 # q by n
  L = length(D)
  
  # get constrained eta
  G = rbind(X,C)
  H_mat = BASMU::get_H_mat(G); # from Unconfounded_BIMA.cpp
  theta_eta = hyperplane_MVN_multiple(G, H_mat, sigma_eta^2*D, matrix(0,n,L) ); # output L by n
  eta = Q %*% theta_eta
  
  zeta_y = sample(-10:10,q)
  zeta_m = matrix(rnorm(p*q),p,q)
  # nu = rep(0,p)
  # nu[sample(1:p, round(p*0.2))] = 0.5
  # nu[beta > quantile(beta,0.9)] = 1
  # ---test 1
  if(is.null(nu)){
    nu_tilde = simulate_round_image(side = sqrt(p), range = c(0,0.5),
                                    center_shift = c(0.3,0.3));nu_tilde = nu_tilde$img
    nu = nu_tilde
  }
  # plot_img(nu_tilde,grids_df)
  # ---test 2
  # nu[sample(1:p, round(p*0.2))] = 0.5
  # nu[beta <= quantile(beta,0.1)] = 1
  # ---end of test
  theta_nu = t(Q)%*%nu
  nu = Q %*% theta_nu
  nu = nu*nu_scale
  
  
  alpha_STGP = Q%*%t(Q)%*%alpha; alpha_STGP = Soft_threshold(alpha_STGP[,1], lambda)
  beta_STGP = Q%*%t(Q)%*%beta; beta_STGP = Soft_threshold(beta_STGP[,1], lambda)
  nu_STGP = Q%*%t(Q)%*%nu; nu_STGP = Soft_threshold(nu_STGP[,1], lambda_nu)
  
  # generate image
  M = alpha_STGP %*% t(X) + zeta_m %*% C + eta + matrix(rnorm(p*n, sd = sigma_M),p,n)
  
  # generate outcome
  Y = t(beta_STGP) %*% M + t(c(gamma, zeta_y)) %*% rbind(t(X),C) + t(nu_STGP) %*% eta  + rnorm(n, sd = sigma_Y)
  res = Y - (t(beta_STGP) %*% M + t(c(gamma, zeta_y)) %*% rbind(t(X),C) + t(nu_STGP) %*% eta)
  betaM = t(beta_STGP) %*% M
  gammaX = t(c(gamma, zeta_y)) %*% rbind(t(X),C)
  nu_eta = t(nu_STGP) %*% eta
  
  # true_params
  true_params = list(beta = beta_STGP, alpha = alpha_STGP, nu = nu_STGP,
                     zeta_y = zeta_y, zeta_m = zeta_m,
                     theta_nu = theta_nu,lambda_nu=lambda_nu,
                     gamma = gamma, sigma_M = sigma_M, sigma_Y = sigma_Y)
  
  return(list(Y = Y, M = t(M), X = X, C = C, eta = eta, theta_eta = theta_eta, nu = nu, res=res,
              betaM = betaM, gammaX = gammaX, nu_eta = nu_eta,
              true_params = true_params))
}


Unconfounded_mediation_data_eta = function(beta, alpha,n,q=2,sigma_M = 0.1, 
                                           Q,D,
                                           Q_eta,D_eta,
                                       sigma_Y = 0.1, rho = 0.5, gamma = 2){
  
  p = length(beta)
  X = rnorm(n)
  C = matrix(rnorm(n*q),q,n)*0.1 # q by n
  L_eta = dim(Q_eta)[2]
  eta = Q_eta%*%matrix(rnorm(L_eta*n)*rep(sqrt(D_eta),n),nrow = L_eta, ncol = n)
  
  zeta_y = sample(-10:10,q)
  zeta_m = matrix(rnorm(p*q),p,q)
  nu = rep(0,p)
  nu[sample(1:p, round(p*0.2))] = 0.5
  nu[beta > quantile(beta,0.9)] = 1
  
  alpha_STGP = Q%*%t(Q)%*%alpha
  beta_STGP = Q%*%t(Q)%*%beta
  nu_STGP = Q%*%t(Q)%*%nu
  
  # generate image
  M = alpha_STGP %*% t(X) + zeta_m %*% C + eta + matrix(rnorm(p*n, sd = sigma_M),p,n)
  
  # generate outcome
  Y = t(beta_STGP) %*% M + t(c(gamma, zeta_y)) %*% rbind(t(X),C) + t(nu_STGP) %*% eta  + rnorm(n, sd = sigma_Y)
  
  # true_params
  true_params = list(beta = beta_STGP[,1], alpha = alpha_STGP[,1], zeta_y = zeta_y, zeta_m = zeta_m,
                     gamma = gamma, sigma_M = sigma_M, sigma_Y = sigma_Y,nu = nu_STGP[,1])
  
  return(list(Y = Y, M = t(M), X = X, C = C, eta = eta, nu = nu, 
              true_params = true_params))
}

scalar_Unconfounded_mediation_data = function(beta_img, alpha_img,n,q=2,
                                              sigma_M = 0.1, sigma_Y = 0.1, 
                                              gamma = 2, nu = -2){
  
  p = length(beta_img)
  
  # generate exposure
  X = rnorm(n)
  C = matrix(rnorm(n*q),q,n)*0.1 # q by n
  U = rnorm(n)
  
  zeta_y = sample(-10:10,q)
  zeta_m = matrix(rnorm(p*q),p,q)
  
  
  
  # generate image
  M = alpha_img %*% t(X) + zeta_m %*% C + t(replicate(p,U)) + 
      matrix(rnorm(p*n, sd = sigma_M),p,n)
  
  # generate outcome
  Y = t(beta_img) %*% M + t(c(gamma, zeta_y)) %*% rbind(t(X),C) + 
      nu * U  + rnorm(n, sd = sigma_Y)
  
  # true_params
  true_params = list(beta = beta_img, alpha = alpha_img, 
    zeta_y = zeta_y, zeta_m = zeta_m, nu = nu,
                     gamma = gamma, sigma_M = sigma_M, sigma_Y = sigma_Y)
  
  return(list(Y = Y, M = t(M), X = X, C = C, U = U,  true_params = true_params))
}

convert_to_long = function(datsim, M_res = NULL){

  n = length(datsim$Y)
  if(is.null(M_res)){
    wide_df = as.data.frame(cbind(id = 1:n,Y = datsim$Y[1,], M = datsim$M,
                                X = datsim$X, C = t(datsim$C)))
  }else{
    wide_df = as.data.frame(cbind(id = 1:n,Y = datsim$Y[1,], M_res = M_res,
                                X = datsim$X, C = t(datsim$C)))
  }
  data_long <- tidyr::gather(wide_df, "loc", "intensity", V3:V402, factor_key=TRUE)
  return(list(datsim_long = data_long, 
    datsim_wide = wide_df))

}

Unconfounded_mediation_data = function(beta_img, alpha_img,n,q=2,q_U = 3,sigma_M = 0.1,
                                       sigma_Y = 0.1, rho = 0.5, gamma = 2){
  
  p = length(beta_img)
  # X = rnorm(n)
  # C = matrix(rnorm(n*q),q,n)*0.1 # q by n
  # U = matrix(runif(n*q_U),q_U,n)*0.1 # q by n
  
  # generate exposure
  Sigma = matrix(rep(rho,(q+q_U+1)^2),q+q_U+1,q+q_U+1)
  diag(Sigma) =  rep(1,q+q_U+1)
  Exposure = mvtnorm::rmvnorm(n, sigma = Sigma)
  X = Exposure[,1]
  C = t(Exposure[,1:q+1])
  U = t(Exposure[,1:q_U+1+q])
  
  zeta_y = sample(-10:10,q)
  zeta_m = matrix(rnorm(p*q),p,q)
  zeta_mU = matrix(rnorm(p*q_U),p,q_U)*5
  zeta_yU = sample(-10:10,q_U)*5
  
  
  # generate image
  M = alpha_img %*% t(X) + zeta_m %*% C + zeta_mU %*% U + matrix(rnorm(p*n, sd = sigma_M),p,n)
  
  # generate outcome
  Y = t(beta_img) %*% M + t(c(gamma, zeta_y)) %*% rbind(t(X),C) + t(zeta_yU) %*% U  + rnorm(n, sd = sigma_Y)
  
  # true_params
  true_params = list(beta = beta_img, alpha = alpha_img, zeta_y = zeta_y, zeta_m = zeta_m,
                     gamma = gamma, sigma_M = sigma_M, sigma_Y = sigma_Y)
  
  return(list(Y = Y, M = t(M), X = X, C = C, U = U, zeta_mU = zeta_mU,zeta_yU = zeta_yU, true_params = true_params))
}

#' simulate_round_image
#' 
#' @import BayesGPfit
#' @export
simulate_round_image = function(center_shift = c(0,0),lambda = 0.1,side = 30, 
                                range = c(0,1)){
  n_sqrt = side
  n = n_sqrt*n_sqrt
  grids = GP.generate.grids(d=2L,num_grids=n_sqrt)
  center = apply(grids,2,mean) + center_shift
  rad = apply(grids,1,function(x){sum((x-center)^2)})
  inv_rad = 2-rad
  inv_rad_ST = Soft_threshold(inv_rad,1.2)
  f_mu = Soft_threshold(log(inv_rad_ST^2+1),lambda)
  
  y = f_mu
  nonzero = y[abs(y)>0]
  a = range[1]; b = range[2]
  nonzero_mapped = (nonzero-min(nonzero))/(max(nonzero)-min(nonzero))*(b-a) + a
  y[abs(y)>0] = nonzero_mapped
  
  grids_df = as.data.frame(grids)
  
  
  return(list(img = y, grids_df = grids_df))
}

#' simulate_triang_image
#' 
#' @import BayesGPfit
#' @export
simulate_triang_image = function(center_shift = c(0,0),side = 30, 
                                 rad = 0.5, radiance_bool = T,
                                 range = c(0,1)){
  n_sqrt = side
  n = n_sqrt*n_sqrt
  grids = GP.generate.grids(d=2L,num_grids=n_sqrt)
  center = apply(grids,2,mean) + center_shift
  
  
  bool1 = 1*(grids[,2] > -rad/2)
  bool2 = 1*(grids[,2] < -sqrt(3)*grids[,1]+rad)
  bool3 = 1*(grids[,2] < sqrt(3)*grids[,1]+rad)
  
  beta = bool1*bool2*bool3
  
  if(radiance_bool){
    radiance = apply(grids,1,function(x){sum((x-center)^2)})
    inv_rad = 2-radiance
    inv_rad_ST = Soft_threshold(inv_rad,1.2)
    beta = beta*inv_rad_ST
  }
  
  
  grids_df = as.data.frame(grids)
  
  
  return(list(img = beta, grids_df = grids_df))
}


# previous help functions -------------------------------------------------

#' STGP_mcmc
#' 
#' @export
STGP_mcmc = function(theta_mcmc_sample,region_idx,basis,lambda){
  M = dim(theta_mcmc_sample)[2]
  S = max(unlist(region_idx))
  num_region = length(region_idx)
  est_mcmc = matrix(NA,nrow = S, ncol = M)
  dd = matrix(unlist(lapply(basis$Phi_Q,dim)),nrow=2)
  L_idx = cumsum(dd[2,])
  L_idx = c(rbind(L_idx,L_idx+1))
  L_idx = matrix(c(1,L_idx[-length(L_idx)]),nrow=2)
  for(l in 1:num_region){
    idx = region_idx[[l]]
    theta = theta_mcmc_sample[L_idx[1,l]:L_idx[2,l],]
    beta = basis$Phi_Q[[l]]%*%theta
    est_mcmc[idx,] = (beta-sign(beta)*lambda)*I(abs(beta)>lambda)
  }
  return(est_mcmc)
}

InclusionMap = function(mcmc_sample, true_beta, thresh = "auto", fdr_target = 0.1,
                        max.iter = 100){
  InclusionProb = 1 - apply(mcmc_sample, 1, function(x){mean(abs(x)==0)})
  true_beta = 1*(true_beta!=0)
  thresh_final = thresh
  fdr=NA
  if(thresh=="auto"){
    thresh = 0.5
    for(i in 1:max.iter){
      mapping = 1*(InclusionProb>thresh)
      fdr = FDR(mapping, true_beta)
      print(paste("fdr=",fdr,"thresh=",thresh))
      if(is.na(fdr)){
        print("fdr=NA, target FDR is too small")
        thresh = thresh/1.1
        mapping = 1*(InclusionProb>thresh)
        fdr = FDR(mapping, true_beta)
        print(paste("Use current fdr=",fdr,"thresh=",thresh))
        break
      }
      if(fdr<=fdr_target){
        thresh_final = thresh
        break
      }
      thresh = thresh*1.1
      if(thresh>1){
        print("New thresh>1, keep thresh at the current value and return result.")
        break
      }
    }
  }else{
    mapping = 1*(InclusionProb>thresh)
  }
  return(list(mapping = mapping, thresh=thresh_final,
              InclusionProb=InclusionProb))
}
FDR = function(active_region, true_region){
  sum(active_region!=0 & true_region==0)/sum(active_region!=0)
}
Precision = function(active_region, true_region){
  mean(I(active_region!=0) == I(true_region!=0))
}
Power = function(active_region, true_region){
  sum(active_region !=0 & true_region!=0)/sum(true_region!=0)
}

#' plot_multi_img
#' 
#' @export
plot_multi_img = function(list_of_image, grids_df, n_img_per_row = 3,
                          col_bar = NULL, layout_matrix=NULL, 
                          labels = NULL,
                          font_size=20,
                          color_mode = "BWR",
                          legend_position = "right"){
  n_img = length(list_of_image)
  if(is.null(col_bar)){
    col_bar = range(unlist(list_of_image))
  }
  n_row_img = ceiling(n_img/n_img_per_row)
  all_image_p = vector("list", n_img)
  if(is.null(labels)){
    labels = names(list_of_image)
  }
  for(i in 1:n_img){
    all_image_p[[i]] <-  local({
      i <- i
      
        # ggplot2::scale_fill_viridis_c(limits = col_bar, oob = scales::squish,name = "Value")
      if(color_mode == "BWR"){
        ggplot(grids_df, aes(x=x1,y=x2)) +
          geom_tile(aes(fill = list_of_image[[i]] )) +
          ggtitle(labels[i])+
          theme(plot.title = element_text(size=font_size),
                legend.text = element_text(size=font_size*0.7))+
          ggplot2::scale_fill_gradient2(low = "blue",
                                        mid="white",
                                        high="red",
                                        midpoint = 0, limits = col_bar,
                                        name = "Value")
      }else{
        ggplot(grids_df, aes(x=x1,y=x2)) +
          geom_tile(aes(fill = list_of_image[[i]] )) +
          ggtitle(labels[i])+
          theme(plot.title = element_text(size=font_size),
                legend.text = element_text(size=font_size*0.7))+
          ggplot2::scale_fill_viridis_c(limits = col_bar, oob = scales::squish,name = "Value")
      }
        
    })
    
  }
  
  ggpubr::ggarrange(plotlist = all_image_p, 
                    ncol=n_img_per_row, 
                    nrow=n_row_img, 
                    common.legend = TRUE, legend=legend_position)
  
  
  
}

plot_multi_img2 = function(list_of_image, grids_df, n_img_per_row = 3,
                          col_bar = c(0,1), byrow = T,
                          labels = NULL,
                          add_legend=F){
  n_img = length(list_of_image)
  n_row_img = ceiling(n_img/n_img_per_row)
  all_image = vector("list", n_img)
  if(is.null(labels)){
    names(all_image) = names(list_of_image)
  }else{
    names(all_image) =  labels
  }
  for(i in 1:n_img){
    img = list_of_image[[i]]
    title = names(list_of_image)[i]
    g = ggplot(grids_df, aes(x=x1,y=x2)) +
      geom_tile(aes(fill = img)) +
      scale_fill_viridis_c(limits = col_bar, oob = scales::squish)+
      ggtitle(title)+
      theme(plot.title = element_text(size=20),legend.text=element_text(size=10))
    if(add_legend){
      g = g+theme(legend.position="right")
    }
      # theme(legend.position="none")
    all_image[[i]] = ggplotGrob(g)
  }
  layout_matrix = matrix(c(1:n_img,rep(NA,n_row_img*n_img_per_row-n_img)), 
                         nrow = n_row_img, ncol = n_img_per_row, byrow = byrow)
  
  # layout_matrix[1:n_img] = 1:n_img
  gridExtra::grid.arrange(
    grobs = all_image,
    layout_matrix = layout_matrix
  )
  
}
