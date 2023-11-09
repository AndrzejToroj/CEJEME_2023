############################### LIBRARIES #####################################
if(!require("mvtnorm")) {install.packages("mvtnorm"); library(mvtnorm)}
if(!require("dplyr")) {install.packages("dplyr"); library(dplyr)}
if(!require("tidyr")) {install.packages("tidyr"); library(tidyr)}
if(!require("readr")) {install.packages("readr"); library(readr)}
if(!require("lubridate")) {install.packages("lubridate"); library(lubridate)}
if(!require("readxl")) {install.packages("readxl"); library(readxl)}
if(!require("rgdal")) {install.packages("rgdal"); library(rgdal)}
if(!require("ggplot2")) {install.packages("ggplot2"); library(ggplot2)}
if(!require("hrbrthemes")) {install.packages("hrbrthemes"); library(hrbrthemes)}
if(!require("gridExtra")) {install.packages("gridExtra"); library(gridExtra)}
if(!require("stringr")) {install.packages("stringr"); library(stringr)}
if(!require("zoo")) {install.packages("zoo"); library(zoo)}
if(!require("RColorBrewer")) {install.packages("RColorBrewer"); library(RColorBrewer)}
if(!require("classInt")) {install.packages("classInt"); library(classInt)}
if(!require("geosphere")) {install.packages("geosphere"); library(geosphere)}
if(!require("rstudioapi")) {install.packages("rstudioapi"); library(rstudioapi)}
if(!require("truncnorm")) {install.packages("truncnorm"); library(truncnorm)}
if(!require("plyr")) {install.packages("plyr"); library(plyr)}
if(!require("nimble")) {install.packages("nimble"); library(nimble)}
if(!require("coda")) {install.packages("coda"); library(coda)}

############################### PARAMETERS ######################################
# the path with main folder, where all functions and result folders are located
#main_path <- "~/Desktop/Master_v2"
#main_path <- getActiveDocumentContext()$path

############################### SETTINGS ######################################

### REGIMES
# 1 - EXPANSION
# 0 - RECESSION
set.seed(42)
setwd(main_path)

############################### FUNCTIONS #####################################

# the function calculates..., the equation number 
make_residuals <- function(Y, 
                           S, 
                           theta, 
                           W) {
  N <- length(theta$omega_d)
  TT <- nrow(Y)
  ones <- matrix(1, nrow = TT, ncol = N)
  M_inv <- diag(N) - theta$rho * W
  epsilon <- (Y %*% t(M_inv) 
              - (ones-S) * kronecker(matrix(1,nrow=TT,ncol=1), t(as.matrix(theta$mu_0)))
              - S * kronecker(matrix(1,nrow=TT,ncol=1), t(as.matrix(theta$mu_1))) )
  return(epsilon)
}

rho_posterior_cond <- function(Y, 
                               S, 
                               theta, 
                               W) {
  N <- length(theta$omega_d)
  TT <- nrow(Y)
  det_M_inv <- det(diag(N) - theta$rho * W)
  epsilon <- make_residuals(Y, S, theta, W)
  post_rho <- TT * log(det_M_inv)
  inv_Omega <- solve(diag(theta$omega_d))
  for (tt in 1:TT) {
    post_rho <- post_rho - 0.5 * as.numeric(t(as.matrix(epsilon[tt, ])) %*% inv_Omega %*% as.matrix(epsilon[tt, ]))
  }
  return(post_rho)
}

loglik_conditional_S <- function(S, 
                                 Y, 
                                 theta, 
                                 W) {
  N <- length(theta$omega_d)
  TT <- nrow(Y)
  Omega <- diag(theta$omega_d)
  Omega_inv <- solve(Omega)
  M_inv <- diag(N) - theta$rho * W
  LL_lead_component <-  (- N/2*log(2*pi) 
                         - 0.5*sum(log(diag(Omega)))  #- 0.5*log(det(Omega))   SIMPLIFIED FOR DIAGONAL MATRIX (for num stability)
                         + log(det(M_inv)))
  epsilon <- make_residuals(Y, S, theta, W)
  
  loglik_cS <- 0
  for(tt in 1:TT) {
    loglik_cS <- loglik_cS + LL_lead_component - 0.5*t(as.matrix(epsilon[tt,])) %*% Omega_inv %*% as.matrix(epsilon[tt,])
  }
  return(loglik_cS)
}

Hamilton_filter <- function(Y, 
                            theta,
                            W) {
  
  N <- length(theta$omega_d)
  TT <- nrow(Y)
  
  prediction_0 <- matrix(NA, nrow = TT, ncol = N)
  prediction_1 <- matrix(NA, nrow = TT, ncol = N)
  update_0 <- matrix(NA, nrow = TT+1, ncol = N)
  update_1 <- matrix(NA, nrow = TT+1, ncol = N)
  density_0 <- matrix(NA, nrow = TT, ncol = N)
  density_1 <- matrix(NA, nrow = TT, ncol = N)
  
  P <- list()
  for (nn in 1:N) {
    P[[nn]] <- matrix(c(theta$p_00[nn], 1-theta$p_00[nn], 1-theta$p_11[nn], theta$p_11[nn]),
                      nrow = 2, byrow = TRUE)
    d_steady <- rep(NA, 2)
    d_steady[1] <- (1-P[[nn]][2,2]) / ((1-P[[nn]][2,2]) + (1-P[[nn]][1,1]))
    d_steady[2] <- (1-P[[nn]][1,1]) / ((1-P[[nn]][2,2]) + (1-P[[nn]][1,1]))
    update_0[1,nn] <- d_steady[1]
    update_1[1,nn] <- d_steady[2]
  }
  for (tt in 1:TT) {
    for (nn in 1:N) {
      prediction_0[tt,nn] <- update_0[tt,nn] * P[[nn]][1,1] + update_1[tt,nn] * P[[nn]][2,1]
      prediction_1[tt,nn] <- update_1[tt,nn] * P[[nn]][2,2] + update_0[tt,nn] * P[[nn]][1,2]
      density_0[tt,nn] <- dnorm(Y[tt,nn],
                                mean = theta$rho * t(as.matrix(W[nn,])) %*% as.matrix(Y[tt,]) + theta$mu_0[nn],
                                sd = sqrt(theta$omega_d[nn]))
      density_1[tt,nn] <- dnorm(Y[tt,nn],
                                mean = theta$rho * t(as.matrix(W[nn,])) %*% as.matrix(Y[tt,]) + theta$mu_1[nn],
                                sd = sqrt(theta$omega_d[nn]))
      temp_0 <- density_0[tt,nn] * prediction_0[tt,nn]
      temp_1 <- density_1[tt,nn] * prediction_1[tt,nn]
      update_0[tt+1,nn] <- temp_0 / (temp_0 + temp_1)
      update_1[tt+1,nn] <- temp_1 / (temp_0 + temp_1)
    }
  }
  update_0 <- update_0[-1,]
  update_1 <- update_1[-1,]
  
  update <- list(p_0 = update_0,
                 p_1 = update_1)
  return(update)
  
}

multiMoveGibbs <- function(theta, 
                           p_Hamilton, 
                           R) {
  N <- length(theta$omega_d)
  TT <- nrow(p_Hamilton$p_0)
  S <- list()
  for (rr in 1:R) {
    RNG_for_S <- matrix(runif(TT * N), nrow = TT, ncol = N)
    S_rr <- matrix(NA, nrow = TT, ncol = N)
    S_rr[TT,] <- ifelse(RNG_for_S[TT,] < p_Hamilton$p_0[TT,], 0, 1)
    for (tt in 1:(TT-1)) {
      p_trans_from0_tt <- ifelse(S_rr[TT-tt+1,]==0,theta$p_00,1-theta$p_00)
      p_trans_from1_tt <- ifelse(S_rr[TT-tt+1,]==0,1-theta$p_11,theta$p_11)
      p0_tt_temp <- p_trans_from0_tt * p_Hamilton$p_0[TT-tt,]
      p1_tt_temp <- p_trans_from1_tt * p_Hamilton$p_1[TT-tt,]
      p0_tt <- p0_tt_temp / (p0_tt_temp + p1_tt_temp)
      S_rr[TT-tt,] <- ifelse(RNG_for_S[TT-tt,] < p0_tt, 0, 1)
    }
    S[[rr]] <- S_rr
  }
  return(S)
}

loglik_uc <- function(theta, 
                      Y, 
                      W, 
                      R) {
  
  p_Hamilton <- Hamilton_filter(Y, theta, W)
  S_simul <- multiMoveGibbs(theta, p_Hamilton, R)
  loglik_cS <- lapply(S_simul, FUN = loglik_conditional_S, Y = Y, theta = theta, W = W)
  return(mean(unlist(loglik_cS)))
  
}

relist <- function(theta_vector, 
                   N) {
  theta <- list(rho = as.vector(theta_vector[1]),
                mu_1 = as.vector(theta_vector[2:(N+1)]),
                mu_0 = as.vector(theta_vector[(N+2):(2*N+1)]),
                omega_d = as.vector(theta_vector[(2*N+2):(3*N+1)]),
                p_00 = as.vector(theta_vector[(3*N+2):(4*N+1)]),
                p_11 = as.vector(theta_vector[(4*N+2):(5*N+1)]))
  return(theta)
}

sample_posterior <- function(initial, 
                             hyperpar, 
                             S, 
                             S0, 
                             S_rho, 
                             S0_rho, 
                             Y, 
                             W) {
  
  #Step 0: imply from arguments and initialize result matrix
  N <- ncol(Y)
  TT <- nrow(Y)
  simulation <- data.frame(matrix(nrow = S, ncol = length(unlist(initial))))
  colnames(simulation) <- names(unlist(initial))
  rownames(simulation) <- 1:S
  simulation[1,] <- unlist(initial)
  
  lowerbound_rho <- 1/min(eigen(W)$values)
  
  #Step 1: hyperparameters for priors
  attach(hyperpar, warn.conflicts = FALSE)
  
  for (ss in 2:S) {

    #Step 2: draw S from MM-Gibbs
    theta_prev <- relist(as.numeric(simulation[ss-1,]), N)
    p_Hamilton <- Hamilton_filter(Y, theta_prev, W)
    S_simul <- multiMoveGibbs(theta_prev, p_Hamilton, 1)[[1]]

    #Steps 3-4: sample p_00 and p_11 
    transitions <- ifelse(S_simul[2:TT,] == 0 & S_simul[1:(TT-1),] == 0, 0,
                          ifelse(S_simul[2:TT,] == 0 & S_simul[1:(TT-1),] == 1, 10,
                                 ifelse(S_simul[2:TT,] == 1 & S_simul[1:(TT-1),] == 0, 1,
                                        ifelse(S_simul[2:TT,] == 1 & S_simul[1:(TT-1),] == 1, 11, NA))))
  
    trans_sum <- apply(data.frame(transitions), 2, FUN = table, simplify = FALSE)#, vars = colnames(transitions))
    rm(transitions)
    n <- as.data.frame(c(0,1,10,11))
    colnames(n) <- c("x")
    
    for (ii in 1:length(trans_sum)) {
      trans_sum[[ii]] <- as.data.frame(trans_sum[[ii]])
      colnames(trans_sum[[ii]]) <- c("x", paste0("freq",ii))
      n <- merge(x = n, y = trans_sum[[ii]], by.x = "x", by.y = "x", all = TRUE)
    }
    n[is.na(n)] <- 0
    rm(trans_sum)
    rownames(n) <- paste0("n", c(0,1,10,11))
    n <- n[,-1]
    
    alpha_posterior <- n + kronecker(matrix(1, nrow = 1, ncol = N), rbind(as.matrix(alpha_prior[1,]), as.matrix(alpha_prior[2,])))
    for (nn in 1:N) {
      p_11_nn <- rbeta(1, alpha_posterior[4,nn], alpha_posterior[3,nn])
      p_00_nn <- rbeta(1, alpha_posterior[1,nn], alpha_posterior[2,nn])
      simulation[ss, 3*N+1+nn] <- p_00_nn
      simulation[ss, 4*N+1+nn] <- p_11_nn
    }
    
    #Step 5: sample sigma^2
    v_posterior <- v_prior + TT
    epsilon <- make_residuals(Y, S_simul, theta_prev, W)
    delta_posterior <- delta_prior + diag(var(epsilon)) * (TT-1)
    for (nn in 1:N) {
      sigma2_nn <- rinvgamma(n=1, shape = v_posterior/2, scale = delta_posterior[nn]/2)
      simulation[ss, 2*N+1+nn] <- sigma2_nn
    }
    
    #Step 6: sample mu
    X_n <- list()
    m_posterior_n <- list()
    M_posterior_n <- list()
    for (nn in 1:N) {
      X_n[[nn]] <- cbind(matrix(1,nrow = TT) - S_simul[,nn], S_simul[,nn])
      sigma2_nn <- simulation[ss, 2*N+1+nn]
      rho <- simulation[ss-1, 1]
      M_posterior_n[[nn]] <- solve(solve(M_prior) + (sigma2_nn^(-1))*t(X_n[[nn]]) %*% X_n[[nn]])
      m_posterior_n[[nn]] <- M_posterior_n[[nn]] %*% (solve(M_prior) %*% m_prior + (sigma2_nn^(-1)) * t(X_n[[nn]]) %*% (as.matrix(Y[,nn] - rho * Y %*% as.matrix(W[nn,]))))
      mu_nn <- rmvnorm(n=1, mean=m_posterior_n[[nn]], sigma = M_posterior_n[[nn]])
      simulation[ss, 1+nn] <- mu_nn[2]
      simulation[ss, N+1+nn] <- mu_nn[1]
    }
    
    #Step 7: rho (M-H)
    simul_rho <- rep(NA, S_rho)
    acceptance_p <- rep(NA, S_rho-1)
    
    theta_now <- relist(as.numeric(simulation[ss,]), N)
    theta_now[[1]] <- as.numeric(simulation[ss-1,1])
    theta_cand <- theta_now
    simul_rho[1] <- as.numeric(theta_prev[1])
    MH_accepted <- FALSE
    
    while(MH_accepted == FALSE) {
      
      if(is.na(acceptance_p[1])) {
        c <- 0.00015 #0.00026 for GDP, 0.00015 for UE
        print(paste0("Initiating c to ", c))
      } else {
        if(mean(acceptance_p) < 0.2) {
          c <- c * 0.66
          print(paste0("Downscaling c to ", c))
        } 
        if(mean(acceptance_p) > 0.4) {
          c <- c * 1.5
          print(paste0("Upscaling c to ", c))
        }
        if(mean(acceptance_p) >= 0.2 && mean(acceptance_p) <= 0.4) {
          MH_accepted <- TRUE
          print(paste0("Terminating procedure with c=", c, " and avg accept p = ", mean(acceptance_p)))
        }
      }
      
      if(MH_accepted == FALSE) {
        for (sr in 2:S_rho) {
          rho_candidate <- rtruncnorm(1, a=lowerbound_rho, b=1, mean = simul_rho[sr-1], sd = c*100)
          theta_cand[[1]] <- rho_candidate
          log_post_now <- rho_posterior_cond(Y, S_simul, theta_now, W)
          log_post_cand <- rho_posterior_cond(Y, S_simul, theta_cand, W)
          acceptance_p[sr-1] <- min(exp(log_post_cand - log_post_now), 1)
          if (runif(1) <= acceptance_p[sr-1]) {
            simul_rho[sr] <- rho_candidate
            theta_now <- theta_cand
          } else {
            simul_rho[sr] <- simul_rho[sr-1]
          }
        }
        print(paste0("Avg acceptance prob: ", mean(acceptance_p)))
      }
    }
    simul_rho <- simul_rho[(S0_rho+1):length(simul_rho)]
    draw_rho <- simul_rho[ceiling(runif(1, min = 0, max = length(simul_rho)))]
    
    # if (ss==S-500){
    #   png(file = paste0("trace_rho_",ss,".png"), width = 300, height = 300)
    #   plot( simul_rho, type="s", xpd=NA, ylab="Parameter", xlab="Sample", las=1)
    #   hist(simul_rho, 50, freq=FALSE, main="", las=1,
    #        xlab="x", ylab="Probability density")
    #   
    #   mh.draws <- mcmc(simul_rho)
    #   png(file = paste0("geweke_plot_",ss,".png"), width = 300, height = 300)
    #   geweke.plot(mh.draws, 
    #               frac1 = 0.1, 
    #               frac2 = 0.5,
    #               nbins=40, 
    #               pvalue=0.05)
    #   dev.off()
    #   
    #   png(file = paste0("autocorr_plot_",ss,".png"), width = 300, height = 300)
    #   autocorr.plot(mh.draws)
    #   dev.off()
      
    #  conv_test<-summary(mh.draws)
    #}
    
    simulation[ss, 1] <- draw_rho
    
    print(paste0("### Iteration ", ss))
    
  }
  output<-list()
  output$simul <- simulation[(S0+1):nrow(simulation),]
  #output$convergence <- conv_test
  #output$chain <- mh.draws
  return(output)
}

#----------------------------- W MATRIX - distance -----------------------------#
matrix_W_distance <- function(map,
                              g = 1,
                              f = distCosine){
  
  distance <- distm(coordinates(map), 
                    fun = f)
  rownames(distance) <- map@data$ID
  colnames(distance) <- map@data$ID
  gamma <- g
  W <- 1 / (distance ^ gamma)
  diag(W) <- 0
  W<- W / as.matrix(rowSums(W)) %*% matrix(1, nrow=1, ncol=nrow(distance))
  return(W)
}

draw_impulse2 <- function(map, N, n, theta, W, ef, r, legend, i,
         n_color=neg_color, p_color=pos_color,
         m_point=middle_point, last_point=max_point){
  
  # impulse <- abs(theta$mu_1- theta$mu_0)
  # #pal[1:m_point] <- brewer.pal(m_point, n_color)
  # #pal <- rev(pal)
  # #pal[(m_point+1):last_point] <- brewer.pal(last_point-m_point, p_color)
  pal <- brewer.pal(last_point, n_color)
  # impulse2 <- as.matrix(rep(0,N))
  # impulse2[i] <- impulse[i]
  map@data$response <- as.vector(ef[,i])
  map@data$bracket <- cut(map@data$response, r)
  spplot(map, "bracket", lwd=0.1, col.regions=pal, colorkey=legend,
         par.settings=list(axis.line = list(col =  'transparent')),
         main=list(label=n[i], cex=0.8, fontfamily="serif"))
}

