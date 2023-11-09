############################### SETTINGS ######################################
# identifiers of dataset
country <- "PL" # PL for Poland, USA for the United States (not in CEJEME article)
variable <- "GDP" # UE for unemployment rate, GDP for Gross Domestic Product
# the path with main folder, where all functions and result folders are located
main_path <- "C:/Andrzej/OneDrive - SGH/moje_papery/MG_MP/replication_pack_CEJEME"

S <- 1350#20000 # a number of iterations in simulation
S0 <- 0#4000 # burn-in for Gibbs
S_rho <- 10000 # a number of iterations for rho simulation
S0_rho <- 3000 # burn-in for rho simulation
chains <- TRUE # if TRUE, two chains will be simulated (rather than one)
continue <- TRUE
if (continue == TRUE) {
  S0 <- 0
}

############################### LIBRARIES #####################################
library(mvtnorm)
library(plyr)
library(dplyr)
library(tidyr)
library(readr)
library(lubridate)
library(readxl)
library(rgdal)
library(stringr)
library(zoo)
library(classInt)
library(geosphere)

############################### PARAMETERS #####################################
print("Setting parameters for simulation")
# the path with saved data (matrix W and datasets)
basic_file <- paste0(main_path, "/basic_data.RData")
# the path to save the result of first chain
save_first_chain <- paste0(main_path, "/post_simul/posteriorA_", country, "_", variable, ".RData")
save_first_chain0 <- paste0(main_path, "/post_simul/posterior0A_", country, "_", variable, ".RData")
# the path to save the result of second chain
save_second_chain <- paste0(main_path, "/post_simul/posteriorB_", country, "_", variable, ".RData")
save_second_chain0 <- paste0(main_path, "/post_simul/posterior0B_", country, "_", variable, ".RData")

############################### SETTINGS ######################################
print("Loading data and functions")
set.seed(42)
setwd(main_path)
load(basic_file)
source("0_MC_MG_HF_conv_functions.R")
dataframe <- paste(country, variable, "ch", sep = "_")
N <- n_regions # n_regions for Poland, n_states for USA
eval(parse(text = paste0("W <- W_", country)))

IDs <- unique(PL_GDP_ch$ID)
set_prior <- matrix(NA, nrow = 5, ncol = N)
ii <- 0
for (id in IDs) {
  ii <- ii+1
  set_prior[1,ii] <- as.numeric(quantile(PL_GDP_ch$Value[PL_GDP_ch$ID == id], probs = 0.2))
  set_prior[2,ii] <- as.numeric(quantile(PL_GDP_ch$Value[PL_GDP_ch$ID == id], probs = 0.8))
  set_prior[3,ii] <- as.numeric(mean(PL_GDP_ch$Value[PL_GDP_ch$ID == id]))
  set_prior[4,ii] <- (set_prior[1,ii]+set_prior[3,ii])/2 
  set_prior[5,ii] <- (set_prior[2,ii]+set_prior[3,ii])/2 
}
set_prior_avg <- rowMeans(set_prior)


settings_list <- list()
settings_list$PL_UE <- list(rho = 0.5, # spatial autoregression parameter
                            mu_1 = rep(-1.3, N), # constant during expansion
                            mu_0 = rep(0.2, N), # constant during recession
                            omega_d = rep(1, N), # variances (already squared)
                            p_00 = rep(0.8, N), # probability of staying in expansion
                            p_11 = rep(0.8, N),
                            alpha_prior = matrix(c(8, 2, 1, 9), nrow = 2, byrow = TRUE),
                            v_prior = 6,
                            delta_prior = 10,
                            m_prior = matrix(c(0.2, -1.3), nrow = 2),
                            M_prior = diag(2))
settings_list$PL_GDP <- list(rho = 0.5,
                             mu_1 = rep(9.6, N),
                             mu_0 = rep(2.9, N),
                             omega_d = rep(1, N), 
                             p_00 = rep(0.8, N),
                             p_11 = rep(0.8, N),
                             alpha_prior = matrix(c(8, 2, 2, 8), nrow = 2, byrow = TRUE),
                             v_prior = 6,
                             delta_prior = 50,
                             #m_prior = matrix(c(1,10), nrow = 2),
                             m_prior = matrix(c(2,9), nrow = 2),
                             #M_prior = diag(2)
                             M_prior = 2*diag(2))
settings_list$USA_GDP <- list(rho = 0.5,
                              mu_1 = rep(4, N),
                              mu_0 = rep(2.3, N),
                              omega_d = rep(1, N),
                              p_00 = rep(0.8, N),
                              p_11 = rep(0.8, N),
                              alpha_prior = matrix(c(8, 2, 1, 9), nrow = 2, byrow = TRUE),
                              v_prior = 6,
                              delta_prior = 2,
                              m_prior = matrix(c(2.3,4), nrow = 2),
                              M_prior = diag(2))
settings_list$USA_UE <- list(rho = 0.5,
                             mu_1 = rep(-0.2, N),
                             mu_0 = rep(0.3, N),
                             omega_d = rep(1, N),
                             p_00 = rep(0.8, N),
                             p_11 = rep(0.8, N),
                             alpha_prior = matrix(c(8, 2, 1, 9), nrow = 2, byrow = TRUE),
                             v_prior = 6,
                             delta_prior = 8,
                             m_prior = matrix(c(0.5,-1), nrow = 2),
                             M_prior = diag(2))



if (dataframe=="PL_UE_ch"){
  Y <- PL_UE_ch
  settings <- settings_list$PL_UE
} else if (dataframe=="PL_GDP_ch"){
  Y <- PL_GDP_ch
  settings <- settings_list$PL_GDP
} else if (dataframe=="USA_GDP_ch"){
  Y <- USA_GDP_ch
  settings <- settings_list$USA_GDP
} else {
  Y <- USA_UE_ch
  Y$Period <- Y$Date
  settings <- settings_list$USA_UE
}

######################### START VALUES & HYPERPARAMETERS #######################
### Numbers used for regimes:
# 1 - EXPANSION
# 0 - RECESSION
# List of start values for parameter
theta0 <- list(rho = settings$rho,
               mu_1 = settings$mu_1,
               mu_0 = settings$mu_0,
               omega_d = settings$omega_d,
               p_00 = settings$p_00,
               p_11 = settings$p_11)
print("Start values for parameters were loaded")
if(continue == TRUE) {
  simulation_path <- paste0(main_path, "/post_simul")
  eval(parse(text = paste0("simulation_file <- paste0(simulation_path, '/posterior0A_", country, "_", variable, ".RData')")))
  load(simulation_file) 
  N_a <- nrow(chainA$simul)
  theta0 <- list(rho = chainA$simul$rho[N_a],
                 mu_1 = chainA$simul[N_a,2:(N+1)],
                 mu_0 = chainA$simul[N_a,(N+2):(2*N+1)],
                 omega_d = chainA$simul[N_a,(2*N+2):(3*N+1)],
                 p_00 = chainA$simul[N_a,(3*N+2):(4*N+1)],
                 p_11 = chainA$simul[N_a,(4*N+2):(5*N+1)])
  print("Start values for parameters overwritten with chain A endpoint")
}

# List of hyperparameters (priors)
hyperpar0 = list(alpha_prior = settings$alpha_prior,
                 v_prior = settings$v_prior,
                 delta_prior = settings$delta_prior,
                 m_prior = settings$m_prior,
                 M_prior = settings$M_prior)
print("List of hyperparameters (priors) were loaded")


########################### RUNNING SIMULATIONS #############################

#------------------------------- MATRIX Y ----------------------------------#
Y <- Y %>% 
  select(c(Name, Period, Value)) %>% 
  pivot_wider(names_from = Name, values_from = Value)
Y <- as.matrix(Y[,-1])

#------------------------------- FIRST CHAIN -------------------------------#
print(paste0("Start of first simulation for ", variable, " in ", country))
start <- Sys.time()
chainA0 <- sample_posterior(initial = theta0, 
                                hyperpar = hyperpar0, 
                                S = S, 
                                S0 = S0, 
                                S_rho = S_rho, 
                                S0_rho = S0_rho, 
                                Y = Y, 
                                W = W)
end <- Sys.time()
print("Simulation 1 took ", end - start)
if(continue==TRUE) {
  colnames(chainA0$simul) <- colnames(chainA$simul)
  chainA$simul <- rbind(chainA$simul, chainA0$simul[-1,])
} else {
  chainA <- chainA0
}
rm(chainA0)
save.image(save_first_chain)
save(chainA, file = save_first_chain0)

#------------------------------- SECOND CHAIN --------------------------------#
if (chains==TRUE){
  load(save_first_chain0)
  theta_0 <- colMeans(chainA$simul)
  theta0 <- relist(theta_0, N)
  if(continue == TRUE) {
    load(save_second_chain0)
    N_b <- nrow(chainB$simul)
    theta0 <- list(rho = chainB$simul$rho[N_b],
                   mu_1 = chainB$simul[N_b,2:(N+1)],
                   mu_0 = chainB$simul[N_b,(N+2):(2*N+1)],
                   omega_d = chainB$simul[N_b,(2*N+2):(3*N+1)],
                   p_00 = chainB$simul[N_b,(3*N+2):(4*N+1)],
                   p_11 = chainB$simul[N_b,(4*N+2):(5*N+1)])
    print("Start values for parameters overwritten with chain B endpoint")
  }
  
  print(paste0("Start of second simulation for ", variable, " in ", country))
  start <- Sys.time()
  chainB0 <- sample_posterior(initial = theta0, 
                                  hyperpar = hyperpar0, 
                                  S = S, 
                                  S0 = S0, 
                                  S_rho = S_rho, 
                                  S0_rho = S0_rho, 
                                  Y = Y, 
                                  W = W)
  end <- Sys.time()
  print("Simulation 2 took ", end - start)
  if(continue==TRUE) {
    colnames(chainB0$simul) <- colnames(chainB$simul)
    chainB$simul <- rbind(chainB$simul, chainB0$simul[-1,])
    
  } else {
    chainB <- chainB0
  }
  rm(chainB0)
  save.image(save_second_chain)
  save(chainB, file = save_second_chain0)
  
}
