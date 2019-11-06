library(rstan)
library(parallel)

rstan_options(auto_write = TRUE)
options(mc.cores = min(2, detectCores()))

mixed_model_data = list(
  Y = as.vector(t(Cmatrix)), #Cmatrix has i -> j  in the j,i spot, so transpose
  D = 2,
  N = nrow(Y),
  sigma_a = 2,
  sigma_b = 2,
  sigma_z = 2,
  sigma_beta0 = 3,
  intercept = TRUE,
  zero_constraint = TRUE,
  self_edges = FALSE,
  dist = 1 #1 is euclidean, 2 is squared euclidean (bc stan only allows numeric data)
  
)

tmp = stan(file = "../latent_net_stan/mixed_model_stan.stan", model_name = "mixed_test", data = mixed_model_data)
