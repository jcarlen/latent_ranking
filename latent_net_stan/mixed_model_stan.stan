// Jane Carlen, November 2019
// Code to implement latent mixed modesls in stan

// TO DO:
// 0) Invariance of positions wrt rotation?
//    Change contstraint to zero-sum? (simplex plus scale parameter?)
//
// 1) distances 
//      i)  euclidean distance
//      ii) squared eucliden distant
//      
// 2) edge distributions:
//      i) binomial/bernoulli
//      ii) negative binomial 
//
// 3) covariates
//      i) nodal covariates
//      ii) edge covariates
//---------------------------------------------------------------------------------------------------------------


// lambda = BX + 
// intercept or not?
// constrain sender and receiver by a_1, b_1 = 0
// directed only
// also constrain z_1? 

// The input data is a vector 'y' of length 'N'.
data {
  int<lower=2> N; //network size
  int<lower=0> D; //latent symmetric dimensions
  //int<lower=0> R; //latent assymmetric dimensions
  //int<lower=0> p; //number of covariates, including intercept
  int Y[N*N]; //observed edges, as vector
  real<lower=0> sigma_beta0; //prior sd of intercept
  int<lower=0,upper=1>intercept; //prior sd of position
  int<lower=0,upper=1>zero_constraint; // let sender_1 = 0, receiver_1 = 0
  int<lower=0,upper=1>pos_constraint; // let 1st - Dth positions be on primary axes?
  int<lower=0,upper=1>self_edges; 
  int<lower=1,upper=2>dist; // 1 is euclidean, 2 is squared euclidean
}

transformed data {
  int N_0 = N - zero_constraint;
  int N_Z = N - pos_constraint*D;
}

// The parameters accepted by the model. Our model
// accepts two parameters 'mu' and 'sigma'.
parameters {
  real beta0; //intercept,if present
  real sender[N_0];
  real receiver[N_0];
  vector[D] Z[N_Z]; //array of size N_Z, contaning vectors with D elements
  real<lower=0> Z_constrained[D*pos_constraint];
  real mu_sender;
  real mu_receiver;
  real mu_Z[D];
  real<lower=0> sigma_sender; //prior sd of sender
  real<lower=0> sigma_receiver; //prior sd of receiver 
  real<lower=0> sigma_z; //prior sd of position
  
}

transformed parameters {
  
  real lambda[N*N] = rep_array(0.0, N*N);
  
  // Calculate lambda
  for (i in 1:N) { 
    for (j in 1:N) { 
      
      if (intercept) {
        lambda[(i-1)*N + j] += beta0; 
      }
      
      if (!zero_constraint) {
         lambda[(i-1)*N + j] +=  sender[i] + receiver[j];
      } else {
          if (i != 1) { lambda[(i-1)*N + j] +=  sender[i - zero_constraint];}
          if (j != 1) { lambda[(i-1)*N + j] +=  receiver[j - zero_constraint];}
      }


      if (dist == 1) { //euclidean
        if (!pos_constraint) {
          lambda[(i-1)*N + j] +=  -distance(Z[i], Z[j]);
        } else {
          
          if (i <= D && j > D) {
            vector[D] fixed = rep_vector(0.0, D);
            fixed[i] = Z_constrained[i];
            lambda[(i-1)*N + j] +=  -distance(fixed, Z[j - D]);
          }
          
          if (j <= D && i > D) {
            vector[D] fixed = rep_vector(0.0, D);
            fixed[j] = Z_constrained[j];
            lambda[(i-1)*N + j] +=  -distance(fixed, Z[i - D]);
          }
            
          if (i <= D && j <= D) {
            lambda[(i-1)*N + j] += -distance([0,Z_constrained[j]], [Z_constrained[i],0]);  
          }
        }
      }

        
      if (dist == 2) { //squared
        if (!pos_constraint) {
          lambda[(i-1)*N + j] +=  -squared_distance(Z[i], Z[j]);
        } else {
          
          if (i <= D && j > D) {
            vector[D] fixed = rep_vector(0.0, D);
            fixed[i] = Z_constrained[i];
            lambda[(i-1)*N + j] +=  -squared_distance(fixed, Z[j - D]);
          }
          
          if (j <= D && i > D) {
            vector[D] fixed = rep_vector(0.0, D);
            fixed[j] = Z_constrained[j];
            lambda[(i-1)*N + j] +=  -squared_distance(fixed, Z[i - D]);
          }
            
          if (i <= D && j <= D) {
            lambda[(i-1)*N + j] += -squared_distance([0,Z_constrained[j]], [Z_constrained[i],0]);  
          }
        }
      }
      
      lambda[(i-1)*N + j] = exp(lambda[(i-1)*N + j]);
      
    }
  }
  
}


model {
  
  sigma_sender ~ scaled_inv_chi_square_lpdf(3,1); //prior sd of sender
  sigma_receiver ~ scaled_inv_chi_square_lpdf(3,1); //prior sd of receiver 
  sigma_z ~ scaled_inv_chi_square_lpdf(sqrt(N),N/8); //prior sd of receiver 
  
  mu_sender ~ normal(0, 3);
  mu_receiver ~ normal(0, 3);
  mu_Z ~ normal(0,2);
  
  if (intercept) { beta0 ~ normal_lpdf(0, sigma_beta0); }
  sender ~ normal_lpdf(mu_sender, sigma_sender);
  receiver ~ normal_lpdf(mu_receiver, sigma_receiver);
  
  for (i in 1:N) {
    if (pos_constraint && i <= D) {
      Z_constrained[i] ~ normal(0, sigma_z);
    } else {
      Z[i - pos_constraint*D] ~ normal_lpdf(mu_Z, sigma_z);
    }
  }
  
  Y ~ poisson_lpmf(lambda);
  
  // penalty to prevent reflection
  // target += -100*Z[10][1];
  
  if (!self_edges) {
    for (i in 1:N) {
      target += -poisson_lpmf(Y[(i-1)*N+i] | lambda[(i-1)*N+i]);
    }
  }

}

