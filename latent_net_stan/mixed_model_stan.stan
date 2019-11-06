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
  real<lower=0> sigma_a; //prior sd of sender
  real<lower=0> sigma_b; //prior sd of receiver 
  real<lower=0> sigma_z; //prior sd of position
  real<lower=0> sigma_beta0; //prior sd of intercept
  int<lower=0,upper=1>intercept; //prior sd of position
  int<lower=0,upper=1>zero_constraint;
  int<lower=0,upper=1>self_edges; 
  int<lower=1,upper=2>dist; // 1 is euclidean, 2 is squared euclidean
}

transformed data {
  int N2 = N - zero_constraint;
}

// The parameters accepted by the model. Our model
// accepts two parameters 'mu' and 'sigma'.
parameters {
  real beta0; //intercept,if present
  real sender[N2];
  real receiver[N2];
  vector[D] Z[N2];
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
        if (!zero_constraint) {
          lambda[(i-1)*N + j] +=  -distance(Z[i], Z[j]);
        } else {
          if (i == 1 && j !=1) { lambda[(i-1)*N + j] +=  -distance(rep_vector(0.0, D), Z[j - zero_constraint]);}
          if (i != 1 && j ==1) { lambda[(i-1)*N + j] +=  -distance(Z[i - zero_constraint], rep_vector(0.0, D));}
          // if i and j both 0 distance is zero
        }
      }
      
      if (dist == 2) { //squared eucliden
        if (!zero_constraint) {
          lambda[(i-1)*N + j] +=  -squared_distance(Z[i], Z[j]);
        } else {
          if (i == 1 && j !=1) { lambda[(i-1)*N + j] +=  -squared_distance(rep_vector(0.0, D), Z[j - zero_constraint]);}
          if (i != 1 && j ==1) { lambda[(i-1)*N + j] +=  -squared_distance(Z[i - zero_constraint], rep_vector(0.0, D));}
          // if i and j both 0 distance is zero
        }
      }
      
      lambda[(i-1)*N + j] = exp(lambda[(i-1)*N + j]);
      //add in covariates
    }
  }
  
}


model {
  if (intercept) { beta0 ~ normal_lpdf(0, sigma_beta0); }
  sender ~ normal_lpdf(mu_a, sigma_a);
  receiver ~ normal_lpdf(mu_b, sigma_b);
  for (i in 1:N2) {
    Z[i] ~ normal_lpdf(0, sigma_z);
  }
  
  Y ~ poisson_lpmf(lambda);
  
  if (!self_edges) {
    for (i in 1:N) {
      target += -poisson_lpmf(Y[(i-1)*N+i] | lambda[(i-1)*N+i]);
    }
  }

}

