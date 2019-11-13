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

// directed only
// also constrain z_1? 

data {
  int<lower=2> N; //network size
  int<lower=0> D; //latent symmetric dimensions
  //int<lower=0> R; //latent assymmetric dimensions
  //int<lower=0> p; //number of covariates, including intercept
  int Y[N*N]; //observed edges, as vector
  real<lower=0> sigma_beta0; //prior sd of intercept
  int<lower=0,upper=1>self_edges; 
  int<lower=1,upper=2>dist; // 1 is euclidean, 2 is squared euclidean
}

// The parameters accepted by the model. Our model
// accepts two parameters 'mu' and 'sigma'.
parameters {
  real beta0; //intercept,if present
  simplex[N] sender_raw;  // https://mc-stan.org/docs/2_20/stan-users-guide/parameterizing-centered-vectors.html
  simplex[N] receiver_raw;
  real<lower=0> sender_scale;
  real<lower=0> receiver_scale;
  //matrix[N,D] Z;
  simplex[N] Z_raw[D]; //D simplex vectors of length N
  positive_ordered[D] Z_scale;
}

transformed parameters {
  
  real lambda[N*N] = rep_array(0.0, N*N);
  vector[N] sender = (sender_raw-rep_vector(inv(N), N))*sender_scale;
  vector[N] receiver = (receiver_raw-rep_vector(inv(N), N))*receiver_scale;
  
  matrix[N,D] Z;
  for (i in 1:D) {
    Z[,i] = Z_raw[i]*Z_scale[i];
  }
      
  // Calculate lambda
  for (i in 1:N) { 
    for (j in 1:N) { 
      
      lambda[(i-1)*N + j] += beta0 + sender[i] + receiver[j];

      if (dist == 1) { //euclidean
          lambda[(i-1)*N + j] +=  -sqrt(squared_distance(Z[i,], Z[j,]) + .0001); //issues with zero distance (gradient at endpoint?)
      }
       
      if (dist == 2) { //squared
         lambda[(i-1)*N + j] +=  -squared_distance(Z[i,], Z[j,]);
       }
      
      lambda[(i-1)*N + j] = exp(lambda[(i-1)*N + j]);
    }
  }
  
}


model {
  
  beta0 ~ normal(0, sigma_beta0);
  sender_raw ~ dirichlet(rep_vector(1, N)); //prior of sender
  receiver_raw ~ dirichlet(rep_vector(1, N)); //prior of sender
  sender_scale ~ normal(0, 5);
  receiver_scale ~ normal(0, 5);
  
  for (i in 1:D) {
      //Z[,i] ~ normal(0,2);
      Z_raw[i] ~ dirichlet(rep_vector(1, N)); 
  }
  
  //Z_scale ~ normal(0,5);
  
  Y ~ poisson(lambda);
  
  if (!self_edges) {
    for (i in 1:N) {
      target += -poisson_lpmf(Y[(i-1)*N+i] | lambda[(i-1)*N+i]);
    }
  }

}


