# Helper likelihood function for latent space rating
#
#    Diag of network matrix representation should be 0 if not counting self-edges (!)
#
#    if est is "MAP" return l(Y|theta) + l(theta) ; if "Y" return l(Y|theta) ; if "theta" return l(theta)
#    if 'object' has parameter values they will be used (not including hyperparameters)
#    family = poisson fits the normal poisson latent space model
#    family = binomial fits a quasi-Stiglery (quasi-symmetric) type model with positions -- off by a constant

library(extraDistr) # for divnchisq
library(gtools) #for logit function
#----------------------------------------------------------------------------------------------------
# 1. llik, llik2 ----------------------------------------------------------------------------------------------------
#   assumes no self loops
#   (for weighting?) add a constant to the data matrix (Y) if it has zeros (should be >=1 so no neg weights)
#       this was found to be advantageous for residuals anyway
#       BUT it seems to compromise the separation of clusters
#   faster to calculate dist (uses c) than call it pair by pair or do pair by pair in R. 
#   long term can write rcpp function for just the distances of interest  
#   do rows or columns in the data have noticeably higher variance?

llik <- function(object=NULL, Y=NULL, sender=NULL, receiver=NULL, beta=NULL,
                 Z=NULL, sender.var = 10, receiver.var = 10, Z.var = 10,
                 beta.var = 9, sender.var.df = 3, receiver.var.df = 3, Z.var.df = NULL,
                 prior.sender.var = 1, prior.receiver.var = 1, prior.Z.var = NULL,
                 est = "MAP", family = "poisson") {
  
  if(is.null(Y)) {Y = object$model$Ym}
  if(is.null(Z)) {Z = object$Z}
  if(is.null(sender)) {sender = object$sender}
  if(is.null(receiver)) {receiver = object$receiver}
  if(is.null(beta)) {beta = object$beta}
  N = nrow(Y) 
  if(is.null(Z.var.df)) {Z.var.df = sqrt(N)}
  if(is.null(prior.Z.var)) { prior.Z.var = N/8}
  
  if(!is.null(object$sender.var)) sender.var = object$sender.var
  if(!is.null(object$receiver.var)) receiver.var = object$receiver.var
  if(!is.null(object$Z.var)) Z.var = object$Z.var
  if(!is.null(object$beta.var)) beta.var = object$beta.var
  
  Z_dist = as.matrix(dist(Z, upper = T))
  
  if (length(sender) == N) {
    #l_lambda = t(receiver + t(sender - Z_dist)) + beta;
    l_lambda = beta + outer(sender, receiver, "+") - Z_dist #slightly faster
    
  } else {
    sr = lapply(1: (length(sender)/N), function(x) {
      outer(sender[ (N*(x-1)+1) : (N*(x))], receiver[(N*(x-1)+1) : (N*(x))], 
          function(s,r) {(s+r)^2})
      })
    sr = sqrt(Reduce ("+", sr))
    l_lambda = beta + sr - Z_dist
  }
  
  if (family == "poisson") {
    lgamma.constant = sum(lgamma(as.vector(Y+1)), na.rm = T) 
    lambda = exp(l_lambda); diag(lambda) = 0
    pY = sum( Y * l_lambda - lambda, na.rm = T) - lgamma.constant
  }

  if (family == "binomial") {
    lambda = inv.logit(l_lambda)
    Yt = Y + t(Y); diag(Yt) =  0
    pY =  sum( Y * log(lambda), na.rm = T) + sum((Yt - Y)*log(1-lambda), na.rm = T)
  }   
  
  if (est == "Y") {return(pY)}
  
  ptheta = log(exp(-beta^2/(2*beta.var)) / sqrt(2*pi*beta.var)) +
    sum(log(exp(-sender^2/(2*sender.var)) / sqrt(2*pi*sender.var))) + 
    sum(log(exp(-receiver^2/(2*receiver.var)) / sqrt(2*pi*receiver.var))) +
    sum(log(exp(-Z^2/(2*Z.var)) / sqrt(2*pi*Z.var))) +
    log(dinvchisq(sender.var, sender.var.df, prior.sender.var)) + 
    log(dinvchisq(receiver.var, receiver.var.df, prior.receiver.var)) + 
    log(dinvchisq(Z.var, Z.var.df, prior.Z.var))
  
  if (est == "theta") {return(ptheta)}
  
  map = pY + ptheta # = p(Y|theta) + p(theta)
  
  if (est == "MAP") {return(map)}
  
}

#for using llik with optim:
llik2 <- function(theta, Y, d, est = "MAP", family = "poisson") {
  
  n = nrow(Y)
  nparam = length(theta)
  
  B = theta[1]
  Z = matrix(theta[2:(1+d*n)], ncol = d, nrow = n)
  #split whatever's left
  a = theta[(2+d*n): (1+d*n/2 + (nparam-1)/2) ]
  b = theta[(2+d*n/2 + (nparam-1)/2) : nparam]
  
  d2 = length(a)/n
  
  if (est == "MAP") {
    sigma2_a = theta[nparam-2]
    sigma2_b = theta[nparam-1]
    sigma2_z = theta[nparam]
  }
  
  #exact update conditioned on graph parameters
  if (est == "MAPe") {
    sigma2_a = (sum(a^2) + sender.var.df*prior.sender.var^2) / (n*d2 + 2 + prior.sender.var)
    sigma2_b = (sum(b^2) + receiver.var.df*prior.receiver.var^2) / (n*d2 + 2 + prior.receiver.var)
    sigma2_z = (sum(Z^2) + Z.var.df*prior.Z.var^2) / (n*d + 2 + prior.Z.var)
  } 

  if (est == "MAPe") {est = "MAP"} # to evaluate likelihood correctly
  
  return(llik(Y=Y, sender = a, receiver = b, beta = B, Z = Z,
              sender.var = sigma2_a, receiver.var = sigma2_b,
              Z.var = sigma2_z, family = family, est = est))
}

#----------------------------------------------------------------------------------------------------
# 2. llik_hat, llik2_hat  -----------------------------------------------------------------------------------------------
#    add a constant to the data matrix (Y) if it has zeros (should be >=1 so no neg weights)
#       this was found to be advantageous for residuals anyway
#    faster to calculate dist (uses c) than call it pair by pair or do pair by pair in R. 
#    long term can write rcpp function for just the distances of interest  
#
#    user rows or columans? do rows or columns in the data have noticeably higher variance?
#----------------------------------------------------------------------------------------------------
llik_hat <- function(object=NULL, Y=NULL, sender=NULL, receiver=NULL, beta=NULL,
                 Z=NULL, sender.var = 10, receiver.var = 10, Z.var = 10,
                 beta.var = 9, sender.var.df = 3, receiver.var.df = 3, Z.var.df = NULL, #N = number of nodes
                 prior.sender.var = 1, prior.receiver.var = 1, prior.Z.var = NULL,
                 est = "MAP", family = "poisson",
                 r = .25, rtype = "constant", Rmin = 5, replace = T, margin = "none", W = NULL) {
  
  if(is.null(Y)) {Y = object$model$Ym}
  if(is.null(Z)) {Z = object$Z}
  if(is.null(sender)) {sender = object$sender}
  if(is.null(receiver)) {receiver = object$receiver}
  if(is.null(beta)) {beta = object$beta}
  N = nrow(Y) 
  if(is.null(Z.var.df)) {Z.var.df = sqrt(N)}
  if(is.null(prior.Z.var)) { prior.Z.var = N/8}
  
  if(!is.null(object$sender.var)) sender.var = object$sender.var
  if(!is.null(object$receiver.var)) receiver.var = object$receiver.var
  if(!is.null(object$Z.var)) Z.var = object$Z.var
  if(!is.null(object$beta.var)) beta.var = object$beta.var
  
  # don't calculate whole z-dist matrix?:
  # actually faster to (since it does it all in C with one call)
  Z_dist = as.matrix(dist(Z, upper = T))
  
  # poisson by default
  
  # Sampling  -----
  # Moved this outside so it doesn't have to be recalculated every time
  # Y = Y + c # add a constant to Y if it has zeros (should be >=1 so no neg weights)
  # diag(Y) = 0 # assume no self-loops
  # lgamma.constant = sum(lgamma(as.vector(Y+1))) # constant for factorial term
  # W = Y*(log(Y)-1)   
  # diag(W) = 0
  
  # Number to sample
  if (rtype == "constant") {R = rep(round((N-1)*r), N)} else {
    R = Rmin + round(apply(Y, 1, var)/sum(apply(Y, 1, var)) * N * (r-Rmin))
  }
  
  if (replace) { #hansen hurwitz only works with sampling by replacement
    if (is.null(W)) { # do hansen hurwitz, weighted by Y -------
      
      if (margin == "row" | margin == "none") { # by row 
        #assign(s, t(sapply(1:N, function(x){sample(N, size = R[x], replace = T, prob = Y[x,]) })))
        pY_hat  = sapply(1:N, function(x){
              s = sample(N, size = R[x], replace = T, prob = Yw[x,]) 
              l_lambda = beta + sender[x]  + receiver[s] - Z_dist[x,s]
              #Yw_row[x] * mean( (l_lambda - exp(l_lambda)/Y[x,s]) ) #y_i' = y_i/M_i = 
              Yw_row[x] * mean( (Y[x,s]*l_lambda - exp(l_lambda))/Yw[x,s]) 
              # unbiased estimator of variance
              # (Yw_row[x]^2/(R[x]*(R[x]-1))) *sum( ( (Y[x,s]*l_lambda - exp(l_lambda))/Yw[x,s] - mean((Y[x,s]*l_lambda - exp(l_lambda))/Yw[x,s]))^2)
              })
        pY_hat = sum(pY_hat)  - lgamma.constant
      }
      
      if (margin == "col" | margin == "none") { # by col
        pY_hat_col  = sapply(1:N, function(x){
             s = sample(N, size = R, replace = T, prob = Yw[,x]) 
             l_lambda = beta + sender[s]  + receiver[x] - Z_dist[x,s]
            
             if( Yw_col[x] * mean( (Y[s,x]*l_lambda - exp(l_lambda))/Yw[s,x] )  < -1.153803e+30 ) {print("-inf")}
             
             Yw_col[x] * mean( (Y[s,x]*l_lambda - exp(l_lambda))/Yw[s,x] ) 
             # unbiased estimator of variance
             # (Yw_col[x]^2/(R[x]*(R[x]-1))) *sum( ( (Y[s,x]*l_lambda - exp(l_lambda))/Yw[s,x] - mean((Y[s,x]*l_lambda - exp(l_lambda))/Yw[s,x]))^2)
             
          })
        pY_hat_col = sum(pY_hat_col)  - lgamma.constant
      }
  
      if (margin == "col") {pY_hat = pY_hat_col}
      if (margin == "none") {pY_hat = mean(pY_hat, pY_hat_col)}
      
    } else {
      # hansen hurwitz, weighted by W  -------
      
      if (margin == "row" | margin == "none") {
        pY_hat  = sapply(1:N, function(x){
          s = sample(N, size = R, replace = T, prob = W[x,]) 
          l_lambda = beta + sender[x]  + receiver[s] - Z_dist[x,s]
          W_row[x] * mean( (Y[x,s]*l_lambda - exp(l_lambda))/W[x,s] ) #y_i' = y_i/M_i = 
          #( 1/(R*(R-1)) ) * sum((( Y[x,s]*l_lambda - exp(l_lambda))/W[x,s] - 
          #     W_row[x]*mean( (Y[x,s]*l_lambda - exp(l_lambda) )/W[x,s]) )^2) #variance
        })
        pY_hat = sum(pY_hat)  - lgamma.constant
      }
      
      if (margin == "col" | margin == "none") {
      
      # by col
        pY_hat_col  = sapply(1:N, function(x){
          s = sample(N, size = R, replace = T, prob = W[,x]) 
          l_lambda = beta + sender[s]  + receiver[x] - Z_dist[x,s]
          W_col[x] * mean( (Y[s,x]*l_lambda - exp(l_lambda))/W[s,x] ) #y_i' = y_i/M_i = 
          #( 1/(R*(R-1)) ) * sum((( Y[s,x]*l_lambda - exp(l_lambda))/W[s,x] - 
          #     W_col[x]*mean( (Y[s,x]*l_lambda - exp(l_lambda) )/W[s,x]) )^2) #variance
        })
        pY_hat_col = sum(pY_hat_col)  - lgamma.constant
      }
      
      if (margin == "col") {pY_hat = pY_hat_col}
      if (margin == "none") {pY_hat = mean(pY_hat, pY_hat_col)}
      
    }
  } 
  
  if (!replace) {
    # ad hoc, without replacement, but slower? var estimate? ----
    
    # Y weights
    if (margin == "row" | margin == "none") {
      pY_hat  = sapply(1:N, function(x){
        Rx = R[x]
        s = sample(N, size = Rx, replace = F, prob = Y[x,]) 
        l_lambda = beta + sender[x]  + receiver[s] - Z_dist[x,s]
        Yw_row2 = Yw_row[x] - c(0,cumsum(Y[x,s[-Rx]]))
        mean(Yw_row2 * (l_lambda - exp(l_lambda)/Y[x,s]) + c(0, cumsum(l_lambda[-Rx]*Y[x,s[-Rx]] - exp(l_lambda[-Rx]) ) ) )
      })
      pY_hat = sum(pY_hat)  - lgamma.constant
    }
    
    if (margin == "col" | margin == "none") {
      pY_hat_col  = sapply(1:N, function(x){
        Rx = R[x]
        s = sample(N, size = Rx, replace = F, prob = Y[,x]) 
        l_lambda = beta + sender[s]  + receiver[x] - Z_dist[x,s]
        Yw_col2 = Yw_col[x] - c(0,cumsum(Y[s[-Rx], x]))
        mean(Yw_col2 * (l_lambda - exp(l_lambda)/Y[s,x]) + c(0, cumsum(l_lambda[-Rx]*Y[s[-Rx], x] - exp(l_lambda[-Rx]) ) ) )
      })
      pY_hat_col = sum(pY_hat_col)  - lgamma.constant
    }
    
    if (margin == "col") {pY_hat = pY_hat_col}
    if (margin == "none") {pY_hat = mean(pY_hat, pY_hat_col)}
  }
    
  # add for binomial?
  # return ----
  if (est == "Y") {return(pY_hat)}
  
  ptheta = log(exp(-beta^2/(2*beta.var)) / sqrt(2*pi*beta.var)) +
    sum(log(exp(-sender^2/(2*sender.var)) / sqrt(2*pi*sender.var))) + 
    sum(log(exp(-receiver^2/(2*receiver.var)) / sqrt(2*pi*receiver.var))) +
    sum(log(exp(-Z^2/(2*Z.var)) / sqrt(2*pi*Z.var))) +
    log(dinvchisq(sender.var, sender.var.df, prior.sender.var)) + 
    log(dinvchisq(receiver.var, receiver.var.df, prior.receiver.var)) + 
    log(dinvchisq(Z.var, Z.var.df, prior.Z.var))
  
  if (est == "theta") {return(ptheta)}
  
  map = pY_hat + ptheta # = p(Y|theta) + p(theta)
  
  if (est == "MAP") {return(map)}
  
}


#for using llik_hat with optim:
llik_hat2 <- function(theta, Y, d, est = "Y", family = "poisson",
                      r = .25, rtype = "constant", Rmin = 5, replace = T, margin = "none",  W = NULL ) {
  
  n=nrow(Y)
  B = theta[1]
  Z = matrix(theta[2:(1+d*n)], ncol = d, nrow = n)
  a = theta[(2+d*n):(n+1+d*n)]
  b = theta[(n+2+d*n):(2*n+1+d*n)]
  
  if (est == "MAP") {
    sigma2_a = theta[2*n+2+d*n]
    sigma2_b = theta[(2*n+3+d*n)]
    sigma2_z = theta[(2*n+4+d*n)]
  }
  
  if (est == "MAPe") {
    sigma2_a = (sum(a^2) + sender.var.df*prior.sender.var^2) / (n + 2 + prior.sender.var)
    sigma2_b = (sum(b^2) + receiver.var.df*prior.receiver.var^2) / (n + 2 + prior.receiver.var)
    sigma2_z = (sum(Z^2) + Z.var.df*prior.Z.var^2) / (n*d + 2 + prior.Z.var)
  } 
  
  if (est == "MAPe") {est = "MAP"} # to evaluate likelihood correctly
  
  return(llik_hat(Y=Y, sender = a, receiver = b, beta = B, Z = Z, 
                  sender.var = 1, receiver.var =1, Z.var = 1,
                  family = family, est = est,
                  r = r, rtype = rtype, replace = T, margin = margin, Rmin = Rmin, W = W)
         )
}


#----------------------------------------------------------------------------------------------------
# 3. gradient: llik_gr, llik_gr_hat -----------------------------------------------------------------------------------------------

llik_gr <- function(theta, Y, d, est = "MAP") {
  
  n = nrow(Y)
  nparam = length(theta)
  
  B = theta[1]
  Z = matrix(theta[2:(1+d*n)], ncol = d, nrow = n)
  #split whatever's left
  a = theta[(2+d*n): (1+d*n/2 + (nparam-1)/2) ]
  b = theta[(2+d*n/2 + (nparam-1)/2) : nparam]
  
  d2 = length(a)/n
  
  if (est == "MAP") {
    sigma2_a = theta[nparam-2]
    sigma2_b = theta[nparam-1]
    sigma2_z = theta[nparam]
  }
  
  #exact update of variance parameters conditioned on graph parameters
  if (est == "MAPe") {
    #optimal values
    sigma2_a = (sum(a^2) + sender.var.df*prior.sender.var^2) / (n*d2 + 2 + prior.sender.var)
    sigma2_b = (sum(b^2) + receiver.var.df*prior.receiver.var^2) / (n*d2 + 2 + prior.receiver.var)
    sigma2_z = (sum(Z^2) + Z.var.df*prior.Z.var^2) / (n*d + 2 + prior.Z.var)
  }
  

  if (length(a == N) ) {

    #lambda = exp(t(b + t(a - Z_dist)) + B); diag(lambda) = 0
    lambda = B + outer(a, b, "+") - Z_dist; diag(lambda) = 0 #slightly faster
    
    da = rowSums(Y) - rowSums(lambda)
    
    db = colSums(Y) - colSums(lambda)
    
  } else {
    
    sr2 = lapply(1: (length(a)/N), function(x) {
      outer(a[ (N*(x-1)+1) : (N*(x))], b[(N*(x-1)+1) : (N*(x))], "+")
    })
    
    sr = sqrt(Reduce ("+", lapply(sr2, "^", 2)))
    sr_inv = 1/sr; diag(sr_inv) = rep(0, N)
    
    lambda = exp(B + sr - Z_dist); diag(lambda) = 0
    
    da = unlist(lapply(sr2, function(x) {rowSums(x*sr_inv*(Y - lambda))}))
      
    db = 
      #l_lambda = t(receiver + t(sender - Z_dist)) + beta;
      l_lambda = B + outer(sender, receiver, "+") - Z_dist #slightly faster
      
    } else {
      
      
      l_lambda = 
  }
  
  #add sampling?
  Z_dist = dist2(Z)
  dist_inv = 1/Z_dist; diag(dist_inv) = rep(0, N) #term1 = cbind(term1, term1)
  tmp1 = (Y + t(Y)) * dist_inv
  tmp2 = dist_inv * (lambda + t(lambda))
  zid_zjd = lapply(1:N, function(x) {t(Z[x,] - t(Z))}) # N, N x d matrices
  tmp3 = as.vector(t( sapply(1:N, function(i) {colSums(tmp2[i,]*zid_zjd[[i]])}) - #d x N -> N x d
                        sapply(1:N, function(i) {colSums(tmp1[i,]*zid_zjd[[i]])}) ))#d x N -> N x d
  
  #return
  if (est == "Y") {
    return(c(sum(Y - lambda),
         tmp3,
         da,
         db))
  }
  
  if (est == "MAP" ) {
    
    return(c(-B/beta.var + sum(Y - lambda),
             
             tmp3 - as.vector(Z)/sigma2_z, # d x N
             
             da - a/sigma2_a,
             
             db - b/sigma2_b,
             
             N*pi*(1/(2*pi*sigma2_a)) + sum(a^2)/(2*sigma2_a^2) +  sender.var.df*prior.sender.var/(2*sigma2_a^2) - (1 + sender.var.df/2)/sigma2_a,
             
             N*pi*(1/(2*pi*sigma2_b)) + sum(b^2)/(2*sigma2_b^2) +  receiver.var.df*prior.receiver.var/(2*sigma2_b^2) - (1 + receiver.var.df/2)/sigma2_b,
             
             N*d*pi*(1/(2*pi*sigma2_z)) + sum(Z^2)/(2*sigma2_z^2) + Z.var.df*prior.Z.var/(2*sigma2_z^2) -
               (1 + Z.var.df/2)/sigma2_z)
             
           )
  }
  
  if (est == "MAPe") {
    return(c(-B/beta.var + sum(Y - lambda),
             
             tmp3 - as.vector(Z)/sigma2_z, # d x N
             
             da - a/sigma2_a,
             
             db - b/sigma2_b)
           )
  }
  
}
# can we pass the same sample to the gradient function?

llik_gr_hat <- function(theta, Y, d, est = "MAP",
                        r = .25, rtype = "constant", Rmin = 5,
                        replace = T, margin = "none", W = NULL ) {
  
  n = nrow(Y)
  
  B = theta[1]
  Z = matrix(theta[2:(1+d*n)], ncol = d, nrow = n)
  a = theta[(2+d*n):(n+1+d*n)]
  b = theta[(n+2+d*n):(2*n+1+d*n)]
  
  if (est == "MAP") {
    sigma2_a = theta[2*n+2+d*n]
    sigma2_b = theta[(2*n+3+d*n)]
    sigma2_z = theta[(2*n+4+d*n)]
  }
  
  #exact update of variance parameters conditioned on graph parameters
  if (est == "MAPe") {
    #optimal values
    sigma2_a = (sum(a^2) + sender.var.df*prior.sender.var^2) / (n + 2 + prior.sender.var)
    sigma2_b = (sum(b^2) + receiver.var.df*prior.receiver.var^2) / (n + 2 + prior.receiver.var)
    sigma2_z = (sum(Z^2) + Z.var.df*prior.Z.var^2) / (n*d + 2 + prior.Z.var)
  }
  
  #add sampling?
  Z_dist = dist2(Z)
  dist_inv = 1/Z_dist; diag(dist_inv) = rep(0, N) #term1 = cbind(term1, term1)
  lambda = exp(t(b + t(a - Z_dist)) + B); diag(lambda) = 0
  tmp1 = (Y + t(Y)) * dist_inv
  tmp2 = dist_inv * (lambda + t(lambda))
  zid_zjd = lapply(1:N, function(x) {t(Z[x,] - t(Z))}) # N, N x d matrices
  tmp3 = as.vector(t( sapply(1:N, function(i) {colSums(tmp2[i,]*zid_zjd[[i]])}) - #d x N -> N x d
                        sapply(1:N, function(i) {colSums(tmp1[i,]*zid_zjd[[i]])}) ))#d x N -> N x d
  
  #return
  if (est == "Y") {
    return(c(sum(Y - lambda),
             tmp3,
             rowSums(Y) - rowSums(lambda),
             colSums(Y) - colSums(lambda)))
  }
  
  if (est == "MAP" ) {
    return(c(-B/beta.var + sum(Y - lambda),
             
             tmp3 - as.vector(Z)/sigma2_z, # d x N
             
             rowSums(Y) - rowSums(lambda) - a/sigma2_a,
             
             colSums(Y) - colSums(lambda) - b/sigma2_b,
             
             N*pi*(1/(2*pi*sigma2_a)) + sum(a^2)/(2*sigma2_a^2) +  sender.var.df*prior.sender.var/(2*sigma2_a^2) - (1 + sender.var.df/2)/sigma2_a,
             
             N*pi*(1/(2*pi*sigma2_b)) + sum(b^2)/(2*sigma2_b^2) +  receiver.var.df*prior.receiver.var/(2*sigma2_b^2) - (1 + receiver.var.df/2)/sigma2_b,
             
             N*d*pi*(1/(2*pi*sigma2_z)) + sum(Z^2)/(2*sigma2_z^2) + Z.var.df*prior.Z.var/(2*sigma2_z^2) - (1 + Z.var.df/2)/sigma2_z)
           
    )
  }
  
  if (est == "MAPe") {
    return(c(-B/beta.var + sum(Y - lambda),
             
             tmp3 - as.vector(Z)/sigma2_z, # d x N
             
             rowSums(Y) - rowSums(lambda) - a/sigma2_a,
             
             colSums(Y) - colSums(lambda) - b/sigma2_b)
    )
  }
  
}

#----------------------------------------------------------------------------------------------------
# pre-calculation (outside the function so not repeated) ----
#Y = as.sociomatrix(movie_net_watch_full, ignore.eval = F, attrname = "ratings_diff")
Y = t(Cmatrix)
#Y = M #note row variance much higher than col variance
c = 1
Yw = Y + c # for Y-based weights, add a constant to Y if it has zeros so no 0 weights
diag(Yw) = 0 # assume no self-loops
lgamma.constant = sum(lgamma(as.vector(Y+1))) # constant for factorial term
W = Yw*(log(Yw)-1) + 5 #add so all weights > 0  
diag(W) = 0
Yw_row = rowSums(Yw)
Yw_col = colSums(Yw)
W_row = rowSums(W)
W_col = colSums(W)

N = n = nrow(Y); d = 2
B = rnorm(1)
Z = rnorm(N*d)
Z_dist = dist2(Y)
#Z = cmdscale(Z_dist, k = d)
#Z = Z/max(abs(Z))
a = rnorm(N)
b = rnorm(N)
sigma2_a = 1
sigma2_b = 1
sigma2_z = 1

beta.var = 9
sender.var.df = 3
receiver.var.df = 3
Z.var.df = sqrt(N)
prior.sender.var = 1
prior.receiver.var = 1
prior.Z.var = N/8

#----------------------------------------------------------------------------------------------------
# run the model
# target for citation with shift(4819, 5061)
# llik(c(latent.srp2.init.res$mkl, sender.var = latent.srp2.init.res$mcmc.mle$sender.var, receiver.var = latent.srp2.init.res$mcmc.mle$receiver.var, Z.var = latent.srp2.init.res$mcmc.mle$Z.var), Y = Y, est= "MAP")
# target for movie_acw with shift (37332.87, 37998.51)
# target for movie_watch based on mle (-565856.5, -569542.9)
# llik(c(latent_acw2_res$mkl, sender.var = latent_acw2_res$mcmc.mle$sender.var, receiver.var = latent_acw2_res$mcmc.mle$receiver.var, Z.var = latent_acw2_res$mcmc.mle$Z.var), Y = Y, est= "MAP")
# llik(latent.watch.mle$start, Y = Y, est = "Y")
#---------------------------------------------------------------------------------------------------
# est Y for just maximizing graph likelihood, MAP for MAP, MAPe for map that uses exact updates of variance 
# paramters conditioned on B,Z,a,b
# For MAPe note variance parameters aren't entered in optim or returned, but calculated in llik2 (and gradient func)
# ---------------------------------------------------------------------------------------------------
# Y ------------------------------------------------------------------------------------------------
#  exact ####
t1 = Sys.time()
theta.bfgs = optim(par = c(B, Z, a, b), 
                   fn = llik2, Y = Y, d = 2, est = "Y", 
                   gr = llik_gr,
                   method = "L-BFGS-B",
                   lower = c(rep(-50, 1 + N*(d+2))),
                   control=list(maxit = 100, fnscale = -1)) #<- need high enough maxit
(Sys.time() - t1)
(Sys.time() - t1)/theta.bfgs$counts[1]
theta.bfgs$counts[1]
theta.bfgs$convergence
llik2(theta.bfgs$par, Y = Y, d = 2, est= "Y")

#  approx - only works with gr function #### 
t1 = Sys.time()
theta.bfgs = optim(par = c(B, Z, a, b), 
                   fn = llik_hat2, Y = Y, d = 2, est = "Y", r = .25, margin = "none",
                   gr = llik_gr_hat,
                   method = "L-BFGS-B",
                   lower = c(rep(-5000, 1 + N*(d+2))),
                   control=list(maxit = 100, fnscale = -1)
                   )
Sys.time() - t1
(Sys.time() - t1)/theta.bfgs$counts[1]
theta.bfgs$counts[1]
theta.bfgs$convergence
llik2(theta.bfgs$par, Y = Y, d = 2, est= "Y")

#plots
model = latent.watch.mle$mle #latent_acw2$mkl #latent.srp2.init.r$mkl #
col1 = as.numeric(as.factor(movie_net_acw_full%v%"genres"))#col2
plot(model$receiver - model$sender,
     theta.bfgs$par[(n+2+d*n):(2*n+1+d*n)] - theta.bfgs$par[(2+d*n):(n+1+d*n)], cex = log(colSums(Y)+2)/5)
cor(model$receiver - model$sender,
     theta.bfgs$par[(n+2+d*n):(2*n+1+d*n)] - theta.bfgs$par[(2+d*n):(n+1+d*n)])
plot(matrix(theta.bfgs$par[2:(1+d*n)], ncol = d, nrow = n), col = col1, pch = 16, cex = 2)


# ---------------------------------------------------------------------------------------------------
# MAP ------------------------------------------------------------------------------------------------
#  exact ####
t1 = Sys.time()
theta.bfgs = optim(par = c(B, Z, a, b, sigma2_a, sigma2_b, sigma2_z), 
                   fn = llik2, Y = Y, d = 2, est = "MAP", 
                   gr = llik_gr, 
                   method = "L-BFGS-B",
                   lower = c(-50, rep(-Inf, N*(d+2)), 0, 0, 0), #-Inf lb on B -> error if too neg
                   control=list(maxit = 2000, fnscale = -1)) 
(Sys.time() - t1)
(Sys.time() - t1)/theta.bfgs$counts[1]
theta.bfgs$counts[1]
theta.bfgs$convergence #ERROR
llik2(theta.bfgs$par, Y = Y, d = 2, est= "Y")
llik2(theta.bfgs$par, Y = Y, d = 2, est= "MAP")

#plots
plot(model$receiver - model$sender,
     theta.bfgs$par[(n+2+d*n):(2*n+1+d*n)] - theta.bfgs$par[(2+d*n):(n+1+d*n)])
cor(model$receiver - model$sender,
    theta.bfgs$par[(n+2+d*n):(2*n+1+d*n)] - theta.bfgs$par[(2+d*n):(n+1+d*n)])
plot(matrix(theta.bfgs$par[2:(1+d*n)], ncol = d, nrow = n), col = col1, pch = 16, cex = 2)

# ---------------------------------------------------------------------------------------------------
# MAPe (better, i think) ------------------------------------------------------------------------------------------------
#  exact ####
t1 = Sys.time()
theta.bfgs = optim(par = c(B, Z, a, b), 
                   fn = llik2, Y = Y, d = 2, est = "MAPe", 
                   gr = llik_gr, 
                   method = "L-BFGS-B",
                   lower = c(rep(-Inf, 1 + N*(d+2))),
                   control=list(maxit = 100, fnscale = -1))
Sys.time() - t1
(Sys.time() - t1)/theta.bfgs$counts[1]
theta.bfgs$counts
theta.bfgs$convergence
llik2(theta.bfgs$par, Y = Y, d = 2, est= "Y")
llik2(theta.bfgs$par, Y = Y, d = 2, est= "MAPe")

#  approx ####
t1 = Sys.time()
theta.bfgs = optim(par = c(B, Z, a, b), 
                   fn = llik_hat2, Y = Y, d = 2, est = "MAPe", r = .25, margin = "none",
                   gr = llik_gr_hat, 
                   method = "L-BFGS-B",
                   lower = c(rep(-Inf, 1 + N*(d+2))),
                   control=list(maxit = 2000, fnscale = -1))
Sys.time() - t1
(Sys.time() - t1)/theta.bfgs$counts[1]
theta.bfgs$counts
theta.bfgs$convergence
llik2(theta.bfgs$par, Y = Y, d = 2, est= "Y")
llik2(theta.bfgs$par, Y = Y, d = 2, est= "MAPe")

#check
B1 = theta.bfgs$par[1]
Z1 = matrix(theta.bfgs$par[2:(1+d*n)], ncol = d, nrow = n)
a1 = theta.bfgs$par[(2+d*n):(n+1+d*n)]
b1 = theta.bfgs$par[(n+2+d*n):(2*n+1+d*n)]
sigma2_a = (sum(a1^2) + sender.var.df*prior.sender.var^2) / (n + 2 + prior.sender.var)
sigma2_b = (sum(b1^2) + receiver.var.df*prior.receiver.var^2) / (n + 2 + prior.receiver.var)
sigma2_z = (sum(Z1^2) + Z.var.df*prior.Z.var^2) / (n*d + 2 + prior.Z.var)

#plots
plot(model$receiver - model$sender, b1 - a1, cex = log(rowSums(Y)+5)/5)
cor(model$receiver - model$sender, b1 - a1)
plot(Z1, col = col1, pch = 16, cex = 2)


# ------------------------------------------------------------------------------------------------
# 4. Multi-dimensional sender and receiver
# Y ------------------------------------------------------------------------------------------------
#  exact ####
t1 = Sys.time()
theta.bfgs = optim(par = c(B, Z, a, b), 
                   fn = llik2, Y = Y, d = 2, est = "Y", 
                   #gr = llik_gr,
                   method = "L-BFGS-B",
                   lower = c(rep(-50, 1 + N*(d+2))),
                   control=list(maxit = 1000, fnscale = -1)) #<- need high enough maxit
(Sys.time() - t1)
(Sys.time() - t1)/theta.bfgs$counts[1]
theta.bfgs$counts[1]
theta.bfgs$convergence
llik2(theta.bfgs$par, Y = Y, d = 2, est= "Y")

theta = theta.bfgs$par
nparam = length(theta)
a = theta[(2+d*n): (1+d*n/2 + (nparam-1)/2) ]
b = theta[(2+d*n/2 + (nparam-1)/2) : nparam]
plot(b[1:47], b[48:94], type = "n")
text(b[1:47], b[48:94], rownames(Cmatrix), cex = 1, col = col2)
#plot(a[1:47], a[48:94], type = "n")
plot(a[1:47], a[48:94], type = "n")
text(a[1:47], a[48:94], rownames(Cmatrix), col= col2, cex = .7)
plot(theta[2:48], theta[49:95], type = "n")
text(theta[2:48], theta[49:95], rownames(Cmatrix), col= col2, cex = 1)

plot(b[1:47] - a[1:47], latent.srp2.init.r$mkl$receiver - latent.srp2.init.r$mkl$sender)
plot(b[48:94] - a[48:94], latent.srp2.init.r$mkl$receiver - latent.srp2.init.r$mkl$sender)
plot(b[1:47]+b[48:94] - a[1:47] - a[48:94], latent.srp2.init.r$mkl$receiver - latent.srp2.init.r$mkl$sender)

#plots
model = latent.watch.mle$mle #latent_acw2$mkl #latent.srp2.init.r$mkl #
col1 = as.numeric(as.factor(movie_net_acw_full%v%"genres"))#col2
plot(model$receiver - model$sender,
     theta.bfgs$par[(n+2+d*n):(2*n+1+d*n)] - theta.bfgs$par[(2+d*n):(n+1+d*n)], cex = log(colSums(Y)+2)/5)
cor(model$receiver - model$sender,
    theta.bfgs$par[(n+2+d*n):(2*n+1+d*n)] - theta.bfgs$par[(2+d*n):(n+1+d*n)])
plot(matrix(theta.bfgs$par[2:(1+d*n)], ncol = d, nrow = n), col = col1, pch = 16, cex = 2)



# ------------------------------------------------------------------------------------------------
# 5. NOTES --------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------
theta.sann = optim(c(B, Z, a, b),
                   Y = Y, d = 2, llik_hat2, method = "SANN",
                   lower = c(rep(-Inf, 1 + N*(d+2)), 0, 0, 0),
                   control=list(maxit = 200, fnscale = -1))

theta.cg = optim(c(a, b, B, Z, sigma2_a, sigma2_b, sigma2_z),
                 Y = Y, d = 2, llik2, method = "CG",
                 lower = c(rep(-Inf, 1 + N*(d+2)), 0, 0, 0),
                 control=list(maxit = 200, fnscale = -1))

theta.bfgs = optim(c(a, b, B, Z, sigma2_a, sigma2_b, sigma2_z), #53 seconds
                   Y = t(Cmatrix), d = 2, llik2, method = "BFGS",
                   lower = c(rep(-Inf, 1 + N*(d+2)), 0, 0, 0),
                   control=list(maxit = 200, fnscale = -1))

theta = unlist(latent.srp2.init.res$mkl)
llik2(theta, Y = Y, d = 2, est = "Y")
llik2(theta.bfgs$par, Y = Y, d = 2, est = "Y")
llik_hat2(theta.bfgs$par, Y = Y, d = 2, r = 1, margin = "none", est = "Y")
llik_hat2(theta.sann$par, Y = Y, d = 2, r = .25, margin = "none", est = "Y") 
llik_hat2(theta.cg$par, Y = Y, d = 2) 


