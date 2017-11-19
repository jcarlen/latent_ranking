# Helper likelihood function for latent space rating
#
#    Diag of network matrix representation should be 0 if not counting self-edges (!)
#
#    if est is "MAP" return l(Y|theta) + l(theta) ; if "Y" return l(Y|theta) ; if "theta" return l(theta)
#    if 'object' has parameter values they will be used (not including hyperparameters)
#    family = poisson fits the poisson latent space model
#    family = binomial fits a quasi-Stiglery (quasi-symmetric) type model with positions -- off by a constant
#    family = poisson_md fits the poisson latent space model with multi-dimensional random sender and receiver effects
#             assume sr dimension same as Z dimension
#
#    for approximate likelihood schemes see llik_approx.R
#
# --------------------------------------------------------------------------------------------------

library(extraDistr) # for divnchisq
library(gtools) #for logit function

#----------------------------------------------------------------------------------------------------
# 1. likelihood: llik, llik2 ----------------------------------------------------------------------------------------------------
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
                 est = "MAP", family = "poisson",
                 u = NULL, v = NULL,
                 u.var = 10, v.var = 10,
                 u.var.df = 3, v.var.df = 3,
                 prior.u.var = 1, prior.v.var = 1
                 ) {
  
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

  if (family == "poisson") {
    #l_lambda = t(receiver + t(sender - Z_dist)) + beta;
    l_lambda = beta + outer(sender, receiver, "+") - Z_dist #slightly faster
    lgamma.constant = sum(lgamma(as.vector(Y+1)), na.rm = T) 
    lambda = exp(l_lambda); diag(lambda) = 0
    pY = sum( Y * l_lambda - lambda, na.rm = T) - lgamma.constant
  }

  if (family == "binomial") {
    #l_lambda = t(receiver + t(sender - Z_dist)) + beta;
    l_lambda = beta + outer(sender, receiver, "+") + - Z_dist #slightly faster
    lambda = inv.logit(l_lambda)
    Yt = Y + t(Y); diag(Yt) =  0
    pY =  sum( Y * log(lambda), na.rm = T) + sum((Yt - Y)*log(1-lambda), na.rm = T)
  }   

  if (family == "poisson_l") { 
    if (ncol(u) == 1) {u = c(u); v = c(v)} #for outer
    l_lambda = beta + outer(sender, receiver, "+") + outer(u, v, "+")/(Z_dist+1) - Z_dist
    lambda = exp(l_lambda); diag(lambda) = 0
    pY = sum( Y * l_lambda - lambda, na.rm = T) - lgamma.constant
  }
  
  #poisson, multiplicative random effects
  if (family == "poisson_m") { 
    l_lambda = beta + outer(sender, receiver, "+") + u %*% t(v) - Z_dist
    lambda = exp(l_lambda); diag(lambda) = 0
    pY = sum( Y * l_lambda - lambda, na.rm = T) - lgamma.constant
  }
  
  #poisson, multiplicative, project
  #if (family == "poisson_m_proj") { }
  
  #poisson, additive, project
  #if (family == "poisson_a_proj") { #multidimensional projection model
  #  d = ncol(Z) #model assumes rank = dimension
  #  sr = lapply(1:d, function(x) {outer(c(u[,x]), c(v[,x]), "+") } )
  #  zj_zi = lapply(1:d, function(x) { t(outer(Z[,x], Z[,x], "-")) } )
  #  sr_zj_zi = lapply(1:d, function(x) {sr[[x]] * zj_zi[[x]]})
  #  l_lamdba = beta + Reduce('+', sr_zj_zi)/Z_dist - Z_dist
  #  lambda = exp(l_lambda); diag(lambda) = 0
  #  pY = sum( Y * l_lambda - lambda, na.rm = T) - lgamma.constant
  #}

  if (est == "Y") {return(pY)}
  
  ptheta = log(exp(-beta^2/(2*beta.var)) / sqrt(2*pi*beta.var)) +
    sum(log(exp(-sender^2/(2*sender.var)) / sqrt(2*pi*sender.var))) + 
    sum(log(exp(-receiver^2/(2*receiver.var)) / sqrt(2*pi*receiver.var))) +
    sum(log(exp(-Z^2/(2*Z.var)) / sqrt(2*pi*Z.var))) +
    log(dinvchisq(sender.var, sender.var.df, prior.sender.var)) + 
    log(dinvchisq(receiver.var, receiver.var.df, prior.receiver.var)) + 
    log(dinvchisq(Z.var, Z.var.df, prior.Z.var))
  
  #if we have u,v terms:
  if (length(grep("poisson_", family)) > 0) {
    ptheta = ptheta + 
             sum(log(exp(-u^2/(2*u.var)) / sqrt(2*pi*u.var))) + 
             sum(log(exp(-v^2/(2*v.var)) / sqrt(2*pi*v.var))) +   
             log(dinvchisq(u.var, u.var.df, prior.u.var)) + 
             log(dinvchisq(v.var, v.var.df, prior.v.var))
  }
  
  if (est == "theta") {return(ptheta)}
  
  map = pY + ptheta # = p(Y|theta) + p(theta)
  
  if (est == "MAP") {return(map)}
  
}

#for using llik with optim:
llik2 <- function(theta, Y, d, R = 1, est = "MAP", family = "poisson") {
  
  n = nrow(Y)
  
  B = theta[1]
  Z = matrix(theta[2:(1+d*n)], ncol = d, nrow = n)
  
  a = theta[(2+d*n) : (1 + n+d*n)]
  b = theta[(n+2+d*n) : (1+(2+d)*n)] 
  u = v = sigma2_u = sigma2_v = NULL
  
  if (length(grep("poisson_", family)) > 0) {
    u = matrix(theta[(2+(2+d)*n) : (1+(2+d+R)*n)], ncol = R)
    v = matrix(theta[(2+(2+d+R)*n): (1+(2+d+2*R)*n)], ncol = R)
  }
  
  if (est == "MAP") {
    nparam = length(theta)
    if (length(grep("poisson_", family)) > 0) {
      sigma2_a = theta[nparam-4]
      sigma2_b = theta[nparam-3]
      sigma2_z = theta[nparam-2]
      sigma2_u = theta[nparam-1]
      sigma2_v = theta[nparam]
    } else {
      sigma2_a = theta[nparam-2]
      sigma2_b = theta[nparam-1]
      sigma2_z = theta[nparam]
    }
  }
  
  #exact update conditioned on graph parameters
  if (est == "MAPe") {
    sigma2_a = (sum(a^2) + sender.var.df*prior.sender.var^2) / (n + 2 + prior.sender.var)
    sigma2_b = (sum(b^2) + receiver.var.df*prior.receiver.var^2) / (n + 2 + prior.receiver.var)
    sigma2_z = (sum(Z^2) + Z.var.df*prior.Z.var^2) / (n*d + 2 + prior.Z.var)
    if (length(grep("poisson_", family)) > 0) {
      sigma2_u = (sum(u^2) + u.var.df*prior.sender.var^2) / (n*R + 2 + prior.u.var)
      sigma2_v = (sum(v^2) + v.var.df*prior.v.var^2) / (n*R + 2 + prior.v.var)
    }
  } 

  if (est == "MAPe") {est = "MAP"} # to evaluate likelihood correctly
  
  return(llik(Y=Y, sender = a, receiver = b, beta = B, Z = Z,
              sender.var = sigma2_a, receiver.var = sigma2_b,
              Z.var = sigma2_z, family = family, est = est, 
              u = u, v = v, u.var = sigma2_u, v.var = sigma2_v))
}

#----------------------------------------------------------------------------------------------------
# 2. gradient: llik_gr  --------------------------------------------------------------------

llik_gr <- function(theta, Y, d, R = 1, est = "MAP", family = "poisson") {
  
  n = nrow(Y)
  
  B = theta[1]
  Z = matrix(theta[2:(1+d*n)], ncol = d, nrow = n)
  
  a = theta[(2+d*n) : (1 + n+d*n)]
  b = theta[(n+2+d*n) : (1+(2+d)*n)] 
  u = v = sigma2_u = sigma2_v = NULL
  
  if (family == "poisson_l" | family == "poisson_m") {
    u = matrix(theta[(2+(2+d)*n) : (1+(2+d+R)*n)], ncol = R)
    v = matrix(theta[(2+(2+d+R)*n): (1+(2+d+2*R)*n)], ncol = R)
  }
  
  if (est == "MAP") {
    nparam = length(theta)
    if (length(grep("poisson_", family)) > 0) {
      sigma2_a = theta[nparam-4]
      sigma2_b = theta[nparam-3]
      sigma2_z = theta[nparam-2]
      sigma2_u = theta[nparam-1]
      sigma2_v = theta[nparam]
    } else {
      sigma2_a = theta[nparam-2]
      sigma2_b = theta[nparam-1]
      sigma2_z = theta[nparam]
    }
  }
  
  #exact update of variance parameters conditioned on graph parameters
  if (est == "MAPe") {
    #optimal values
    sigma2_a = (sum(a^2) + sender.var.df*prior.sender.var^2) / (n + 2 + prior.sender.var)
    sigma2_b = (sum(b^2) + receiver.var.df*prior.receiver.var^2) / (n + 2 + prior.receiver.var)
    sigma2_z = (sum(Z^2) + Z.var.df*prior.Z.var^2) / (n*d + 2 + prior.Z.var)
    if (length(grep("poisson_", family)) > 0) {
      sigma2_u = (sum(u^2) + u.var.df*prior.sender.var^2) / (n*R + 2 + prior.u.var)
      sigma2_v = (sum(v^2) + v.var.df*prior.v.var^2) / (n*R + 2 + prior.v.var)
    }
  }
  
  Z_dist = dist2(Z)
  dist_inv = 1/Z_dist; diag(dist_inv) = 0
  
  if (family == "poisson") {
    lambda = exp(B + outer(a, b, "+") - Z_dist); diag(lambda) = 0
    Y_l = Y - lambda
    #dz:
    tmp1 = (-Y_l - t(Y_l)) * dist_inv
    zid_zjd = lapply(1:N, function(x) {t(Z[x,] - t(Z))}) 
    dz = as.vector(t( sapply(1:N, function(i) {colSums(tmp1[i,]*zid_zjd[[i]])})))
  }
  
  if (family == "poisson_l") { 
    u = c(u); v = c(v) #for outer
    uv = outer(u, v, "+")
    dist_inv1 = 1/(Z_dist + 1); diag(dist_inv1) = 0
    lambda = exp(B + outer(a, b, "+") + uv*dist_inv1 - Z_dist); diag(lambda) = 0
    Y_l = Y - lambda
    #dz:
    tmp1 = dist_inv * ( (1 + uv * dist_inv1^2)*(Y_l) + (1 + t(uv) * dist_inv1^2)*(t(Y_l)) )
    zid_zjd = lapply(1:N, function(x) {t(Z[x,] - t(Z))})
    dz = as.vector(t( sapply(1:N, function(i) {colSums(tmp1[i,]*-zid_zjd[[i]])}))) 
    tmp = dist_inv1*(Y_l)
    du = rowSums(tmp)
    dv = colSums(tmp)
  }

  if (family == "poisson_m") { 
    lambda = exp(B + outer(a, b, "+") + u %*% t(v) - Z_dist); diag(lambda) = 0
    Y_l = Y - lambda
    #dz:
    tmp1 = (-Y_l - t(Y_l)) * dist_inv 
    zid_zjd = lapply(1:N, function(x) {t(Z[x,] - t(Z))})
    dz = as.vector(t( sapply(1:N, function(i) {colSums(tmp1[i,]*zid_zjd[[i]])})))
    du = as.vector(sapply(1:R, function(i) {Y_l%*%v[,i]})) 
    dv = as.vector(sapply(1:R, function(i) {t(Y_l)%*%u[,i]})) 
  }

  da = rowSums(Y_l)
  db = colSums(Y_l)
  dB = sum(da)
  
  if (est == "MAP" | est == "MAPe") {
    dB = dB - B/beta.var
    dz = dz - as.vector(Z)/sigma2_z
    da = da - a/sigma2_a
    db = db - b/sigma2_b
    if (family == "poisson_l" | family == "poisson_m") {
      du = du - c(u)/sigma2_u
      dv = dv - c(v)/sigma2_v
    }
  }

  if (family == "poisson") {
    #return
    if (est == "Y" | est == "MAPe") {
      return(c(dB, dz, da, db))
    }
  
    if (est == "MAP" ) {
    
     return(c(dB, dz, da, db,
             
             N*pi*(1/(2*pi*sigma2_a)) + sum(a^2)/(2*sigma2_a^2) +  sender.var.df*prior.sender.var/(2*sigma2_a^2) - (1 + sender.var.df/2)/sigma2_a,
             
             N*pi*(1/(2*pi*sigma2_b)) + sum(b^2)/(2*sigma2_b^2) +  receiver.var.df*prior.receiver.var/(2*sigma2_b^2) - (1 + receiver.var.df/2)/sigma2_b,
             
             N*d*pi*(1/(2*pi*sigma2_z)) + sum(Z^2)/(2*sigma2_z^2) + Z.var.df*prior.Z.var/(2*sigma2_z^2) -
               (1 + Z.var.df/2)/sigma2_z
             ))
    }
  }
 
  if (family == "poisson_l" | family == "poisson_m") {
    
    if (est == "Y" | est == "MAPe") {
      return(c(dB, dz, da, db, du, dv))
    }
    
    if (est == "MAP" ) {
      
      return(c(dB, dz, da, db, du, dv,
               
               N*pi*(1/(2*pi*sigma2_a)) + sum(a^2)/(2*sigma2_a^2) +  sender.var.df*prior.sender.var/(2*sigma2_a^2) - (1 + sender.var.df/2)/sigma2_a,
               
               N*pi*(1/(2*pi*sigma2_b)) + sum(b^2)/(2*sigma2_b^2) +  receiver.var.df*prior.receiver.var/(2*sigma2_b^2) - (1 + receiver.var.df/2)/sigma2_b,
               
               N*d*pi*(1/(2*pi*sigma2_z)) + sum(Z^2)/(2*sigma2_z^2) + Z.var.df*prior.Z.var/(2*sigma2_z^2) -
                 (1 + Z.var.df/2)/sigma2_z,
               
               N*R*pi*(1/(2*pi*sigma2_u)) + sum(u^2)/(2*sigma2_u^2) +  u.var.df*prior.u.var/(2*sigma2_u^2) - (1 + u.var.df/2)/sigma2_u,
             
               N*R*pi*(1/(2*pi*sigma2_v)) + sum(v^2)/(2*sigma2_v^2) +  v.var.df*prior.v.var/(2*sigma2_v^2) - (1 + v.var.df/2)/sigma2_v
             ))
    }
  }
    
  }

# -------------------------------------------------------------------------------------------------
# 3. NOTES --------------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------------------
#
# other optimization methods
#
# theta.sann = optim(c(B, Z, a, b),
#                    Y = Y, d = 2, llik_hat2, method = "SANN",
#                    lower = c(rep(-Inf, 1 + N*(d+2)), 0, 0, 0),
#                    control=list(maxit = 200, fnscale = -1))
# 
# theta.cg = optim(c(a, b, B, Z, sigma2_a, sigma2_b, sigma2_z),
#                  Y = Y, d = 2, llik2, method = "CG",
#                  lower = c(rep(-Inf, 1 + N*(d+2)), 0, 0, 0),
#                  control=list(maxit = 200, fnscale = -1))
# 
# theta.bfgs = optim(c(a, b, B, Z, sigma2_a, sigma2_b, sigma2_z), #53 seconds
#                    Y = t(Cmatrix), d = 2, llik2, method = "BFGS",
#                    lower = c(rep(-Inf, 1 + N*(d+2)), 0, 0, 0),
#                    control=list(maxit = 200, fnscale = -1))
# 
# theta = unlist(latent.srp2.init.res$mkl)
# llik2(theta, Y = Y, d = 2, est = "Y")
# llik2(theta.bfgs$par, Y = Y, d = 2, est = "Y")
# llik_hat2(theta.bfgs$par, Y = Y, d = 2, r = 1, margin = "none", est = "Y")
# llik_hat2(theta.sann$par, Y = Y, d = 2, r = .25, margin = "none", est = "Y") 
# llik_hat2(theta.cg$par, Y = Y, d = 2) 


