# quick quasi-Newton for Latent Space Ranking Model
#
# make sure diag of Y is 0 if not counting self-edges
#
# Implementation notes/questions:
#    better to update after each z or all z's at once? (seems better all z's at once, i.e. simultaneous updates, they all move from old z location together)
#    stata j sometimes sneaking away. probably because it has far lowest citation count
#    likelihood levels off after certain #runs, why? yet in terms of correlation to mcmc value fit may be getting better
#    improve ad-hoc tuning/Z search to be more dynamic for any input data set
#    Bayesian priors mean we shouldn't have to center paramters, but could it speed fitting?
#      without bayesian, non-identifiability between sender, receiver, beta (+- constant)
#    how sensitivite to hyperpriors? 
#    note behavior of low connectivity nodes, e.g. if something only sends to 1 other things in network,
#       positions aren't reliable, can tail away from rest of nodes and make visualization worse
#       it's sender and/or receiver coef is conflated with position#
# To DO: 
#    Allow for different variances for z's in different dimensions? Seems realistic.
#    Long term: write the *binomial* quasi-Newton algo

# library ####
library(ergm)
library(latentnet)
library(geoR) # for divnchisq
library(gtools) #for logit function

# log-likelihood function, off by a constant from: ####
#    if est is "MAP" return l(Y|theta) + l(theta) ; if "Y" return l(Y|theta) ; if "theta" return l(theta)
#    if object has parameter values they will be used (not including hyperparameters)
#    family = poisson fits the normal poisson latent space model
#    family = binomial fits a quasi-Stiglery (quasi-symmetric) type model with positions
#    [removed] family = poisson.c fits a constricted poisson position model (T_ij = ~Y_ij + ~Y_ji) by fitting upper tri and letting lower tri be the remainder. Note this is odd because it doesn't demand that E(Cij) = T - E(Cji) given the parameters (positions + random effects + intercept), just that E(C_ji) is a reminder not a glm prediction. This also means the last sender param has no info.  
#              better way to do this with penalty for missing the total instread of sharp constraint (which probably can't be realized )]

llik <- function(object=NULL, Y=NULL, sender=NULL, receiver=NULL, beta=NULL,
                 Z=NULL, sender.var = 10, receiver.var = 10, Z.var = 10,
                 beta.var = 9, sender.var.df = 3, receiver.var.df = 3, Z.var.df = NULL, #N = number of nodes
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
       l_lambda = t(receiver + t(sender - Z_dist)) + beta; 
    
       if (family == "poisson") {
          lambda = exp(l_lambda); diag(lambda) = 0
          pY = sum( Y * l_lambda - lambda)
       }
       
       # if (family == "poisson.c") {
       #   lambda = exp(l_lambda); diag(lambda) = 0
       #   Tmatrix = Y + t(Y); diag(Tmatrix) = 0
       #   lambda[lambda > Tmatrix] = Tmatrix[lambda > Tmatrix]
       #   lambda[lower.tri(lambda)] = (Tmatrix - t(lambda))[lower.tri(lambda)]
       #   l_lambda = log(lambda); diag(l_lambda) = NA
       #   pY = sum( Y * l_lambda - lambda, na.rm = T)
       # }
       
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


# quasi-Netwon ####

lsqn <- function(Y, N=nrow(Y), D = 2, runs = 10, tol = .01, #Y is graph, N = number of nodes
                 # hyperparameters - using defaults from ergmm
                 v_a = 3, v_b = 3, v_z = sqrt(N),
                 s2_a = 1, s2_b = 1, s2_z = N/8,
                 sigma2_B = 9,
                 # prior initial values (will be updated)
                 Z.init = "MDS",
                 RE.init = "rnorm",
                 Z.init.user = NULL, 
                 sigma2_a = 10, sigma2_b = 10, sigma2_z = 10,
                 stepsize.init.a = 1, stepsize.init.b = 1, 
                 stepsize.init.B = 1, stepsize.init.z = 1,
                 noSelfEdges = 1,
                 jmax = 50, epsilon = 1e-10) {
    
    r = 0
    maxllik = 0
    while (r <= runs) {
    while (r == 0) {
  
      # initialize postions (Z) ####
      # multidimenstional scaling (MDS), random normal, or user-specified
      if(Z.init == "rnorm") {
        Z = matrix(rnorm(D*N), ncol = 2)
      } else if (Z.init == "user") {
        Z = Z.init.user
      } else { #MDS is default initializer
        Z_dist = as.matrix(dist(Y))
        Z = cmdscale(Z_dist, k = D)
      }
      
      # Standardize Z to center at origin, [0,1] range
      Z = scale(Z, scale = F) #translate to origin
      Z = Z/max(abs(Z)) # shrink - improve this?
      z = t(Z)
      Z_dist = as.matrix(dist(Z), upper = T)
      
      # initialize a, b, B ####
      a = rnorm(N) #rep(0, N) # #sender
      b = rnorm(N) #rep(0, N) # #receiver
      B = 0 #intercept
      
      # latentnet type initialization
      if (RE.init == "latentnet") {
        a = logit( (rowSums(Y!=0) + 1)/(N-1+2) ) - (1/N) * sum(logit( (rowSums(Y!=0) + 1) / (N - 1 + 2)) )
        b = logit( (colSums(Y!=0) + 1)/(N-1+2) ) - (1/N) * sum(logit( (colSums(Y!=0) + 1) / (N - 1 + 2)) )
        sigma2_a = var(a)
        sigma2_b = var(b)
        B = ( 1/(N*(N-1)) * sum(Y>mean(Y)) + mean(Z_dist))
        sigma2_z = var(as.vector(z))
      }
      
      r = r+1
  
    }

    # Update  Z, sigma2_z, B, a, sigma2_a, b, sigma2_b  ####
    # sigma updates are closed form. Others by coordinate ascent

    # Z ####
  
    # - init ####
    Z_dist = as.matrix(dist(Z, upper = T)) #just in case
    stepsize.z = matrix(stepsize.init.z, D, N)
    zid_zjd = lapply(1:N, function(x) {t(Z[x,] - z)}) # N, N x d matrices
    dist_inv = 1/Z_dist; diag(dist_inv) = rep(0, N) #term1 = cbind(term1, term1)
    yij_yji = Y + t(Y) # each entry is y_ij + y_ji
      tmp1 = yij_yji * dist_inv
    exp_term = exp(sweep(t(b - Z_dist + B), 1, a, "+"));  # = t(b + t(a-Z_dist) +B)
    exp_term = exp_term + t(exp_term) # each entry is ij + ji
    diag(exp_term) = 0
      tmp2 = dist_inv * exp_term
    
    #first deriv wrt z_id:
    diff_z =  sapply(1:N, function(i) {colSums(tmp2[i,]*zid_zjd[[i]])}) - #d x N
              sapply(1:N, function(i) {colSums(tmp1[i,]*zid_zjd[[i]])}) - #d x N
              z/sigma2_z # d x N
      
    zsign = sign(diff_z)   
    
    #second deriv wrt z_id:
    tmp3 = dist_inv^(3) * (yij_yji  - exp_term*(1 + Z_dist))
    diff_z2 = t(t(sapply(1:N, function(i) {colSums(tmp3[i,]*(zid_zjd[[i]])^2)})) + #N x d
              colSums(dist_inv * (-yij_yji + exp_term)) #N, added to above by column
              - 1/sigma2_z)
      
    zsign2 = sign(diff_z2)   
      
    # - update ####
    for (i in sample(N)) { #update one dimension at a time? randomize order? fix first point to prevent translations?
          #print(i)
          for (d in sample(D)) {
            
                j = 0
                while(abs(diff_z[d,i]) > tol && j <= jmax) { #
                br = FALSE
                
                #consider a new zid:
                znew_i = z[,i]; znew_i[d]= znew_i[d] + stepsize.z[d,i]*zsign[d,i]*-zsign2[d,i]
                tmpdist = sqrt(colSums((znew_i - t(Z))^2)); if (noSelfEdges) {tmpdist[i] = 1}
                
                #recalculate zdiff for that zid
                diff_z[d,i] = - sum(yij_yji[i,] * 1/tmpdist * (znew_i[d] - Z[,d])) +
                              sum(1/tmpdist * (znew_i[d] - Z[,d]) * 
                              (exp(B + a[i] + b - tmpdist) + exp(B + a + b[i] - tmpdist))) -
                              znew_i[d]/sigma2_z #yij_yji is zero if i = j so don't worry about Z not updated to znew on that line
                
                #compre to z orig, did we cross a 0?
                if (sign(diff_z[d,i]) != zsign[d,i]) {
                   s = z[d,i] + seq(0, 1, length.out = 11) * stepsize.z[d,i]*zsign[d,i]*-zsign2[d,i]
                   tmpdiff = sapply(s, function(x) {
                                      znew_i[d]= x
                                      tmpdist = sqrt(colSums((znew_i - t(Z))^2)); if (noSelfEdges) {tmpdist[i] = 1}
                                      return(- sum(yij_yji[i,] * 1/tmpdist * (znew_i[d] - Z[,d])) +
                                            sum(1/tmpdist * (znew_i[d] - Z[,d]) * 
                                            (exp(B + a[i] + b - tmpdist) + exp(B + a + b[i] - tmpdist))) -
                                            znew_i[d]/sigma2_z)})
                   #update lower case z (not upper case Z yet)
                   z[d,i] = s[which.max(sign(tmpdiff)!=zsign[d,i])]
                   #Z[i,d] = z[d,i]
                   stepsize.z[d,i] = stepsize.z[d,i]/10
                   if (stepsize.z[d,i] < epsilon) {br = TRUE}
                   diff_z[d,i] = tmpdiff[which.max(sign(tmpdiff)!=zsign[d,i])]
                   zsign[d,i] = sign(diff_z[d,i]) 
                      znew_i = z[,i]
                      tmpdist = sqrt(colSums((znew_i - t(Z))^2)); if (noSelfEdges) {tmpdist[i] = 1}
                      tmpexp = exp(B + a[i] + b - tmpdist) + exp(B + a + b[i] - tmpdist)
                      tmpexp[i] = 0
                      tmp4 = 1/tmpdist^(3) * (yij_yji[i,]  - tmpexp*(1 + tmpdist))
                  diff_z2[d,i] = sum(tmp4*(znew_i[d] - t(Z)[d,])^2) +
                               sum(1/tmpdist * (-yij_yji[i,] + tmpexp)) - 1/sigma2_z
                  zsign2[d,i] = sign(diff_z2[d,i]) 
                } else {
                    if (stepsize.z[d,i] <= max(abs(z))) { stepsize.z[d,i] = stepsize.z[d,i] * 2
                    } else {
                        #print("had to look in both directions")
                        #look in both directions
                        s = z[d,i] + seq(-1, 1, length.out = 101) * stepsize.z[d,i]
                        s = s[-51]
                        tmpdiff = sapply(s, function(x) {
                                      znew_i[d]= x
                                      tmpdist = sqrt(colSums((znew_i - t(Z))^2)); if (noSelfEdges) {tmpdist[i] = 1}
                                      return(- sum(yij_yji[i,] * 1/tmpdist * (znew_i[d] - Z[,d])) +
                                            sum(1/tmpdist * (znew_i[d] - Z[,d]) * 
                                            (exp(B + a[i] + b - tmpdist) + exp(B + a + b[i] - tmpdist))) -
                                            znew_i[d]/sigma2_z)})
                        
                        stepsize.z[d,i] = abs(s[which.min(abs(tmpdiff))] - z[d,i]) #don't let it pick itself
                        z[d,i] = s[which.min(abs(tmpdiff))] #update lower case z (not upper case Z yet)
                        diff_z[d,i] = tmpdiff[which.min(abs(tmpdiff))]
                        zsign[d,i] = sign(diff_z[d,i]) 
                              znew_i = z[,i]
                              tmpdist = sqrt(colSums((znew_i - t(Z))^2)); if (noSelfEdges) {tmpdist[i] = 1}
                              tmpexp = exp(B + a[i] + b - tmpdist) + exp(B + a + b[i] - tmpdist)
                              tmpexp[i] = 0
                              tmp4 = 1/tmpdist^(3) * (yij_yji[i,]  - tmpexp*(1 + tmpdist))
                        diff_z2[d,i] = sum(tmp4*(znew_i[d] - t(Z)[d,])^2) +
                            sum(1/tmpdist * (-yij_yji[i,] + tmpexp)) - 1/sigma2_z
                        zsign2[d,i] = sign(diff_z2[d,i]) 
                    }
                  }
                if(br) {break}
                j = j+1}
          }
    }
    
    z  = z - rowMeans(z)
    Z = t(z)
    Z_dist = as.matrix(dist(Z, upper = T))
  
    # sigma2_z ####
    sigma2_z = (sum(z^2) + v_z*s2_z^2) / (N*d + 2 + v_z) #z has length N*d
  
    #B: ####
    stepsize.B = stepsize.init.B
    lambdamat = exp(sweep(t(b - Z_dist + B), 1, a, "+"))
    diff_B = sum( Y - lambdamat) + sum(diag(lambdamat))*noSelfEdges - B/sigma2_B #in
    #second deriv always negative -> concave
    Bsign = sign(diff_B)
    
    while (abs(diff_B) > tol) {
      Bnew = B + stepsize.B*Bsign
      lambdamat = exp(sweep(t(b - Z_dist + Bnew), 1, a, "+"))
      diff_B = sum( Y - lambdamat) + sum(diag(lambdamat))*noSelfEdges - Bnew/sigma2_B
      if (sign(diff_B) != Bsign) { #look in this range
        s = B + seq(0, 1, length.out = 11) * stepsize.B*Bsign
        tmp = sapply(s, function(B) {
          lambdamat = exp(sweep(t(b - Z_dist + B), 1, a, "+"))
          sum( Y - lambdamat) + sum(diag(lambdamat))*noSelfEdges - B/sigma2_B
        })
        B = s[which.max(sign(tmp)!=Bsign)]
        stepsize.B = stepsize.B/10
        diff_B = tmp[which.max(sign(tmp)!=Bsign)]
        Bsign = sign(diff_B)
      } else {stepsize.B = stepsize.B * 2}
    }
    
    # a ####
    #init
    stepsize.a = rep(stepsize.init.a, N)
    lambdamat = exp(sweep(t(b - Z_dist + B), 1, a, "+"))
    diff_a = rowSums(Y) - rowSums(lambdamat) + diag(lambdamat) - a/sigma2_a #i,j entry is i to j
    asign = sign(diff_a)
    
    #go
    while (max(abs(diff_a)) > tol) {
      #anew = a + stepsize.a*asign
      #lambdamat = exp(sweep(t(b - Z_dist + B), 1, anew, "+"))
      #diff_a = rowSums(Y) - rowSums(lambdamat) + diag(lambdamat) - anew/sigma2_a #i,j entry is i to 
      for (i in 1:N) {
        while (abs(diff_a[i]) > tol) {
          anew = a[i] + stepsize.a[i]*asign[i]
          lambdavec = exp(B + anew + b - Z_dist[i,])
          lambdavec[i] = lambdavec[i] - lambdavec[i]*noSelfEdges #if necessary, remove self edge
          diff_a[i] = sum( Y[i,]) -
            sum(lambdavec) - 
            anew/sigma2_a
          if (sign(diff_a[i]) != asign[i]) {
            s = a[i] + seq(0, 1, length.out = 11) * stepsize.a[i]*asign[i]
            if (noSelfEdges) {
              tmp = sum( Y[i,]) - 
                sapply(s, function(a) {sum(exp(B + a + b - Z_dist[i,])[-i]) + a/sigma2_a})
            } else { 
              tmp = sum( Y[i,]) - 
                sapply(s, function(a) {sum(exp(B + a + b - Z_dist[i,])) + a/sigma2_a})  
            }
            a[i] = s[which.max(sign(tmp)!=asign[i])]
            stepsize.a[i] = stepsize.a[i]/10
            diff_a[i] = tmp[which.max(sign(tmp)!=asign[i])]
            asign[i] = sign(diff_a[i])
          } else {stepsize.a[i] = stepsize.a[i] * 2}
        }
      }
    }
    
    #sigma2_a ####
    sigma2_a = (sum(a^2) + v_a*s2_a^2) / (N + 2 + v_a)
    
    # b ####
    #init
    stepsize.b = rep(stepsize.init.b, N)
    lambdamat = exp(sweep(t(b - Z_dist + B), 1, a, "+"))
    diff_b = colSums(Y) - colSums(lambdamat) + diag(lambdamat) - b/sigma2_b #i,j entry is i to j
    bsign = sign(diff_b)
    
    #go
    while (max(abs(diff_b)) > tol) {
      #bnew = b + stepsize.b*bsign
      #lambdamat = exp(sweep(t(bnew - Z_dist + B), 1, a, "+"))
      #diff_b = colSums(Y) - colSums(lambdamat) + diag(lambdamat) - bnew/sigma2_b #i,j entry is i to 
      for (i in 1:N) {
        while (abs(diff_b[i]) > tol) {
          bnew = b[i] + stepsize.b[i]*bsign[i]
          lambdavec = exp(B + bnew + a - Z_dist[,i])
          lambdavec[i] = lambdavec[i] - lambdavec[i]*noSelfEdges #if necessary, remove self edge
          diff_b[i] = sum( Y[,i]) -
            sum(lambdavec) - 
            bnew/sigma2_b
          if (sign(diff_b[i]) != bsign[i]) {
            s = b[i] + seq(0, 1, length.out = 11) * stepsize.b[i]*bsign[i]
            if (noSelfEdges) {
              tmp = sum( Y[,i]) - 
                sapply(s, function(b) {sum(exp(B + a + b - Z_dist[,i])[-i]) + b/sigma2_b})
            } else { 
              tmp = sum( Y[,i]) - 
                sapply(s, function(a) {sum(exp(B + a + b - Z_dist[,i])) + b/sigma2_b})  
            }
            b[i] = s[which.max(sign(tmp)!=bsign[i])]
            stepsize.b[i] = stepsize.b[i]/10
            diff_b[i] = tmp[which.max(sign(tmp)!=bsign[i])]
            bsign[i] = sign(diff_b[i])
          } else {stepsize.b[i] = stepsize.b[i] * 2}
        }
      }
    }
    
    # sigma2_b ####
    sigma2_b = (sum(b^2) + v_b*s2_b^2) / (N + 2 + v_b)
    
    # likelihood ####
    currentllik = llik(Y=Y, sender = a, receiver = b, beta = B, Z = t(z), sender.var = sigma2_a,
                       receiver.var = sigma2_b, Z.var = sigma2_z, beta.var = sigma2_B)
    print(currentllik)
    if (currentllik  > maxllik) {
        maxllik = currentllik
        map = list(Z = scale(Z, scale =F), sender = a, receiver = b, beta = B,
        sender.var = sigma2_a, receiver.var = sigma2_b, Z.var = sigma2_z,  #or just recalculate these given other params
        diff_Z = diff_z, diff_sender = diff_a, diff_receiver = diff_b, diff_Z2 = diff_z2,
        stepsize.Z = stepsize.z, llik = maxllik)
    }
    
    # ending tasks ####
    r = r+1
    
}         #shift positions to origin mean before returning
    return(list(last = list(Z = scale(Z, scale = F), sender = a, receiver = b, beta = B,
                sender.var = sigma2_a, receiver.var = sigma2_b, Z.var = sigma2_z,
                diff_Z = diff_z, diff_sender = diff_a, diff_receiver = diff_b, diff_Z2 = diff_z2,
                stepsize.Z = stepsize.z, llik = currentllik),
                map = map))
}


# predict network based on quasi-Newton fit
# either MAP point estimate or random draw|MAP
predict.lsqn <-function(model, type = "MAP", names = NULL) {
  N = nrow(model$map$Z)
  Z_dist = as.matrix(dist(model$map$Z, upper = TRUE), N, N)
  lambda = exp(t(model$map$receiver + t(model$map$sender - Z_dist)) + model$map$beta); diag(lambda) = NA
  if (!is.null(names)) {row.names(lambda) = names}
  if (type == "MAP") return(lambda) else {
    if (type == "rpois") {
      if (!is.null(names)) {row.names(lamdda) = names}
        return(matrix(rpois(N^2, lambda), N, N))
    }
      
  }
  
}


#BIC estimation? ####
samp1 = rep(0,1000)
N = 47
for (i in 1:length(samp1)) {
    sender.var = rinvchisq(1, sender.var.df, prior.sender.var)
    receiver.var = rinvchisq(1, receiver.var.df, prior.receiver.var)
    Z.var = rinvchisq(1, sqrt(N), N/8)
  object = list(
    beta = rnorm(1, 0, sd = sqrt(beta.var)),
    sender = rnorm(N, 0, sqrt(sender.var)),
    receiver = rnorm(N, 0, sqrt(receiver.var)),
    Z = matrix(rnorm(N*d, 0, sqrt(Z.var)), nrow = N, ncol = d)
  )
  samp1[i] = llik(object, Y, est = "Y")
}

#DIC estimation? ####

-2* ((mean(latent.srp0$sample$lpY) + lgamma.constant) -
       llik(latent.srp0$mkl, Z = matrix(0, 47, 2), Y, est = "Y"))
       #(latnet.srp0$mcmc.mle$lpY + lgamma.constant)

-2* ((mean(latent.srp2$sample$lpY) + lgamma.constant) -
       llik(latent.srp2$mkl, Y, est = "Y"))
       #(latnet.srp2$mcmc.mle$lpY + lgamma.constant)
