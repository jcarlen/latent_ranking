# Code for running examples and analyzing results for latent space models
# Comparing MCMC algorithm to quasi-Newton

library(scales) #for plot colors
library(vegan) #for procrustes

# Citation Example ####

Y = t(Cmatrix)
model_cite = lsqn(Y, runs = 50, tol = .01, RE.init = "latentnet")
#best results with these settings ^
#better to do more initializations with fewer runs than vice versa
model_cite2 = lsqn(Y, runs = 50, tol = .0001,
                   Z.init = "user", Z.init.user = model_cite$map$Z) #best results with these settings

# Movie Example ####
  # setup ####
movie = readRDS("~/Documents/citation/latent_ranking/movie_output.RDS")
movie_net_acw = movie$movie_net_acw
movie_net_acw_full = movie$movie_net_acw_full

W = as.sociomatrix(movie_net_watch_full, ignore.eval = F, attrname = "ratings_diff")
#M = as.sociomatrix(movie_net_acw, ignore.eval = F, attrname = "ratings_diff")
M = as.sociomatrix(movie_net_acw_full, ignore.eval = F, attrname = "ratings_diff")

# try whole movie networks, restricted to those with at least one in and out tie, 3664
# Y = movie_rating[(colSums(movie_rating!=0)!=0 & rowSums(movie_rating!=0)!=0),
#                 (colSums(movie_rating!=0)!=0 & rowSums(movie_rating!=0)!=0)]
# this is prohibitively slow in R without parallelization

rownames(M) = movie_net_acw_full%v%"titles"
N = nrow(M)
sigma2_B = 9
v_a = 3; v_b = 3; v_z = sqrt(N)
s2_a = 1; s2_b = 1; s2_z = N/8

  # run ####
model_acw_full = lsqn(M, N = nrow(M), D = 2, runs = 100, tol = .1) #takes a few minutes
#model_acw_full.rnorm = lsqn(M, N = nrow(M), D = 2, runs = 100, Z.init = "rnorm", tol = .1)
model_watch = lsqn(W, N = nrow(W), D = 2, runs = 50, tol = .1) #14:40

  # run mcmc ####

test = ergmm(movie_net_acw_full ~ euclidean(d=2) + rreceiver + rsender, 
             response = "ratings_diff", family = "Poisson.log",
             control = control.ergmm(pilot.runs = 4,  burnin = 50000, 
                                     interval = 500, sample.size = 10000,
                                     optim.method = "CG"),
             seed = 123, verbose = 3,
             tofit = c("mle"),
             user.start = list(sender = a, receiver = b, sender.var = sigma2_a, receiver.var = sigma2_b, 
                               beta = B, Z = Z, Z.var = sigma2_z))

Sys.time() #9:09
latent_acw2 = ergmm(movie_net_acw_full ~ euclidean(d=2) + rreceiver + rsender, 
                    response = "ratings_diff", family = "Poisson.log",
                    control = control.ergmm(pilot.runs = 4,  burnin = 5000000, 
                                            interval = 500, sample.size = 10000), seed = 123, verbose = 3)
Sys.time()
latent_acw2.0 = ergmm(movie_net_acw_full ~ euclidean(d=2) + rreceiver + rsender, 
                      response = "ratings_diff", family = "Poisson.log",
                      control = control.ergmm(pilot.runs = 4,  burnin = 0, 
                                              interval = 500, sample.size = 10000), seed = 123)

Sys.time()
model_acw_full = lsqn(M, N = nrow(M), D = 2, runs = 100, tol = .1) #takes a few minutes
Sys.time()
model_acw_full2 = lsqn(M, N = nrow(M), D = 2, runs = 100, tol = .001) #takes a few minutes
Sys.time()

latent_acw2.init.qn = ergmm(movie_net_acw_full ~ euclidean(d = 2) + rsender + rreceiver,
                            response = "ratings_diff", family = "Poisson.log",
                            control = control.ergmm(pilot.runs = 4,  burnin = 0, 
                                                    interval = 500, sample.size = 10000), seed = 123,
                            user.start = with(model_acw_full$map, list(Z=Z, sender=sender, receiver=receiver,
                                                                       beta = beta, sender.var = sender.var, receiver.var = receiver.var,
                                                                       Z.var = Z.var))
)
Sys.time()
  # random initialization values ####
N = 128; n = 128; D= 2; d = 2
v_z = sqrt(N); s2_z = N/8
Z_dist = as.matrix(dist(M))
Z = cmdscale(Z_dist, k = D)
Z = Z/max(abs(Z))
Z_dist = as.matrix(dist(Z, upper = T))
sigma2_z = var(as.vector(Z))
# latentnet type initialization - worse than random for sender, receiver, beta
#a = logit( (rowSums(Y!=0) + 1)/(N-1+2) ) - (1/N) * sum(logit( (rowSums(Y!=0) + 1) / (N - 1 + 2)) )
#b = logit( (colSums(Y!=0) + 1)/(N-1+2) ) - (1/N) * sum(logit( (colSums(Y!=0) + 1) / (N - 1 + 2)) )
#B = ( 1/(N*(N-1)) * sum(Y>mean(Y)) + mean(Z_dist))
a = rnorm(N) #rep(0, N) # #sender
b = rnorm(N) #rep(0, N) # #receiver
sigma2_a = var(a)
sigma2_b = var(b)
B = 0 #intercept


# Check Results - Either Example, likelihood ####

model_qn = model_cite2 #model_acw_full
model_mcmc = latent.srp2 #latent_acw2 
#col2 <- alpha(color1[cutree(journals.cluster, h = 0.6)], 0.6) 
col2 <- as.factor(movie_net_acw_full%v%"genres")
# Check against last run
Z = model_qn$Z; a = model_qn$a; b = model_qn$b; B = model_qn$B
# Check against maximum likelihood achieved
Z = model_qn$Z.map; a = model_qn$abest; b = model_qn$bbest; B = model_qn$Bbest
sigma2_a = model_qn$sigma2_abest; sigma2_b = model_qn$sigma2_bbest; sigma2_z = model_qn$sigma2_zbest
z = t(Z)
Z_dist = as.matrix(dist(Z, upper = T))

# Compare likelihoods for quasi-Newton, mcmc.mle, and mcmc starts (in particular burnin.start found via optimization methods)
llik(Y, model_qn$abest, model_qn$bbest, model_qn$Bbest, model_qn$Z.map, model_qn$sigma2_abest, model_qn$sigma2_bbest, model_qn$sigma2_zbest) 

llik(Y, model_mcmc$mcmc.mle$sender, model_mcmc$mcmc.mle$receiver, model_mcmc$mcmc.mle$beta,
     model_mcmc$mcmc.mle$Z, model_mcmc$mcmc.mle$sender.var, model_mcmc$mcmc.mle$receiver.var,
     model_mcmc$mcmc.mle$Z.var, est = "MLE") 

llik(Y, model_mcmc$mkl$sender, model_mcmc$mkl$receiver,
     model_mcmc$mkl$beta, model_mcmc$mkl$Z, model_mcmc$mcmc.mle$sender.var,
     model_mcmc$mcmc.mle$receiver.var, model_mcmc$mcmc.mle$Z.var, est = "MLE")

llik(Y, model_mcmc$mle$sender, model_mcmc$mle$receiver,
     model_mcmc$mle$beta, model_mcmc$mle$Z,   1,1,1 , est = "MLE")

#sum(model_mcmc$mle$sender^2 + v_a*s2_a^2) / (N + 2 + v_a),
#sum(model_mcmc$mle$receiver^2 + v_b*s2_b^2) / (N + 2 + v_b),
#sum(model_mcmc$mle$Z^2 + v_z*s2_z^2) / (N*d + 2 + v_z))

llik(Y, model_mcmc$burnin.start$sender, model_mcmc$burnin.start$receiver,
     model_mcmc$burnin.start$beta, model_mcmc$burnin.start$Z,
     model_mcmc$burnin.start$sender.var, model_mcmc$burnin.start$receiver.var, model_mcmc$burnin.start$Z.var) 

llik(Y, model_mcmc$sampling.start$sender, model_mcmc$sampling.start$receiver, model_mcmc$sampling.start$beta,
     model_mcmc$sampling.start$Z,model_mcmc$sampling.start$sender.var,
     model_mcmc$sampling.start$receiver.var, model_mcmc$sampling.start$Z.var) 

#llik(Y, model_mcmc$sample$sender[10000,], model_mcmc$sample$receiver[10000,],
#     model_mcmc$sample$beta[10000], model_mcmc$sample$Z[10000,1:47,1:2],
#     model_mcmc$sample$sender.var[10000], model_mcmc$sample$receiver.var[10000], model_mcmc$sample$Z.var[10000]) 

Z_dist = as.matrix(dist(model_mcmc$mle$Z, upper = T))
lambda = exp(t(model_mcmc$mle$receiver + t(model_mcmc$mle$sender - Z_dist)) + model_mcmc$mle$beta); 
sum( Y * (t(model_mcmc$mle$beta + t(model_mcmc$mle$sender - Z_dist)) + model_mcmc$mle$beta) - lambda) + sum(diag(lambda))

  # Correlations, etc ####
    # correlation in sender and receiver terms ####
with(model_qn$map, cor(sender, model_mcmc$mkl$sender))
with(model_qn$map, cor(receiver, model_mcmc$mkl$receiver))

    # compare overall ratings - more meaningful than sender, receiver on own ####
rating1 = cbind(as.vector(model_qn$map$receiver) - as.vector(model_qn$map$sender), as.vector(model_mcmc$mkl$receiver) -as.vector(model_mcmc$mkl$sender))
cor(rating1[,1], rating1[,2])
plot(rating1, type = "n")
text(rating1, as.character(1:nrow(model_mcmc$model$Ym)), srt = -45, cex = .5) #rownames(Y)

    # compare lambdas ####
Z_dist = as.matrix(dist(model_qn$map$Z, upper = T))
with(model_qn, lambdavec = t(receiver + t(sender - Z_dist)) + beta); diag(lambdavec) = NA
lambdavec = as.vector(lambdavec)[!is.na(as.vector(lambdavec))]
lambdavec.mcmc.mle = t(model_mcmc$mkl$receiver+
                         t(model_mcmc$mkl$sender - as.matrix(dist(model_mcmc$mkl$Z, upper= T)))) +
model_mcmc$mkl$beta; diag(lambdavec.mcmc.mle) = NA
lambdavec.mcmc.mle = as.vector(lambdavec.mcmc.mle)[!is.na(as.vector(lambdavec.mcmc.mle))]
cor(exp(lambdavec), exp(lambdavec.mcmc.mle))
cor(lambdavec, lambdavec.mcmc.mle)
plot(exp(lambdavec), exp(lambdavec.mcmc.mle))

  # Plots ####
par(mfrow = c(2,2))
plot(model_qn$map$Z, col = col2, pch = 20, cex = 2, main = "grad_ascent") #col = terrain.colors(nlevels(col2))[col2]
text(model_qn$map$Z, rownames(Y), cex = .5)
#legend("bottomright", legend = levels(col2), fill  = terrain.colors(nlevels(col2)))
plot(model_mcmc$mcmc.mle$Z, col = col2, pch = 20, cex = 2, main = "mcmc.mle")
text(model_mcmc$mcmc.mle$Z, rownames(Y), cex = .5)
plot(procrustes(model_mcmc$mcmc.mle$Z, model_qn$map$Z, scale = F)$Yrot, col = col2, pch = 20, cex = 2, main = "grad_ascent procrustes to mcmc.mle") 
text(procrustes(model_mcmc$mcmc.mle$Z, model_qn$map$Z, scale = F)$Yrot, rownames(Y), cex = .5)
plot(procrustes(model_mcmc$mcmc.mle$Z, model_mcmc$burnin.start$Z, scale = F)$Yrot,
     col = col2, pch = 20, cex = 2, main = "burnin.start procrustes to mcmc.mle",
     xlab = NULL, ylab = NULL) 
text(procrustes(model_mcmc$mcmc.mle$Z, model_mcmc$burnin.start$Z, scale = F)$Yrot, rownames(Y), cex = .5)

# compare positions with procrustes test
protest(model_mcmc$mkl$Z, model_qn$map$Z)
protest(model_mcmc$mkl$Z, model_mcmc$burnin.start$Z)
protest(model_mcmc$mkl$Z, model_mcmc$sampling.start$Z)

# More models and likelihoods, comparisons for time ####
latent.srp2 = ergmm(Cnet ~ euclidean(d = 2) + rsender + rreceiver,
                    response = "citations", family = "Poisson.log", seed = 30,
                    tofit = c("mcmc", "mkl", "procrustes", "mle"),
                    control = ergmm.control(burnin = 3000000, interval = 500,
                                            sample.size = 5000, mle.maxit = 100),
                    verbose = 0)

#27 minutes with 3000000, interval = 500, sample.size = 5000
llik(Y, latent.srp2$mkl$sender, latent.srp2$mkl$receiver, latent.srp2$mkl$beta, latent.srp2$mkl$Z, latent.srp2$mcmc.mle$sender.var, latent.srp2$mcmc.mle$receiver.var, latent.srp2$mcmc.mle$Z.var, est = "MAP")
#29182
#27 minutes with 5000000, interval = 1000, sample.size = 5000, same likelihood
#8 minutes with 1000000, interval = 500, sample.size = 1000, but less good estimate
#19  minutes with 3000000, interval = 500, sample.size = 1000

#####
#####
#####
# How is directly optimization done (for initialization) in latentnet? ####
# The starting paramters are sent to stats::optim function for numerical optimization with "L-BFGS-B" https://en.wikipedia.org/wiki/Limited-memory_BFGS
# the first round finds a map, the last finds mkl (assuming "mkl" on) for the positions (3 fewer paramters than first round)
# Hoff also does, but uses "SANN" or "BFGS" in optim:
# http://www.stat.washington.edu/people/pdhoff/Code/hoff_raftery_handcock_2002_jasa/dist.r
# http://www.stat.washington.edu/people/pdhoff/Code/hoff_raftery_handcock_2002_jasa/mcmc.r


#initialization
N = 47; n = 47; D= 2; d = 2; Y = t(Cmatrix)
v_z = sqrt(N); s2_z = N/8
Z_dist = as.matrix(dist(Y))
Z = cmdscale(Z_dist, k = D)
Z = Z/max(abs(Z))
Z_dist = as.matrix(dist(Z, upper = T))
sigma2_z = var(as.vector(Z))

a = rnorm(N) #rep(0, N) # #sender
b = rnorm(N) #rep(0, N) # #receiver
sigma2_a = var(a)
sigma2_b = var(b)
B = 0 #intercept

llik2 <- function(theta, Y, d, n=nrow(Y), family = "poisson", est = "MAP") {
  a = theta[1:n]
  b = theta[(n+1):(2*n)]
  B = theta[(2*n+1)]
  Z = matrix(theta[(2*n+2):(2*n+1+d*n)], ncol = d, nrow = n)
  sigma2_a = theta[2*n+2+d*n]
  sigma2_b = theta[(2*n+3+d*n)]
  sigma2_z = theta[(2*n+4+d*n)]
  return(llik(Y=Y, sender = a, receiver = b, beta = B, Z = Z,
              sender.var = sigma2_a, receiver.var =sigma2_b,
              Z.var = sigma2_z, family = family, est = est))
}

llik2(c(a, b, B, Z, sigma2_a, sigma2_b, sigma2_z), Y, d = 2)
#other optimization methods:

# i checked that this high enough enough maxit
# standard random normal initiation better than latentnet type
# now optim method is implemented in ergmm so any of these
#  methods can be used
theta.sann = optim(c(a, b, B, Z, sigma2_a, sigma2_b, sigma2_z),
                   Y = Y, d = 2, llik2, method = "SANN",
                   lower = c(rep(-Inf, 189), 0, 0, 0),
                   control=list(maxit = 200, fnscale = -1))
theta.cg = optim(c(a, b, B, Z, sigma2_a, sigma2_b, sigma2_z),
                 Y = Y, d = 2, llik2, method = "CG",
                 lower = c(rep(-Inf, 189), 0, 0, 0),
                 control=list(maxit = 200, fnscale = -1))
theta.bfgs = optim(c(a, b, B, Z, sigma2_a, sigma2_b, sigma2_z),
                   Y = Y, d = 2, llik2, method = "BFGS",
                   #lower = c(rep(-Inf, 189), 0, 0, 0),
                   control=list(maxit = 200, fnscale = -1))

llik2(theta.sann$par, Y = Y, d = 2) # ~29180-200 #28895.01 with latentnet type sender and reciever init
llik2(theta.bfgs$par, Y = Y, d = 2) # ~29180-200 #28895.01 with latentnet type sender and reciever init
llik2(theta.cg$par, Y = Y, d = 2) # ~29180-200 #28895.01 with latentnet type sender and reciever init

par(mfrow = c(3,1))
plot(model_mcmc$mcmc.mle$Z, col = col2, pch = 16)
plot(procrustes(model_mcmc$mkl$Z, matrix(theta.sann$par[(2*n+2):(2*n+1+d*n)],ncol = d, nrow = n), scale = F)$Yrot,
     col = col2, pch = 16)
plot(procrustes(model_mcmc$mkl$Z, matrix(theta.bfgs$par[(2*n+2):(2*n+1+d*n)],ncol = d, nrow = n), scale = F)$Yrot,
     col = col2, pch = 16)

plot(a, model_mcmc$burnin.start$sender)
plot(theta$par[1:47], model_qn$abest)
plot(theta$par[1:47], model_mcmc$burnin.start$sender)
plot(theta$par[1:47], model_mcmc$mcmc.mle$sender)

plot(b, model_mcmc$burnin.start$receiver)
plot(theta$par[48:94], model_qn$bbest)
plot(theta$par[48:94], model_mcmc$burnin.start$receiver)
plot(theta$par[48:94], model_mcmc$mcmc.mle$receiver)

plot(cmdscale(as.matrix(dist(Y)), k = D), col = col2, pch = 16)
plot(matrix(theta$par[(2*n+2):(2*n+1+d*n)], ncol = d, nrow = n), col = col2, pch = 16)

plot(model_mcmc$mcmc.mle$Z, col = col2, pch = 16)
tmpz = matrix(theta.sann$par[(2*n+2):(2*n+1+d*n)], ncol = d, nrow = n)
plot(tmpz, col = col2, pch = 16)

# Hoff uses his mlpY in optim function to optimize negative log prob of the graph (for binary graph),
# but not the whole posterior likelihood. To get initially optimal Z and intercept.


#####
#####
#####
# OTHER: ####
  # experimenting with deriv wrt distance; doesn't work ####
for (i in 2:128) {
  print(-Y[i,j] - Y[j,i] + 
          exp(B + a[i] + b[j] - Z_dist[i,j]) + 
          exp(B + a[j] + b[i] - Z_dist[i,j]))
}

  # more checks. these should equal 0 for ga model: ####
#  sum(-.5*1/sigma2_a + a^2/2*sigma2_a^-2) + .5*v_a*s2_a*sigma2_a^-2 - (1 + v_a/2)/sigma2_a
#  sum(-.5*1/sigma2_b + b^2/2*sigma2_b^-2) + .5*v_b*s2_b*sigma2_b^-2 - (1 + v_b/2)/sigma2_b

