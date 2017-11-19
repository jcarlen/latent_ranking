# 0. Setup---------------------------------------------------------------------------------------------------
source(llik.R)
#----------------------------------------------------------------------------------------------------
# 1. pre-calculation (outside the function so not repeated) ----
#Y = as.sociomatrix(movie_net_watch_full, ignore.eval = F, attrname = "ratings_diff")
#Y = t(Cmatrix)
Y = M #note row variance much higher than col variance
lgamma.constant = sum(lgamma(as.vector(Y+1))) # constant for factorial term

N = n = nrow(Y)
#d = 2
#R = 1
B = rnorm(1)
#Z = rnorm(N*d)
Z_dist = dist2(Y)
Z = cmdscale(Z_dist, k = d)
Z = as.vector(Z/max(abs(Z)))
a = rnorm(N)
b = rnorm(N)
u = matrix(rnorm(N*R), ncol = R)
v = matrix(rnorm(N*R), ncol = R)
sigma2_a = 1
sigma2_b = 1
sigma2_z = 1
sigma2_u = 1
sigma2_v = 1

beta.var = 9
sender.var.df = 3
receiver.var.df = 3
Z.var.df = sqrt(N)
prior.sender.var = 1
prior.receiver.var = 1
prior.Z.var = N/8

u.var.df = 3
v.var.df = 3
prior.u.var = 1
prior.v.var = 1

# if approximate likelihood with weights ####
c = 1
Yw = Y + c # for Y-based weights, add a constant to Y if it has zeros so no 0 weights
diag(Yw) = 0 # assume no self-loops
W = Yw*(log(Yw)-1) + 5 #add so all weights > 0  
diag(W) = 0
Yw_row = rowSums(Yw)
Yw_col = colSums(Yw)
W_row = rowSums(W)
W_col = colSums(W)

#----------------------------------------------------------------------------------------------------
# est Y for just maximizing graph likelihood, MAP for MAP, MAPe for map that uses exact updates of variance 
# paramters conditioned on B,Z,a,b
# For MAPe note variance parameters aren't entered in optim or returned, but calculated in llik2 (and gradient func)
#
# 
# target for citation with shift(4819, 5061)
# llik(c(latent.srp2.init.res$mkl, sender.var = latent.srp2.init.res$mcmc.mle$sender.var, receiver.var = latent.srp2.init.res$mcmc.mle$receiver.var, Z.var = latent.srp2.init.res$mcmc.mle$Z.var), Y = Y, est= "MAP")
# target for movie_acw with shift (37332.87, 37998.51)
# target for movie_watch based on mle (-565856.5, -569542.9)
# llik(c(latent_acw2_res$mkl, sender.var = latent_acw2_res$mcmc.mle$sender.var, receiver.var = latent_acw2_res$mcmc.mle$receiver.var, Z.var = latent_acw2_res$mcmc.mle$Z.var), Y = Y, est= "MAP")
# llik(latent.watch.mle$start, Y = Y, est = "Y")
#----------------------------------------------------------------------------------------------------
# 2. Y ------------------------------------------------------------------------------------------------
#  exact ####
t1 = Sys.time()
theta.bfgs = optim(par = c(B, Z, a, b), 
                   fn = llik2, Y = Y, d = 2, est = "Y", 
                   gr = llik_gr,
                   method = "L-BFGS-B",
                   lower = c(rep(-50, 1 + N*(d+2))),
                   control=list(maxit = 2000, fnscale = -1)) #<- need high enough maxit
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
# 3. MAP ------------------------------------------------------------------------------------------------
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
# 4. MAPe (better, i think) ------------------------------------------------------------------------------
#  exact ####
t1 = Sys.time()
theta.bfgs = optim(par = c(B, Z, a, b), 
                   fn = llik2, Y = Y, d = 2, est = "MAPe", 
                   gr = llik_gr, 
                   method = "L-BFGS-B",
                   lower = c(rep(-Inf, 1 + N*(d+2))),
                   control=list(maxit = 1000, fnscale = -1))
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
# 5. Multi-dimensional sender and receiver ------------------------------------------------------------
family = "poisson_l"
d = 2
R = 1
# Y ------------------------------------------------------------------------------------------------

# depending on the model, might want to remove StataJ since it has only 4 in-edges -> instability to estimate 
# two receiver paramters and position
tmp_l  #poisson_l, no gradient
llik2(tmp_l, Y = Y, d = 2, R = 1, est= "Y", family = "poisson_l")

t1 = Sys.time()
theta.bfgs = optim(par = c(B, Z, a, b, u, v), 
                   fn = llik2, Y = Y, d = d, R = R, est = "Y", family = family,
                   gr = llik_gr,
                   method = "L-BFGS-B",
                   lower = c(rep(-50, 1 + N*(d+2*R+2))), #upper = c(rep(20, 1 + N*(d+2*R+2))),
                   control=list(maxit = 5000, fnscale = -1)) #<- need high enough maxit
(Sys.time() - t1)
(Sys.time() - t1)/theta.bfgs$counts[1]
theta.bfgs$counts[1]
theta.bfgs$convergence
llik2(theta.bfgs$par, Y = Y, d = d, R = R, est= "Y", family = family)

# --------------------------------------------------------------------------------------------
# MAP ------------------------------------------------------------------------------------------------
t1 = Sys.time()
theta.bfgs = optim(par = c(B, Z, a, b, u, v, sigma2_a, sigma2_b, sigma2_z, sigma2_u, sigma2_v), 
                   fn = llik2, Y = Y, d = 2, R = 1, est = "MAP", family = family,
                   gr = llik_gr, 
                   method = "L-BFGS-B",
                   lower = c(rep(-50, 6 + N*(d+2*R+2))), #-Inf lb on B -> error if too neg
                   control=list(maxit = 2000, fnscale = -1)) 
(Sys.time() - t1)
(Sys.time() - t1)/theta.bfgs$counts[1]
theta.bfgs$counts[1]
theta.bfgs$convergence #ERROR
llik2(theta.bfgs$par, Y = Y, d = 2, est = "Y", family = family)
llik2(theta.bfgs$par, Y = Y, d = 2, est = "MAP", family = family)

# ------------------------------------------------------------------------------------------------
# MAPe ------------------------------------------------------------------------------------------------

t1 = Sys.time()
theta.bfgs = optim(par = c(B, Z, a, b, u, v), 
                   fn = llik2, Y = Y, d = 2, R = 1, est = "MAPe", family = family,
                   gr = llik_gr, 
                   method = "L-BFGS-B",
                   lower = c(rep(-50, 6 + N*(d+2*R+2))), #-Inf lb on B -> error if too neg
                   control=list(maxit = 2000, fnscale = -1)) 
(Sys.time() - t1)
(Sys.time() - t1)/theta.bfgs$counts[1]
theta.bfgs$counts[1]
theta.bfgs$convergence #ERROR
llik2(theta.bfgs$par, Y = Y, d = 2, R = 1, est = "Y", family = family)
llik2(theta.bfgs$par, Y = Y, d = 2, R = 1, est = "MAPe", family = family)

# ------------------------------------------------------------------------------------------------

# plots ####
theta = theta.bfgs$par
B = theta[1]
Z = matrix(theta[2:(1+d*n)], ncol = d, nrow = n)
Z_dist = dist2(Z)
a = theta[(2+d*n) : (1 + n+d*n)]
b = theta[(n+2+d*n) : (1+(2+d)*n)] 
u = theta[(2+(2+d)*n) : (1+(2+d+R)*n)] 
v = theta[(2+(2+d+R)*n): (1+(2+d+2*R)*n)]
#lambda = exp(B + outer(a, b, "+") + u %*% t(v) - Z_dist); diag(lambda) = 0
#sd(Y - lambda)

col1 = as.numeric(as.factor(movie_net_acw_full%v%"genres"))
model = latent_acw2$mkl
# model = latent.srp2.init.r$mkl
# col1 = col2

par(mfrow = c(1,2)); par(mar = rep(2.5,4))
plot(Z, col = col1, pch = 16, cex = 2, xlim = c(-2,2), ylim = c(-2,2))
text(Z, rownames(Y), cex = .6, col = col1)
plot(model$Z, col = col1, pch = 16, cex = 2, xlim = c(-2,2), ylim = c(-2,2))
text(model$Z, rownames(Y), cex = .6, col = col1)


par(mfrow = c(1,3))
cor(b - a, model$receiver - model$sender)
plot(b - a, model$receiver - model$sender, pch = 16, col = col1, type = "n")
text(b - a, model$receiver - model$sender, rownames(Y), col = col1, srt = -45, cex = .5)

cor(v - u, model$receiver - model$sender)
plot(v - u, model$receiver - model$sender, type = "n")
text(v - u, model$receiver - model$sender, rownames(Y), col = col1, srt = -45, cex = .5)

cor(b + v - a - u, model$receiver - model$sender)
plot(b + v - a - u, model$receiver - model$sender, type = "n")
text(b + v - a - u, model$receiver - model$sender, rownames(Y), col = col1, srt = -45, cex = .5)

dist_inv1 = 1/(dist2(Z)+1); diag(dist_inv1) = 0
u2 = rowMeans(u*dist_inv1) 
v2 = rowMeans(v*dist_inv1) #remember dist_inv1 symmetric
cor(b + v2 - a - u2, model$receiver - model$sender) # new ranking
cor(b + v - a - u, model$receiver - model$sender)

plot(rep(0, N), u, col = col1, type = "n"); text(rep(0, N), u, rownames(Y), col = col1)
points(v, rep(0, N), col = col1, type = "n"); text(v, rep(0, N), rownames(Y), col = col1, srt = 90)

text(u, v, rownames(Y), col = col1)
plot(b + v - a - u, model$mkl$receiver - model$mkl$sender)

plot(matrix(u, ncol= R), pch = 16, col = col1)
points(matrix(v, ncol= R), pch = 10, col = col1)
circplot2(Y, U = matrix(u, ncol= R), V = matrix(v, ncol= R), pscale = .8,
          rcol = col1, ccol = col1, lty = 1, Ymin = 3, lwd = .01)
circplot2(Y, U = matrix(u, ncol= R), V = matrix(v, ncol= R), pscale = .8,
          rcol = col1, ccol = col1, lty = 1, Ymin = 3, lwd = .01, 
          visNetwork = T, group = movie_net_acw_full%v%"genres")
circplot2(Y, U = acw.ame2$U, V = acw.ame2$V, pscale = .8,
          rcol = col1, ccol = col1, lty = 1, Ymin = 3, lwd = .01)
          