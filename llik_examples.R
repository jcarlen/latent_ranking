# 0. Setup---------------------------------------------------------------------------------------------------
source(llik.R)
library(extraDistr)
#----------------------------------------------------------------------------------------------------
family = "poisson"
d = 2
R = 1
p = 1
# 1. Initialize ----
#Y = as.sociomatrix(movie_net_watch_full, ignore.eval = F, attrname = "ratings_diff")
#Y = t(Cmatrix)
#Y = M #note row variance much higher than col variance
#Y = as.matrix(movie_rating_ra)
Y = as.matrix(movie_rating_r); row.names(Y) = movie_net_r_full%v%"titles"

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

#    if approximate likelihood, weights ####
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
(Sys.time() - t1)/theta.poisson$counts[1]
theta.poisson$counts[1]
theta.poisson$convergence
llik2(theta.poisson$par, Y = Y, d = 2, est= "Y")

#  approx - only works with gr function #### 
t1 = Sys.time()
theta.poisson = optim(par = c(B, Z, a, b), 
                   fn = llik_hat2, Y = Y, d = 2, est = "Y", r = .25, margin = "none",
                   gr = llik_gr_hat,
                   method = "L-BFGS-B",
                   lower = c(rep(-5000, 1 + N*(d+2))),
                   control=list(maxit = 100, fnscale = -1)
)
Sys.time() - t1
(Sys.time() - t1)/theta.poisson$counts[1]
theta.poisson$counts[1]
theta.poisson$convergence
llik2(theta.poisson$par, Y = Y, d = 2, est= "Y")

# ---------------------------------------------------------------------------------------------------
# 3. MAP ------------------------------------------------------------------------------------------------
#  exact ####
t1 = Sys.time()
theta.poisson = optim(par = c(B, Z, a, b, sigma2_a, sigma2_b, sigma2_z), 
                   fn = llik2, Y = Y, d = d, est = "MAP", family = "poisson",
                   gr = llik_gr, 
                   method = "L-BFGS-B",
                   lower = c(-50, rep(-Inf, N*(d+2)), 0, 0, 0), #-Inf lb on B -> error if too neg
                   control=list(maxit = 2000, fnscale = -1)) 
(Sys.time() - t1)
(Sys.time() - t1)/theta.poisson$counts[1]
theta.poisson$counts[1]
theta.poisson$convergence #ERROR
llik2(theta.poisson$par, Y = Y, d = d, est= "Y")
llik2(theta.poisson$par, Y = Y, d = d, est= "MAP")

# ---------------------------------------------------------------------------------------------------
# 4. MAPe (better, i think) ------------------------------------------------------------------------------
#  exact ####
t1 = Sys.time()
theta.poisson = optim(par = c(B, Z, a, b), 
                   fn = llik2, Y = Y, d = 2, est = "MAPe", 
                   gr = llik_gr, 
                   method = "L-BFGS-B",
                   lower = c(rep(-Inf, 1 + N*(d+2))),
                   control=list(maxit = 1000, fnscale = -1))
Sys.time() - t1
(Sys.time() - t1)/theta.poisson$counts[1]
theta.poisson$counts
theta.poisson$convergence
llik2(theta.poisson$par, Y = Y, d = 2, est= "Y")
llik2(theta.poisson$par, Y = Y, d = 2, est= "MAPe")

#  approx ####
t1 = Sys.time()
theta.poisson = optim(par = c(B, Z, a, b), 
                   fn = llik_hat2, Y = Y, d = 2, est = "MAPe", r = .25, margin = "none",
                   gr = llik_gr_hat, 
                   method = "L-BFGS-B",
                   lower = c(rep(-Inf, 1 + N*(d+2))),
                   control=list(maxit = 2000, fnscale = -1))
Sys.time() - t1
(Sys.time() - t1)/theta.poisson$counts[1]
theta.poisson$counts
theta.poisson$convergence
llik2(theta.poisson$par, Y = Y, d = 2, est= "Y")
llik2(theta.poisson$par, Y = Y, d = 2, est= "MAPe")

# ------------------------------------------------------------------------------------------------
# 5. Plots -----
theta.poisson = theta.poisson$par
B.poisson = theta.poisson[1]
Z.poisson = matrix(theta.poisson[2:(1+d*n)], ncol = d, nrow = n)
Z_dist.poisson = dist2(Z.poisson)
dist_inv.poisson = 1/Z_dist.poisson; diag(dist_inv.poisson) = 0
a.poisson = theta.poisson[(2+d*n) : (1 + n+d*n)]
b.poisson = theta.poisson[(n+2+d*n) : (1+(2+d)*n)] 
#lambda = exp(B + outer(a, b, "+") + outer(u, v, "+")*dist_inv1 - Z_dist); diag(lambda) = 0
#sd(Y - lambda)

net1 = movie_net_ra_full
col1 = as.numeric(as.factor(net1%v%"genres"))
#col1 = as.numeric(as.factor(movie_net_acw_full%v%"genres"))
#model = latent_acw2$mkl
#model = latent.srp2.init.r$mkl
#col1 = col2

  
# positions
par(mfrow = c(1,2)); par(mar = rep(2.5,4))
plot(Z.poisson, col = col1, pch = 16, cex = 1)
text(Z.poisson, net1%v%"titles", cex = .5, col = col1)
plot(model$Z, col = col1, pch = 16, cex = 2)
text(model$Z, rownames(Y), cex = .6, col = col1)

#global ratings 
cor(b.poisson - a.poisson, net1%v%"stars")
plot(b.poisson - a.poisson, net1%v%"stars", col = col1, pch = 16)
text(b.poisson - a.poisson, net1%v%"stars", net1%v%"titles", col = col1, cex = .5)
     
plot(b.poisson - a.poisson, model$receiver - model$sender,
     pch = 16, col = col1, cex = 1, main = "global")
text(b.poisson - a.poisson, model$receiver - model$sender, rownames(Y),
     col = col1, srt =-30, cex = .5, pos = 2)


# ------------------------------------------------------------------------------------------------
# 6. Multi-dimensional sender and receiver ------------------------------------------------------------
family = "poisson_m"
d = 2
R = 2
p = 1
# Initialize ----
#Y = as.sociomatrix(movie_net_watch_full, ignore.eval = F, attrname = "ratings_diff")
#Y = t(Cmatrix)
#Y = M #note row variance much higher than col variance
Y = as.matrix(movie_rating_ra); rownames(Y) = movie_net_ra_full%v%"titles"
#Y = as.matrix(movie_rating_r); rownames(Y) = movie_net_r_full%v%"titles"

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


# Y ------------------------------------------------------------------------------------------------

# depending on the model, might want to remove StataJ since it has only 4 in-edges -> instability to estimate 
# two receiver paramters and position
tmp_l  #poisson_l, no gradient
llik2(tmp_l, Y = Y, d = 2, R = 1, est= "Y", family = "poisson_l")
tmp_l2  #poisson_l, no gradient, p = 2
llik2(tmp_l2, Y = Y, d = 2, R = 1, p = 2, est= "Y", family = "poisson_l")

t1 = Sys.time()
theta.bfgs = optim(par = c(B, Z, a, b, u, v), 
                   fn = llik2, Y = Y, d = d, R = R, p = p, est = "Y", family = family,
                   gr = llik_gr,
                   method = "L-BFGS-B",
                   lower = c(rep(-50, 1 + N*(d-2+2*(R-1)+2)), rep(0, 2*N)), #upper = c(rep(20, 1 + N*(d+2*R+2))),
                   control=list(maxit = 5000, fnscale = -1)) #<- need high enough maxit
(Sys.time() - t1)
(Sys.time() - t1)/theta.bfgs$counts[1]
theta.bfgs$counts[1]
theta.bfgs$convergence
llik2(theta.bfgs$par, Y = Y, d = d, R = R, est= "Y", family = family, p = p)

# --------------------------------------------------------------------------------------------
# MAP ------------------------------------------------------------------------------------------------
t1 = Sys.time()
if (d == 0) {Z = NULL; sigma2_z = NULL}
theta.bfgs = optim(par = c(B, Z, a, b, u, v, sigma2_a, sigma2_b, sigma2_z, sigma2_u, sigma2_v), 
                   fn = llik2, Y = Y, d = d, R = R, p = p, est = "MAP", family = family,
                   gr = llik_gr, 
                   method = "L-BFGS-B",
                   lower = c(rep(-20, 1 + N*(d+2*R+2)), rep(0, 5)),
                   #lower = c(rep(-50, 1 + N*(d+2*(R-1)+2)), rep(0, 2*N), rep(-50, 5)), #-Inf lb on B -> error if too neg
                   control=list(maxit = 2000, fnscale = -1)) 
(Sys.time() - t1)
(Sys.time() - t1)/theta.bfgs$counts[1]
theta.bfgs$counts[1]
theta.bfgs$convergence #ERROR
llik2(theta.bfgs$par, Y = Y, d = d, R = R, p = p, est = "Y", family = family)
llik2(theta.bfgs$par, Y = Y, d = d, R = R, p = p, est = "MAP", family = family)

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
dist_inv = 1/Z_dist; diag(dist_inv) = 0
a = theta[(2+d*n) : (1 + n+d*n)]
b = theta[(n+2+d*n) : (1+(2+d)*n)] 
u = theta[(2+(2+d)*n) : (1+(2+d+R)*n)] 
v = theta[(2+(2+d+R)*n): (1+(2+d+2*R)*n)]
U = matrix(u, ncol = R)
V = matrix(v, ncol = R)
#lambda = exp(B + outer(a, b, "+") + outer(u, v, "+")*dist_inv1 - Z_dist); diag(lambda) = 0
#sd(Y - lambda)

col1 = as.numeric(as.factor(movie_net_ra_full%v%"genres"))
#col1 = as.numeric(as.factor(movie_net_acw_full%v%"genres"))
model = latent_acw2.coef$start
#model = latent_acw2$mkl
#model = latent.srp2.init.r$mkl
#col1 = col2
dist_inv1 = 1/(dist2(Z)+1)^p; diag(dist_inv1) = NA
u2 = rowMeans(u*dist_inv1, na.rm = T) 
v2 = rowMeans(v*dist_inv1, na.rm = T) #remember dist_inv1 symmetric

# positions
par(mfrow = c(1,2)); par(mar = rep(2.5,4))
plot(Z, col = col1, pch = 16, cex = 1)
text(Z, movie_net_ra_full%v%"titles", cex = .5, col = col1)
plot(model$Z, col = col1, pch = 16, cex = 2)
text(model$Z, rownames(Y), cex = .6, col = col1)

#ratings, without text:
par(mfrow = c(1,3))
plot(b - a, model$receiver - model$sender, pch = 16, col = col1, cex = 4)
plot(v2 - u2, model$receiver - model$sender, pch = 16, col = col1, cex = 4)
plot(b + v2 - a - u2, model$receiver - model$sender, pch = 16, col = col1, cex = 4,
     xlim = c(-3.5, 4.8), ylim = c(-2, 2.6))
text(b + v2 - a - u2, model$receiver - model$sender, rownames(Y),
     col = col1, srt = -40, cex = .8, pos = 4, offset = 1.3)

#more ratings ####
#global
cor(b - a, model$receiver - model$sender)
plot(b - a, model$receiver - model$sender, pch = 16, col = col1, cex = 1, main = "global")
     #xlim = c(-4.5, 3.5), ylim = c(-2, 3))
text(b - a, model$receiver - model$sender, rownames(Y),
              col = col1, srt =-30, cex = .5, pos = 2)
#abline(lm( (model$receiver - model$sender)[-c( 49 , 65 ,102, 61 )] ~
 #            (b-a)[-c( 49 , 65 ,102, 61 )]), col = "grey")

#local
#cor(v - u, model$receiver - model$sender)
#plot(v - u, model$receiver - model$sender, type = "n")
#text(v - u, model$receiver - model$sender, rownames(Y), col = col1, srt = -45, cex = 1)

# or rowMeans((v-u)*dist_inv1)
cor(v2 - u2, model$receiver - model$sender)
plot(v2 - u2, model$receiver - model$sender, pch = 16, col = col1, cex = 1, main = "local")
     #xlim = c(-2, 2), ylim = c(-4, 3))
text(v2 - u2, model$receiver - model$sender, rownames(Y),#paste0(round(movie_net_acw_full%v%"stars", 1)),
     col = col1, srt = -20, cex = .5, pos = 4)
abline(lm( (model$receiver - model$sender)[-c( 49 , 65 ,102, 61 )] ~
             (v2-u2)[-c( 49 , 65 ,102, 61 )]), col = "grey")
#text(v2 - u2, model$receiver - model$sender, as.character(1:nrow(Y)), col = col1, srt = -45, cex = .7, pos = 2)

# global + local
cor(b + v2 - a - u2, model$receiver - model$sender) # new ranking
#cor(b + v2 - a - u2, movie_net_acw_full%v%"stars") 
cor(model$receiver - model$sender, movie_net_acw_full%v%"stars") 
plot(b + v2 - a - u2, model$receiver - model$sender, pch = 16, col = col1, cex = 1)
     #xlim = c(-3, 3), ylim = c(-4, 3.5))
text(b + v2 - a - u2, model$receiver - model$sender, rownames(Y),
     #paste0(round(movie_net_acw_full%v%"stars", 1)),
     col = col1, srt = -45, cex = .5, pos = 2)

plot(b + v2 - a - u2, net1%v%"stars", pch = 16, col = col1, cex = 1)
#xlim = c(-3, 3), ylim = c(-4, 3.5))
text(b + v2 - a - u2, net1%v%"stars", net1%v%"titles",
     #paste0(round(movie_net_acw_full%v%"stars", 1)),
     col = col1, srt = 0, cex = .5, pos = 2)


#text(b + v2 - a - u2, model$receiver - model$sender, , col = col1, srt = -45, cex = 1, pos = 2)
plot((b + v2 - a - u2)[movie_net_ra_full%v%"genres" == "Romance"],
     (movie_net_ra_full%v%"stars")[movie_net_ra_full%v%"genres" == "Romance"],
     col = col1[movie_net_ra_full%v%"genres" == "Romance"], pch = 16)
text((b + v2 - a - u2)[movie_net_ra_full%v%"genres" == "Romance"],
     (movie_net_ra_full%v%"stars")[movie_net_ra_full%v%"genres" == "Romance"],
     rownames(Y)[movie_net_ra_full%v%"genres" == "Romance"],
     col = col1[movie_net_ra_full%v%"genres" == "Romance"], cex = .5)

tmp = inner_join(data.frame(id = movie_net_r_full%v%"titles", rating = b2 - a2),
                 data.frame(id = movie_net_ra_full%v%"titles", rating = b - a), by = id)
a2 = theta.bfgs$par[(2+d*n) : (1 + n+d*n)]
b2 = theta.bfgs$par[(n+2+d*n) : (1+(2+d)*n)] 
plot((b + v2 - a - u2)[movie_net_ra_full%v%"genres" == "Romance"|
                         movie_net_ra_full%v%"genres" == "Horror|Romance"][-9],
     (b2 - a2),
     col = col1[movie_net_ra_full%v%"genres" == "Romance"], pch = 16)

plot(b - a, b + v2 - a - u2, col = col1[movie_net_ra_full%v%"genres" == "Romance"], pch = 16)
          
plot(model$sender[order(model$sender)], pch = 16, col = col1[order(model$sender)])
plot(model$receiver[order(model$receiver)], pch = 16, col = col1[order(model$receiver)])
cor(u2, model$sender)
plot(u2, model$sender, pch = 16, col = col1)
#text(u2, model$sender, as.character(1:nrow(Y)), col = col1, srt = -45, cex = .7, pos = 2)
text(u2, model$sender, rownames(Y), col = col1, srt = -45, cex = .7, pos = 2)
plot(v2, model$receiver, pch = 16, col = col1)
#text(v2, model$receiver, as.character(1:nrow(Y)), col = col1, srt = -45, cex = .7, pos = 2)
text(v2, model$receiver, rownames(Y), col = col1, srt = -45, cex = .7, pos = 2)
plot(v2/(abs(b) + abs(v2)), pch = 16, col = col1)

# multidimesional, for poisson_m ####
plot(matrix(u, ncol= R), pch = 16, col = col1)
points(matrix(v, ncol= R), pch = 10, col = col1)
circplot2(Y, U = matrix(u, ncol= R), V = matrix(v, ncol= R), pscale = .8,
          rcol = col1, ccol = col1, lty = 1, Ymin = 3, lwd = .01)
circplot2(Y, U = matrix(u, ncol= R), V = matrix(v, ncol= R), pscale = .8,
          rcol = col1, ccol = col1, lty = 1, Ymin = 3, lwd = .01, 
          visNetwork = T, group = movie_net_acw_full%v%"genres")
circplot2(Y, U = acw.ame2$U, V = acw.ame2$V, pscale = .8,
          rcol = col1, ccol = col1, lty = 1, Ymin = 3, lwd = .01)

# visNetwork -----
net1 = Cnet#movie_net_ra_full
vc = r1

nodes <- data.frame(id = 1:length(net1%v%"vertex.names"), 
                    label = "",
                    #label = movie_net_acw%v%"vertex.names",
                    title = paste(net1%v%"vertex.names"),#, round(net1%v%"stars",2)), #for selection dropdown
                    value = vc, 
                    group = net1%v%"group")
                    #group = net1%v%"genres") #conveys node size, on a scale from 0-1

edges <- data.frame(from=data.frame(as.edgelist(net1))$X1, 
                    to=data.frame(as.edgelist(net1))$X2)
#value = movie_net_acw%e%"avg_diff")

net_igraph <- igraph::graph.data.frame(edges, directed=TRUE, vertices=nodes)

visNetwork(nodes, edges, main = "Network", submain = "latent space positions") %>%
  visIgraphLayout(acw_igraph, layout = "layout.norm",
                  layoutMatrix = Z, type = "full") %>%
  #visOptions(highlightNearest = list(enabled =TRUE, degree = 1)) %>%
  visNodes(#color = list(highlight = list(background = "black")),
    shape = "dot", #shape = "text"
    scaling = list(min = 5, max = 30, label = list(max=20)),
    labelHighlightBold = TRUE,
    borderWidth = .5,
    borderWidthSelected = 2,
    font= list(size = 20, color = "black", strokeWidth = 1)) %>%
  visEdges(shadow = FALSE,
           scaling = list(min = 1, max = 10),
           #selectionWidth = "function(x) {acw_avg_diff}",
           hidden = FALSE,
           #scaling = list(min = 1, max = 10),
           arrows = list(to = list(enabled = FALSE, scaleFactor = 5),
                         middle = list(enabled = TRUE, scaleFactor = 1)),
           color = list(inherit = TRUE, highlight = "red", opacity = 0.04), #, color = "white"
           smooth = list(type = "curvedCW")) %>%
  visOptions(highlightNearest = list(enabled = T, degree = 0.5),
             nodesIdSelection = FALSE, 
             selectedBy = "title") %>%
  visLegend(enabled = TRUE, useGroups=TRUE) %>%
  visInteraction(selectConnectedEdges = TRUE)
