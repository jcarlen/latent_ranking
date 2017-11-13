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


else {
  
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