#############
# This function is taken from pfa:pfa.gwas()
#############

pfa = function (Z, Sigma, t, Kmax, reg = "L1", e = 0.05, gamma, K, 
                plot = "-log",  m = 0.95) {
				
  ## Vladimir Vutov added the additional parameter m.
  Z <- as.vector(Z)
  Sigma <- as.matrix(Sigma)
  p <- length(Z)
  SD <- sqrt(diag(Sigma))
  Z <- Z/SD
  Sigma <- diag(1/SD) %*% Sigma %*% diag(1/SD)
  pca <- svd(Sigma, nu = 0, nv = Kmax)
  lambda <- pca$d
  eigvec <- pca$v
  
  if (missing(K)) {
    K = 1
    while (K < Kmax & sqrt(sum(lambda[(K + 1):length(lambda)]^2)) >= 
           e * sum(lambda)) K = K + 1
  }

  sqrt_lambda <- as.matrix(sqrt(lambda[1:K]))
  b <- as.matrix(eigvec[, 1:K])
  for (i in 1:K) {
    b[, i] <- b[, i] * sqrt_lambda[i]
  }
  if (reg == "L1") {
    W.hat <- rq(Z ~ b - 1, 0.5)$coef
  }
  else if (reg == "L2") {
    o = order(abs(Z))
    Zperm = Z[o]
    Lperm = as.matrix(b[o, ])
    Z.reduce = Zperm[1:(p * m)]
    L.reduce = as.matrix(Lperm[1:(p * m), ])
    W.hat = lsfit(x = L.reduce, y = Z.reduce, intercept = F)$coef
  }
  rs <- rowSums(b^2)
  inv_a <- sqrt(((1 - rs) + abs(1 - rs))/2)
  bW.est <- b %*% (W.hat)
  P <- 2 * (1 - pnorm(abs(Z)))
  sort <- sort(P, index.return = TRUE)
  index <- sort$ix
  P <- sort$x
  
  if (missing(gamma)) 
    gamma <- as.numeric(quantile(P, probs = 0.4))
  p0.est <- min(p, sum(P > gamma)/(1 - gamma))
  t.default <- TRUE
  
  if (!missing(t)) {
    # if (t == "pval") { # Vladimir Vutov comments this section, on 05/12/2021
    # 
    #   t <- P
    #  t.default <- FALSE
    # }
    
    if (is.numeric(t)) 
      t.default = (sum(t >= 0) + sum(t <= 1) < 2 * length(t))
  }
  
  if (t.default) {
    logt.l <- max(min(log(P)), log(0.00000000000001))
    logt.u <- max(log(P))
    grid <- (logt.u - logt.l) * seq(from = 0.01, to = 1, 
                                    by = 0.025) * 0.5 + 0.85 * logt.l + 0.15 * logt.u
    t <- exp(grid)
  }
  FDPt <- Vt <- Rt <- rep(0, length(t))
  for (l in 1:length(t)) {
    P1 <- 2 * (1 - pnorm(abs(Z)))
    Rt[l] <- sum(P1 <= t[l])
    a <- rep(0, p)
    for (j in 1:p) {
      qtl <- qnorm(t[l]/2)
      if (inv_a[j] > 0) {
        a[j] <- pnorm((qtl + bW.est[j])/inv_a[j]) + 
          pnorm((qtl - bW.est[j])/inv_a[j])
      }
      else {
        a[j] <- as.numeric(abs(bW.est[j]) > abs(qtl))
      }
    }
    Vt[l] <- min(sum(a), Rt[l])
    if (Rt[l] == 0) {
      FDPt[l] <- 0
    }
    else {
      FDPt[l] <- Vt[l]/Rt[l]
    }
  }
  adj.P <- as.vector(rep(0, p))
  for (j in 1:p) {
    if (inv_a[j] > 0) {
      adj.P[j] <- 2 * (1 - pnorm(abs(Z[j] - bW.est[j])/inv_a[j]))
    }
    else {
      adj.P[j] <- as.numeric(abs(Z[j] - bW.est[j]) == 0)
    }
  }
  sort <- sort(adj.P, index.return = TRUE)
  adj.index <- sort$ix
  adj.P <- sort$x
  Pvals <- data.frame(p.value = P, Index = index)
  adjPvals <- data.frame(p.value = adj.P, Index = adj.index)
  if (t.default) {
    FDPvals <- data.frame(minus.logt = -log(t), rejects = Rt, 
                          false.rejects = Vt, FDP = FDPt)
  }
  else {
    FDPvals <- data.frame(t = t, rejects = Rt, false.rejects = Vt, 
                          FDP = FDPt)
  }
  results <- list(Pvalue = Pvals, adjPvalue = adjPvals, FDP = FDPvals, 
                  pi0 = p0.est/p, K = K, Lamba = lambda, sigma = NULL, Zval = Z, Sigma = Sigma)
  
  class(results) <- "FDPresults"
   if (plot == "-log") {
     par(mfrow = c(2, 2))
     hist(P, main = "Histogram of p-values", xlab = "p-values")
     plot(-log(t), Rt, xlab = "-log(t)", ylab = "", main = "Number of total rejections", 
          type = "o")
     plot(-log(t), Vt, xlab = "-log(t)", ylab = "", main = "Number of estimated false rejections", 
          type = "o")
     plot(-log(t), FDPt, xlab = "-log(t)", ylab = "", main = "Estimated FDP", 
          type = "o")
   }
   else if (plot == "linear") {
     par(mfrow = c(2, 2))
     hist(P, main = "Histogram of p-values", xlab = "p-values")
     plot(t, Rt, xlab = "t", ylab = "", main = "Number of total rejections", 
          type = "o")
     plot(t, Vt, xlab = "t", ylab = "", main = "Number of estimated false rejections", 
          type = "o")
     plot(t, FDPt, xlab = "t", ylab = "", main = "Estimated FDP", 
          type = "o")
   }
   else if (plot == "log") {
     par(mfrow = c(2, 2))
     hist(P, main = "Histogram of p-values", xlab = "p-values")
     plot(log(t), Rt, xlab = "log(t)", ylab = "", main = "Number of total rejections", 
          type = "o")
     plot(log(t), Vt, xlab = "log(t)", ylab = "", main = "Number of estimated false rejections", 
          type = "o")
     plot(log(t), FDPt, xlab = "log(t)", ylab = "", main = "Estimated FDP", 
          type = "o")
   }
  
  return(results)
}
