mmm = function(X., y.){
    
    n <-  dim(X.)[1]; p <-  dim(X.)[2]
    
    if(is.vector(y.)){ n_class <-  length(unique(y.))}
    
    if(is.matrix(y.)){ n_class <-  dim(y.)[2]} 
    
    Estimates <- matrix(NA, nrow = n_class - 1, ncol = p) 
    
    Psi = c()
    
    for(j in 1:p){
       
        mult.model <- VGAM::vglm(y. ~ X.[, j], multinomial())
        
        tmp_est <- coef(mult.model, matrix = TRUE)
        
        for(v in 1:dim(tmp_est)[2]){
            
            Estimates[v, j] = tmp_est[2, v]
        }
        
        B <-  getScore(object  = mult.model, Z = cbind(1, X.[,j])) 
        
          
         V1 <- vcov(mult.model)
             
        pattern <- pat_vcov(n_cl = n_class)
        
        V = V1[, c(pattern)]
        
        V = V[c(pattern),]
        
        MMM  = (V %*% t(B)) * sqrt(n)
        
        row_odd <- seq_len(nrow(MMM)) %% 2  
        
        tmp_MMM <- MMM[row_odd == 0,]
        
        Psi <- rbind(Psi, tmp_MMM)
        
    }
    return(list(Betas = Estimates,
                Psi = Psi,
                obs = n, 
                num_p = p))
}

################################################################################

marg_cov = function(psi, n, p){

    Sigma <- matrix(rep(0, p * p), nrow = p)

    tmp_psi = vector()
    
    for(i in 1:p){
        for(j in 1:p){
            
            if(i == j ){
                Sigma[i, i] = (t(psi[i,]) %*% psi[i,])/(n)
            }
            if( i < j ){
                
                Sigma[i, j] = Sigma[j, i] =  (t(psi[i,]) %*% psi[j,])/(n)
                
            }
        }
        
    }
    
    return(Sigma)
}

multi_pfa <- function(X, y, tval = list(), reg = "L1", K = list(), m = 0.95, dispersion = FALSE){
    
    if(is.vector(y) && length(unique(y)) < 3){
        stop("The outcome variable should be multicategorical!")}
    
    if(is.vector(y) && length(unique(y)) > 6){
        stop("The recent version of the package is limitied up to six classes.")}
    
    if(is.matrix(y) && dim(y)[2] < 3){
        stop("The outcome variable should be multicategorical!")}
    
    if(is.matrix(y) && dim(y)[2] > 6){
        stop("The recent version of the package is limitied up to six classes.")}
    
    marginal_models <- mmm(X. = X, y. = y)
    
    if((is.vector(y) && length(unique(y)) == 3) ||
        (is.matrix(y) && dim(y)[2] == 3)){
        
        if(!is.empty.list(tval)){  t1 = tval[[1]]; t2 = tval[[2]]} else {
            t1 <- rlang::missing_arg(); t2 <- rlang::missing_arg()}
      
        if(!is.empty.list(K)){  k1 = K[[1]]; k2 = K[[2]]} else {
            k1 <- rlang::missing_arg(); k2 <- rlang::missing_arg()}
        
        row_odd <- seq_len(nrow(marginal_models$Psi)) %% 2  
        
        s1 <- marg_cov(psi = marginal_models$Psi[row_odd == 1,],
                       n = marginal_models$obs, p = marginal_models$num_p)
        
        re_1 <- pfa(Z = marginal_models$Betas[1,], 
                    Sigma = s1, t = t1, Kmax = marginal_models$obs, 
                    reg = reg, K = k1)
        
        re_1$Pvalue <- data.frame(marginal.estimates = marginal_models$Betas[1,][re_1$Pvalue$Index], 
                                     p.value = re_1$Pvalue$p.value, 
                                       Index = re_1$Pvalue$Index)
        
        colnames(re_1$Pvalue)[1] <- "marginal.estimates(log(mu[,1]/mu[,3])"
        
        re_1$Sigma <- s1
        
        ########################################################################
        
        s2 <- marg_cov(psi = marginal_models$Psi[row_odd == 0,],
                       n = marginal_models$obs, p = marginal_models$num_p)
        
        re_2 <- pfa(Z = marginal_models$Betas[2,], 
                    Sigma = s2, t = t2, Kmax = marginal_models$obs, 
                    reg = reg, K = k2)
        
        re_2$Pvalue <- data.frame(marginal.estimates = marginal_models$Betas[2,][re_2$Pvalue$Index], 
                                  p.value = re_2$Pvalue$p.value, Index = re_2$Pvalue$Index)
       
        colnames(re_2$Pvalue)[1] <- "marginal.estimates(log(mu[,2]/mu[,3])"
        
        re_2$Sigma <- s2
        
        
        return(list(fdp_class.1_class.3 = re_1, fdp_class.2_class.3 = re_2)) 
    }
    
    if((is.vector(y) && length(unique(y)) == 4) ||
       (is.matrix(y) && dim(y)[2] == 4)){
        
        if(!is.empty.list(tval)){  t1 = tval[[1]]; t2 = tval[[2]]} else {
            t1 <- rlang::missing_arg(); t2 <- rlang::missing_arg()
            t3 <- rlang::missing_arg()}
        
        if(!is.empty.list(K)){  k1 = K[[1]]; k2 = K[[2]]} else {
            k1 <- rlang::missing_arg(); k2 <- rlang::missing_arg()
            k3 <- rlang::missing_arg()}
        
        
        s1 <- marg_cov(psi = marginal_models$Psi[seq(from = 1, 
                                                      to  = dim(marginal_models$Psi)[1],
                                                      by  = 3),],
                       n = marginal_models$obs, p = marginal_models$num_p)
        
        re_1 <- pfa(Z = marginal_models$Betas[1,], 
                    Sigma = s1, t = t1, Kmax = marginal_models$obs, 
                    reg = reg,  K = k1)
        
        re_1$Pvalue <- data.frame(marginal.estimates = marginal_models$Betas[1,][re_1$Pvalue$Index], 
                                  p.value = re_1$Pvalue$p.value, Index = re_1$Pvalue$Index)
      
        colnames(re_1$Pvalue)[1] <- "marginal.estimates(log(mu[,1]/mu[,4])"
        
        re_1$Sigma <- s1
        
        ########################################################################
        
        s2 <- marg_cov(psi = marginal_models$Psi[seq(from = 2, 
                                                     to = dim(marginal_models$Psi)[1],
                                                     by = 3),],
                       n = marginal_models$obs, p = marginal_models$num_p)
        
        
        re_2 <- pfa(Z = marginal_models$Betas[2,], 
                    Sigma = s2, t = t2, Kmax = marginal_models$obs, 
                    reg = reg,  K = k2)
        
        re_2$Pvalue <- data.frame(marginal.estimates = marginal_models$Betas[2,][re_2$Pvalue$Index], 
                                  p.value = re_2$Pvalue$p.value, Index = re_2$Pvalue$Index)
        
        colnames(re_2$Pvalue)[1] <- "marginal.estimates(log(mu[,2]/mu[,4])"
        
        re_2$Sigma <- s2
        
        ########################################################################
        
        s3 <- marg_cov(psi = marginal_models$Psi[seq(from = 3, 
                                                     to =dim(marginal_models$Psi)[1],
                                                     by = 3),],
                       n = marginal_models$obs, p = marginal_models$num_p)
        
        re_3 <- pfa(Z = marginal_models$Betas[3,], 
                    Sigma = s3, t = t3, Kmax = marginal_models$obs, 
                    reg = reg,  K = k3)
        
        re_3$Pvalue <- data.frame(marginal.estimates = marginal_models$Betas[3,][re_3$Pvalue$Index], 
                                  p.value = re_3$Pvalue$p.value, Index = re_3$Pvalue$Index)
     
        colnames(re_3$Pvalue)[1] <- "marginal.estimates(log(mu[,3]/mu[,4])"
        
        re_3$Sigma <- s3
        
        
        return(list(fdp_class.1_class.4 = re_1, 
                    fdp_class.2_class.4 = re_2, 
                    fdp_class.3_class.4 = re_3)) 
        
    }
    
    if((is.vector(y) && length(unique(y)) == 5) ||
       (is.matrix(y) && dim(y)[2] == 5)){
        
        if(is.list(t)){  t1 = t[[1]]; t2 = t[[2]]; t3 = t[[3]]; t4 = t[[3]]}
        
        if(is.list(K)){  k1 = K[[1]]; k2 = K[[2]]; k3 = K[[3]]; k4 = K[[3]]}
        
        s1 <- marg_cov(psi = marginal_models$Psi[seq(from = 1, 
                                                     to  = dim(marginal_models$Psi)[1],
                                                     by  = 4),],
                       n = marginal_models$obs, p = marginal_models$num_p)
        
        re_1 <- pfa(Z = marginal_models$Betas[1,], 
                    Sigma = s1, t = t1, Kmax = marginal_models$obs, 
                    reg = reg,  K = k1)
        
        re_1$Pvalue <- data.frame(marginal.estimates = marginal_models$Betas[1,][re_1$Pvalue$Index], 
                                  p.value = re_1$Pvalue$p.value, Index = re_1$Pvalue$Index)
        
        colnames(re_1$Pvalue)[1] <- "marginal.estimates(log(mu[,1]/mu[,5])"
        
        re_1$Sigma <- s1
        
        ########################################################################
        
        s2 <- marg_cov(psi = marginal_models$Psi[seq(from = 2, 
                                                     to = dim(marginal_models$Psi)[1],
                                                     by = 4),],
                       n = marginal_models$obs, p = marginal_models$num_p)
        
        
        re_2 <- pfa(Z = marginal_models$Betas[2,], 
                    Sigma = s2, t = t2, Kmax = marginal_models$num_p, 
                    reg = reg,  K = k2)
        
        re_2$Pvalue <- data.frame(marginal.estimates = marginal_models$Betas[2,][re_2$Pvalue$Index], 
                                  p.value = re_2$Pvalue$p.value, Index = re_2$Pvalue$Index)
        
        colnames(re_2$Pvalue)[1] <- "marginal.estimates(log(mu[,2]/mu[,5])"
        
        re_2$Sigma <- s2
        
        ########################################################################
        
        s3 <- marg_cov(psi = marginal_models$Psi[seq(from = 3, 
                                                     to =dim(marginal_models$Psi)[1],
                                                     by = 4),],
                       n = marginal_models$obs, p = marginal_models$num_p)
        
        re_3 <- pfa(Z = marginal_models$Betas[3,], 
                    Sigma = s3, t = t3, Kmax = marginal_models$obs, 
                    reg = reg,  K = k3)
        
        re_3$Pvalue <- data.frame(marginal.estimates = marginal_models$Betas[3,][re_3$Pvalue$Index], 
                                  p.value = re_3$Pvalue$p.value, Index = re_3$Pvalue$Index)
        
        colnames(re_3$Pvalue)[1] <- "marginal.estimates(log(mu[,3]/mu[,5])"
        
        re_3$Sigma <- s3
        
        ########################################################################
        
        s4 <- marg_cov(psi = marginal_models$Psi[seq(from = 4, 
                                                     to =dim(marginal_models$Psi)[1],
                                                     by = 4),],
                       n = marginal_models$obs, p = marginal_models$num_p)
        
        re_4 <- pfa(Z = marginal_models$Betas[4,], 
                    Sigma = s4, t = t4, Kmax = marginal_models$obs, 
                    reg = reg,  K = k4)
        
        re_4$Pvalue <- data.frame(marginal.estimates = marginal_models$Betas[4,][re_4$Pvalue$Index], 
                                  p.value = re_4$Pvalue$p.value, Index = re_4$Pvalue$Index)
        
        colnames(re_4$Pvalue)[1] <- "marginal.estimates(log(mu[,4]/mu[,5])"
        
        re_4$Sigma <- s4
        
        return(list(fdp_class.1_class.5 = re_1, 
                    fdp_class.2_class.5 = re_2, 
                    fdp_class.3_class.5 = re_3,
                    fdp_class.4_class.5 = re_4)) 
        
    }
    
}




