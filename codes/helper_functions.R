    
    getScore = function(object, Z){
        
        responseRes <- residuals(object, type = "response")
        
        numfac <- dim(responseRes)[2]
        
        mu_res <- responseRes[, 1:(numfac - 1)] 
        
        scoreVectors <- c()
        
        for(col in 1:(numfac - 1) ){
            
            tmpRest <- mu_res[, col] * Z
            
            scoreVectors <- cbind(scoreVectors, tmpRest)
            
        }
        
        scoreVectors = scoreVectors 
        
        return(scoreVectors)
    }
    
  
    pat_vcov <- function(n_cl = n_class){
        
        if(n_cl == 3){pattern <- c(1, 3, 2, 4)}
        
        if(n_cl == 4){pattern <- c(1, 4, 2, 5, 3, 6)}
        
        if(n_cl == 5){pattern <- c(1, 5, 2, 6, 3, 7, 4, 8)} 
        
        return(pattern)
    }
    
    ################################
    # Efron re-scaling correction
    ########################
    Ahat <- function(x0, z){
        N  <- length(z)
        Y0 <- sum((z < x0) & (z > -x0));
        P0 <- 2 * pnorm(x0) - 1;
        Q0 <- sqrt(2) * x0 * dnorm(x0);
        P0hat <- Y0 / N;
        return((P0-P0hat) / Q0);
    }

    ####################################
    getMargQuants <- function(y, X){
        
        p_numb <- dim(X)[2]
        
        cl_1.pval <- c(); cl_2.pval <- c()   
        
        cl_1.zval <- c(); cl_2.zval <- c()   
        
        for(j in 1:p_numb){
        
            g <- vglm(y ~ X[,j], family = multinomial())
           
            tmp_sum <- summary(g)
            
            cl_1.pval <- rbind(cl_1.pval, tmp_sum@coef3[3, 4])
            
            cl_2.pval <- rbind(cl_2.pval, tmp_sum@coef3[4, 4])
            
            cl_1.zval <- rbind(cl_1.zval, tmp_sum@coef3[3, 3])
            
            cl_2.zval <- rbind(cl_2.zval, tmp_sum@coef3[4, 3])
            
        }
        return(list(pval1 = cl_1.pval, pval2 = cl_2.pval,
                    zval1 = cl_1.zval, zval2 = cl_2.zval))
    }
    ###################################
    
    fdr.control <- function(p_values){
        
        
        BH_proc <- p.adjust(p_values, method = "BH")
        BY_proc <- p.adjust(p_values, method = "BY")
        
        rejects_BH  <- c(length(BH_proc[BH_proc <= 0.05]),
                         length(BH_proc[BH_proc <= 0.1]), 
                         length(BH_proc[BH_proc <= 0.15]))
        
        rejects_BY <- c(length(BY_proc[BY_proc <= 0.05]),
                        length(BY_proc[BY_proc <= 0.1]), 
                        length(BY_proc[BY_proc <= 0.15]))
        
        return(list(rejects_BH = rejects_BH, rejects_BY = rejects_BY))
    }
    
    fdr.control.efron <- function(z_values){
        
        efron <-  1 + round(Ahat(x0 = 1, z = z_values), 2)
        
        zval.efron = z_values / efron 
        
        p.val_efron = 2 * (1 - pnorm(abs(zval.efron)))
        
        BH_efron <- p.adjust(p.val_efron, method = "BH")
        
        BY_efron <- p.adjust(p.val_efron, method = "BY")
        
        rejects_BH_Efron <- c(length(BH_efron[BH_efron <= 0.05]),
                              length(BH_efron[BH_efron <= 0.1]), 
                              length(BH_efron[BH_efron <= 0.15]))
        
        rejects_BY_Efron <- c(length(BY_efron[BY_efron <= 0.05]),
                              length(BY_efron[BY_efron <= 0.1]), 
                              length(BY_efron[BY_efron <= 0.15]))
        
        
        return(list(BH_Efron = rejects_BH_Efron, BY_Efron = rejects_BY_Efron))
        
    }