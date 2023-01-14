
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

overdispersion <- function(object){
    
    ovds <- sum(residuals(object, type = "pearson")^2)/(object@df.residual)
    
    ovds
}


vcov_disp <- function(object, ovds){
    
    vcov.tilde <- ovds * vcov(object)
    
    vcov.tilde
}
pat_vcov <- function(n_cl = n_class){
    
    if(n_cl == 3){pattern <- c(1, 3, 2, 4)}
    
    if(n_cl == 4){pattern <- c(1, 4, 2, 5, 3, 6)}
    
    if(n_cl == 5){pattern <- c(1, 5, 2, 6, 3, 7, 4, 8)} 
    
    return(pattern)
}

summary_results <- function(pfa_obj, tmp_seed, cut_offs, p1){
    
    count_p1 <- c(); count_adj_p1 <- c(); p_val_p1 <- c(); adj_p1 <- c()
    
    adj_total <- c(); adj_p <- c()
    
    
    for(i in 1:length(cut_offs)){
        
        p_val_p1 <-  dim(pfa_obj$Pvalue[pfa_obj$Pvalue$p.value <= cut_offs[i] &
                                            pfa_obj$Pvalue$Index %in% 1:p1,])[1]
        
        adj_p1   <-  dim(pfa_obj$adjPvalue[pfa_obj$adjPvalue$p.value <= cut_offs[i] &
                                               pfa_obj$adjPvalue$Index %in% 1:p1,])[1]
        
        adj_p <-  dim(pfa_obj$adjPvalue[pfa_obj$adjPvalue$p.value <= cut_offs[i],])[1]
        
        
        count_p1 <- rbind(count_p1, p_val_p1)
        
        count_adj_p1 <- rbind(count_adj_p1, adj_p1)
        
        adj_total <- rbind(adj_total, adj_p)
    }
    res = data.frame('seed'    = tmp_seed, 
                     't'        = cut_offs,
                     'R'        = pfa_obj$FDP$rejects, 
                     'V'        = pfa_obj$FDP$false.rejects,
                     'po'       = pfa_obj$pi0,
                     'FDP'      = pfa_obj$FDP$FDP, 
                     'K'        = pfa_obj$K, 
                     'Gen_P1'   = count_p1,
                     'Adj_P'    = adj_total,
                     'Adj_P1'   = count_adj_p1)
    
    res
}

