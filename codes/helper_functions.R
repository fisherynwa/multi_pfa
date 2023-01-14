
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
