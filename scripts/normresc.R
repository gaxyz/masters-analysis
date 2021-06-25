# Define scaling function
normalScaling <- function( x ){
    
    # Compute robust estimates of mean and sd
    fit <- MASS::rlm(x~1)
    mu <- fit$coefficients[1]
    s <- fit$s
    y_hat <- (x - mu) / s
    return(y_hat)
    
}