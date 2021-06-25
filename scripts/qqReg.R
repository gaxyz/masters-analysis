# function for computing scaled values
qReg <- function(x, probs = seq(0,1-0.01,0.01),df){
    # Compute empirical quantiles
    eq <- quantile(x, probs=probs)
    # Compute theoretical quantiles
    tq  <- qchisq(probs, df=df ) 
    # Fit a robust linear model
    fit  <- MASS::rlm(tq~eq, maxit=100)
    # Get coefficients
    b <- fit$coefficients[1]
    a <- fit$coefficients[2]
    # Return scaled values
    y_hat <- x*a + b
    return(y_hat)
}