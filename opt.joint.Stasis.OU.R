opt.joint.Stasis.OU<-function (y, gg, cl = list(fnscale = -1), pool = TRUE, meth = "L-BFGS-B", hess = FALSE) 
{
  if (pool) {
    y <- pool.var(y, ret.paleoTS = TRUE)
  }
  small <- 1e-08
    p0st <- mle.Stasis(sub.paleoTS(y, ok = gg == 1))
    names(p0st) <- c("theta_st", "omega")
    p0OU<- mle.GRW(sub.paleoTS(y, ok = gg == 2))

    K <- 6
  
  #Checking if te variance parameter in the RW and Stasis is below critical value, and if TRUE, multiply it with 100
  if (p0st["omega"] <= small) 
    p0st["omega"] <- 100 * small
  if (p0OU["vstep"] <= small) 
    p0OU["vstep"] <- 100 * small
  
    #prepare initial guesses of OU parameters 
    halft <- (sub.paleoTS(y, ok = gg == 2)$tt[length(sub.paleoTS(y, ok = gg == 2)$tt)] - sub.paleoTS(y, ok = gg == 2)$tt[1])/4
    p0 <- c(p0OU["vstep"]/10, sub.paleoTS(y, ok = gg == 2)$mm[length(sub.paleoTS(y, ok = gg == 2)$mm)], log(2)/halft)
    names(p0) <- c("vstep", "theta_OU", "alpha")
    
  #make vector out of estimated model parameters
  p0 <- c(p0st, p0)
  
  #define lower bounds for the L-BFGS-B method
  # three parameters in the unbiased ramdom walk and 4 parameters for the biased random walk
  ll<- c(NA, small, small, NA, small)
  
  # If user wants to use a differen method than L-BFGS-B, then the lower boun is defined as -Inf
  if (meth != "L-BFGS-B") 
    ll <- -Inf
  
    #Here the multivariate parameter estimation routine start:
    w <- try(optim(p0, fn = logL.joint.Stasis.OU, gg = gg, 
                   method = meth, lower = ll, control = cl, hessian = hess, 
                   y = y), silent = TRUE)
    if (class(w) == "try-error") {
      cl <- list(fnscale = -1, parscale = c(1, 10, 100))
      w <- try(optim(p0, fn = logL.joint.Stasis.OU, gg = gg, 
                     method = meth, lower = ll, control = cl, hessian = hess, 
                     y = y), silent = TRUE)
    }
  
  if (class(w) == "try-error") {
    wc <- as.paleoTSfit(logL = NA, parameters = NA, modelName = paste("Stasis", 
                                                                      "OU", sep = "-"), method = "Joint", K = K, n = length(y$mm), 
                        se = NULL)
    return(wc)
  }
  if (hess) 
    w$se <- sqrt(diag(-1 * solve(w$hessian)))
  else w$se <- NULL
  wc <- as.paleoTSfit(logL = w$value, parameters = w$par, modelName = paste("Stasis", 
                                                                            "OU", sep = "-"), method = "Joint", K = K, n = length(y$mm), 
                      se = w$se)
  return(wc)
}