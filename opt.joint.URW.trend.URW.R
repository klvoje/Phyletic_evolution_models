opt.joint.URW.trend.URW<-function (y, gg, cl = list(fnscale = -1), pool = TRUE, meth = "L-BFGS-B", hess = FALSE) 
{
  if (pool) {
    y <- pool.var(y, ret.paleoTS = TRUE)
  }
  small <- 1e-08
    p0rw_1 <- mle.URW(sub.paleoTS(y, ok = gg == 1))
    p0dt<- mle.GRW(sub.paleoTS(y, ok = gg == 2))
    p0rw_2<- mle.URW(sub.paleoTS(y, ok = gg == 3))
    
    K <- 7
  
  #Checking if te variance parameter in the RW and Stasis is below critical value, and if TRUE, multiply it with 100
  if (p0rw_1["vstep"] <= small) 
    p0rw_1["vstep"] <- 100 * small
  if (p0dt["vstep"] <= small) 
    p0dt["vstep"] <- 100 * small
  if (p0rw_2["vstep"] <= small) 
    p0rw_2["vstep"] <- 100 * small
    
  
    # Define ancestral state
    p0anc <- y$mm[1] 
    names(p0anc) <- "anc" 
   
  
    #make vector out of estimated model parameters
  p0 <- c(p0anc, p0rw_1, p0dt, p0rw_2)
  
  #define lower bounds for the L-BFGS-B method
  # three parameters in the unbiased ramdom walk and 4 parameters for the biased random walk
  ll<- c(NA, small, NA, small, small)
  
  # If user wants to use a differen method than L-BFGS-B, then the lower boun is defined as -Inf
  if (meth != "L-BFGS-B") 
    ll <- -Inf
  
    #Here the multivariate parameter estimation routine start:
    w <- try(optim(p0, fn = logL.joint.URW.trend.URW, gg = gg, 
                   method = meth, lower = ll, control = cl, hessian = hess, 
                   y = y), silent = TRUE)
    if (class(w) == "try-error") {
      cl <- list(fnscale = -1, parscale = c(1, 10, 100))
      w <- try(optim(p0, fn = logL.joint.URW.trend.URW, gg = gg, 
                     method = meth, lower = ll, control = cl, hessian = hess, 
                     y = y), silent = TRUE)
    }
  
  if (class(w) == "try-error") {
    wc <- as.paleoTSfit(logL = NA, parameters = NA, modelName = paste("URW", 
                                                                      "Trend","URW", sep = "-"), method = "Joint", K = K, n = length(y$mm), 
                        se = NULL)
    return(wc)
  }
  if (hess) 
    w$se <- sqrt(diag(-1 * solve(w$hessian)))
  else w$se <- NULL
  wc <- as.paleoTSfit(logL = w$value, parameters = w$par, modelName = paste("URW", 
                                                                            "Trend","URW", sep = "-"), method = "Joint", K = K, n = length(y$mm), 
                      se = w$se)
  return(wc)
}