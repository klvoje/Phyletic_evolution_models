opt.joint.Mult.OU<-function (yl, cl = list(fnscale = -1), pool = TRUE,  hess = FALSE) 
{
  if (class(yl) == "paleoTS" || class(yl[[1]]) != "paleoTS") 
    stop("opt.joint.Mult() is only for a list of multiple paleoTS sequences\n")
  nseq <- length(yl)
  if (pool) {
    for (i in 1:nseq) yl[[i]] <- pool.var(yl[[i]], ret.paleoTS = TRUE)
  }
 
    ll <- c(rep(NA, nseq), 1e-10, rep(NA, nseq), 1e-08)
    uu <- c(rep(NA, nseq), NA, rep(NA, nseq), NA)
    anc0 <- sapply(yl, FUN = function(x) x$mm[1])
    w0 <- sapply(yl, FUN = mle.GRW)
    theta0<-sapply(yl, FUN = function(y) y$mm[length(y$mm)])
    halft <- (yl[[1]]$tt[length(yl[[1]]$tt)] - yl[[1]]$tt[1])/4
    p0 <- c(anc0, mean(w0[2,])/10, theta0, halft)
    if (p0[nseq + 1] <= 0) 
      p0[nseq + 1] <- 1e-07
    names(p0) <- c(paste("anc", 1:nseq, sep = ""), "vstep", paste("theta", 1:nseq, sep = ""), "alpha")
    K <- nseq + 5
 
  
  if (is.null(cl$ndeps)) 
    cl$ndeps <- rep(1e-08, length(p0))
  w <- try(optim(p0, fn = logL.joint.multiv.OU, method = "L-BFGS-B", 
                 lower = ll, upper = uu, control = cl, hessian = hess, yl = yl), silent = TRUE)
  
  n <- sum(sapply(yl, FUN = function(x) length(x$mm)))
  if (hess) 
    w$se <- sqrt(diag(-1 * solve(w$hessian)))
  else w$se <- NULL
  model<-"OU"
  wc <- as.paleoTSfit(logL = w$value, parameters = w$par, modelName = 
                        paste(model,".Multivariate OU", sep = ""), method = "Joint", K = K, n = n, se = w$se)
  
  
  
  return(wc)
  
}
