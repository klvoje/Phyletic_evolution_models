opt.joint.Mult.accel_decel<-function (yl, cl = list(fnscale = -1), model = c("accel", "decel"), pool = TRUE,  hess = FALSE) 
{
  if (class(yl) == "paleoTS" || class(yl[[1]]) != "paleoTS") 
    stop("opt.joint.Mult() is only for a list of multiple paleoTS sequences\n")
  nseq <- length(yl)
  if (pool) {
    for (i in 1:nseq) yl[[i]] <- pool.var(yl[[i]], ret.paleoTS = TRUE)
  }
  model <- match.arg(model)
  
  if (model == "accel") {
    ll <- c(rep(NA, nseq), 0, +0.000001)
    uu <- c(rep(NA, nseq), NA, NA)
    anc0 <- sapply(yl, FUN = function(x) x$mm[1])
    P0 <- sapply(yl, FUN = mle.URW)
    p0 <- c(anc0, 1e-07, 1)
    if (p0[nseq + 1] <= 0) 
      p0[nseq + 1] <- 1e-07
    names(p0) <- c(paste("anc", 1:nseq, sep = ""), "vstep", 
                   "r")
    K <- nseq + 2
  }
  
  if (model == "decel") {
    ll <- c(rep(NA, nseq), 0, NA)
    uu <- c(rep(NA, nseq), NA, -0.00001)
    anc0 <- sapply(yl, FUN = function(x) x$mm[1])
    P0 <- sapply(yl, FUN = mle.URW)
    p0 <- c(anc0, 1e-07, -1)
    if (p0[nseq + 1] <= 0) 
      p0[nseq + 1] <- 1e-07
    names(p0) <- c(paste("anc", 1:nseq, sep = ""), "vstep", 
                   "r")
    K <- nseq + 2
  }
  
  if (is.null(cl$ndeps)) 
    cl$ndeps <- rep(1e-08, length(p0))
  w <- try(optim(p0, fn = logL.joint.multiv.accel_decel, method = "L-BFGS-B", 
                 lower = ll, upper = uu, control = cl, hessian = hess, yl = yl), silent = TRUE)
  
  n <- sum(sapply(yl, FUN = function(x) length(x$mm)))
  if (hess) 
    w$se <- sqrt(diag(-1 * solve(w$hessian)))
  else w$se <- NULL
  wc <- as.paleoTSfit(logL = w$value, parameters = w$par, modelName = 
                        paste(model,".Multivariate accel_decel", sep = ""), method = "Joint", K = K, n = n, se = w$se)
  
  
  
  return(wc)
  
}
