fit.Stasis.OU<-function (y, minb = 7, pool = TRUE, silent = FALSE, hess = FALSE) 
{
  ns <- length(y$mm)
  ng <- 2
  
  ##define shift points (matrix):
  GG <- shifts(ns, ng, minb = minb)
  #Define number of shift points:
  nc <- ncol(GG)
  if (!silent) 
    cat("Total # hypotheses: ", nc, "\n")
  #Create empty list
  wl <- list()
  #create array with length = to switch points and where every entry ) -Inf
  logl <- array(-Inf, dim = nc)
  
  #start loop for estimating maximum likelihood parameters for each data set defined by the switch points
  for (i in 1:nc) {
    if (!silent) 
      cat(i, " ")
    #defines which data point in the time series that belong to each of the two sets
    gg <- shift2gg(GG[, i], ns)
        w <- opt.joint.Stasis.OU(y, gg, pool = pool, hess = hess)
    logl[i] <- w$logL
    wl[[i]] <- w
  }
  
  if (!silent) 
    cat("\n")
  winner <- which.max(logl)
  ww <- wl[[winner]]
  ss <- GG[, winner]
  names(ss) <- paste("shift", 1:(ng - 1), sep = "")
  ww$parameters <- append(ww$parameters, ss)
  ww$all.logl <- logl
  ww$GG <- GG
  return(ww)
}