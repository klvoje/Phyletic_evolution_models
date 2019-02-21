logL.joint.URW.trend.URW<-function (p, y, gg) 
{
  #These parameters need to match the order of p0
  anc<-p[1]
  vstep_1 <- p[2]
  mstep_trend <- p[3]
  vstep_trend <- p[4]
  vstep_2 <- p[5]

  n <- length(y$mm)
  
  #define a vector based on length of first and second RW in the time series
  rw_1.seg <- which(gg == 1)
  rw_2.seg <- which(gg == 3)
  dt.seg <- which(gg == 2)
  
  #Extract the time vector for the random walks and trend:
  tt.rw_1 <- y$tt[rw_1.seg]
  tt.rw_2 <- y$tt[rw_2.seg]
  tt.trend <- y$tt[dt.seg] - y$tt[dt.seg[1] - 1]

  #create a vector for the expected means for the whole time series
  M <- c(rep(anc, sum(gg==1)), (anc + mstep_trend*tt.trend), rep(tail((anc + mstep_trend*tt.trend),1), sum(gg==3)))
  M <- unname(M)

  #Create variance-covariance matrices for the three models
  VVrw_1 <- vstep_1 * outer(tt.rw_1, tt.rw_1, FUN = pmin)
  VVtrend <- vstep_trend * outer(tt.trend, tt.trend, FUN = pmin)
  VVrw_2 <- vstep_2 * outer(tt.rw_2, tt.rw_2, FUN = pmin)
  
  #Create empty n*n matrix
  VVtot <- array(0, dim = c(n, n))
  
  #Create final variance matrix by combining the variance matrices from stasis and RW 
  VVtot[rw_1.seg, rw_1.seg] <- VVrw_1
  VVtot[dt.seg, dt.seg] <- VVtrend
  VVtot[rw_2.seg, rw_2.seg] <- VVrw_2
  
  #Add population variance to the diagonal of the variance matrix 
  diag(VVtot) <- diag(VVtot) + y$vv/y$nn
  
  #
  S <- dmnorm(y$mm, mean = M, varcov = VVtot, log = TRUE)
  return(S)
}