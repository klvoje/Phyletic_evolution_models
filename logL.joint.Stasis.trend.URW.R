logL.joint.Stasis.trend.URW<-function (p, y, gg) 
{
  #These parameters need to match the order of p0
  theta <- p[1]
  omega <- p[2]
  mstep <- p[3]
  vstep <- p[4]
  vstep_URW <- p[5]

  n <- length(y$mm)
  
  #define a vector based on length of first and second RW in the time series
  st.seg <- which(gg == 1)
  dt.seg <- which(gg == 2)
  rw.seg <- which(gg == 3)
  
  
  #Extract the time vector for the trend:
  tt.trend <- y$tt[dt.seg] - y$tt[dt.seg[1] - 1]
  tt.rw <- y$tt[rw.seg]
  
  #create a vector for the expected means for the whole time series
  M <- c(rep(theta, length(st.seg)), theta + mstep * tt.trend, rep(tail((theta+ mstep * tt.trend),1), sum(gg==3)))
  M <- unname(M)

  #Create variance-covariance matrices for the three models
  VVst <- diag(omega, nrow = length(st.seg))
  VVtrend <- vstep * outer(tt.trend, tt.trend, FUN = pmin)
  VVrw <- vstep_URW * outer(tt.rw, tt.rw, FUN = pmin)
  
  #Create empty n*n matrix
  VVtot <- array(0, dim = c(n, n))
  
  #Create final variance matrix by combining the variance matrices from stasis and RW 
  VVtot[st.seg, st.seg] <- VVst
  VVtot[dt.seg, dt.seg] <- VVtrend
  VVtot[rw.seg, rw.seg] <- VVrw
  
  #Add population variance to the diagonal of the variance matrix 
  diag(VVtot) <- diag(VVtot) + y$vv/y$nn
  
  #
  S <- dmnorm(y$mm, mean = M, varcov = VVtot, log = TRUE)
  return(S)
}