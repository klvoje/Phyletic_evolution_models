logL.joint.URW.trend.Stasis<-function (p, y, gg) 
{
  #These parameters need to match the order of p0
  anc <- p[1]
  vstep_URW <- p[2]
  mstep <- p[3]
  vstep <- p[4]
  theta <- p[5]
  omega <- p[6]
  n <- length(y$mm)
  
  #define a vector based on length of first and second RW in the time series
  rw.seg <- which(gg == 1)
  dt.seg <- which(gg == 2)
  st.seg <- which(gg == 3)
  
  #Extract the time vector for the trend:
  tt.rw <- y$tt[rw.seg]
  tt.trend <- y$tt[dt.seg] - y$tt[dt.seg[1] - 1]
  
  #create a vector for the expected means for the whole time series
  M <- c(rep(anc, length(rw.seg)), anc + mstep * tt.trend, rep(theta, length(st.seg)))
  M <- unname(M)

  #Create variance-covariance matrices for the three models
  VVrw <- vstep_URW * outer(tt.rw, tt.rw, FUN = pmin)
  VVtrend <- vstep * outer(tt.trend, tt.trend, FUN = pmin)
  VVst <- diag(omega, nrow = length(st.seg))
  
  #Create empty n*n matrix
  VVtot <- array(0, dim = c(n, n))
  
  #Create final variance matrix by combining the variance matrices from stasis and RW 
  VVtot[rw.seg, rw.seg] <- VVrw
  VVtot[dt.seg, dt.seg] <- VVtrend
  VVtot[st.seg, st.seg] <- VVst
  
  #Add population variance to the diagonal of the variance matrix 
  diag(VVtot) <- diag(VVtot) + y$vv/y$nn
  
  #
  S <- dmnorm(y$mm, mean = M, varcov = VVtot, log = TRUE)
  return(S)
}