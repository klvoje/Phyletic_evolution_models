logL.joint.URW.URW.URW<-function (p, y, gg) 
{
  #These parameters need to mach the order of p0
  anc<-p[1]
  vs_1 <- p[2]
  vs_2 <- p[3]
  vs_3 <- p[4]
  
  n <- length(y$mm)
  #creates a vector that repeates the ancestral value (first sample mean in sequence) 
  M <- rep(anc, n)
  M <- unname(M)
  
  #define a vector based on length of first and second RW in the time series
  rw_1.seg <- which(gg == 1)
  rw_2.seg <- which(gg == 2)
  rw_3.seg <- which(gg == 3)
  
  #Extract the time vector for the random walks:
  tt.rw_1 <- y$tt[rw_1.seg]
  tt.rw_2 <- y$tt[rw_2.seg]  
  tt.rw_3 <- y$tt[rw_3.seg]


  #Multiply variance parameter (vstep) from RW with outer product of time*time 
  VVrw_1 <- vs_1 * outer(tt.rw_1, tt.rw_1, FUN = pmin)
  VVrw_2 <- vs_2 * outer(tt.rw_2, tt.rw_2, FUN = pmin)
  VVrw_3 <- vs_3 * outer(tt.rw_3, tt.rw_3, FUN = pmin)
  
  #Create empty n*n matrix
  VVtot <- array(0, dim = c(n, n))
  
  #Create final variance matrix by combining the variance matrices from stasis and RW 
  VVtot[rw_1.seg, rw_1.seg] <- VVrw_1
  VVtot[rw_2.seg, rw_2.seg] <- VVrw_2
  VVtot[rw_3.seg, rw_3.seg] <- VVrw_3
  
  #Add population variance to the diagonal of the variance matrix 
  diag(VVtot) <- diag(VVtot) + y$vv/y$nn
  
  #
  S <- dmnorm(y$mm, mean = M, varcov = VVtot, log = TRUE)
  return(S)
}