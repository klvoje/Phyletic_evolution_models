logL.joint.accel_decel<-function (p, y) 
{
  anc <- p[1]
  vs <- p[2]
  r <- p[3]
  n <- length(y$mm)
  VV <- vs * outer((exp(r*y$tt) -1)/r, (exp(r*y$tt) -1)/r, FUN = pmin)
  diag(VV) <- diag(VV) + y$vv/y$nn
  M <- rep(anc, n)
  S <- dmnorm(y$mm, mean = M, varcov = VV, log = TRUE)
  return(S)
}