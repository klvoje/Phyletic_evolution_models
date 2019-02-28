plot.OU<-function (x, w, np = 500) 
{
  ee <- ii <- array(dim = np)
  mn <- w$modelName
  mp <- w$par
  x0 <- ifelse(w$method == "AD", x$mm[1], w$par["anc"])
  ttp <- seq(x$tt[1], x$tt[length(x$tt)], length.out = np)
  
  
  
  ee <- mp["theta"] * (1 - exp(-mp["alpha"] * ttp)) + mp["anc"] * exp(-mp["alpha"] * ttp)
  vv <- (mp["vstep"]/(2 * mp["alpha"])) * (1 - exp(-2 * mp["alpha"] * ttp))
  
  
  res <- list(tt = ttp, ee = ee, ll = ee - 1.96 * sqrt(vv), 
              uu = ee + 1.96 * sqrt(vv))
  return(res)
}