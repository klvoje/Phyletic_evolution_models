logL.joint.multiv.OU<-function (p, yl) 
{
  Smult <- 0
  nseq <- length(yl)
  for (i in 1:nseq) {
    
    Smult <- Smult + logL.joint.OU(p = c(p[i], p[nseq +1], p[i+(nseq +1)], p[length(p)]), yl[[i]])
    
  }
  return(Smult)
}