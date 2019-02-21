logL.joint.multiv.accel_decel<-function (p, yl) 
{
  Smult <- 0
  nseq <- length(yl)
  for (i in 1:nseq) {
     
    Smult <- Smult + logL.joint.accel_decel(p = c(p[i], p[nseq +1], p[nseq + 2]), yl[[i]])
    
  }
  return(Smult)
}


