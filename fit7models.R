fitallmodels<-function (y, silent = FALSE, ...) 
{
  args <- list(...)
  check.var <- TRUE
  if (length(args) > 0) 
    if (args$pool == FALSE) 
      check.var <- FALSE
  if (check.var) {
    tv <- test.var.het(y)
    pv <- round(tv$p.value, 0)
    wm <- paste("Sample variances not equal (P = ", pv, "); consider using argument pool=FALSE", 
                collapse = "")
    if (pv <= 0.05) 
      warning(wm)
  }
    m1 <- opt.joint.GRW(y, ...)
    m2 <- opt.joint.URW(y, ...)
    m3 <- opt.joint.Stasis(y, ...)
    m4 <- opt.joint.OU(y, ...)
    m5 <- opt.joint.decel(y, ...)
    
  mc <- compareModels(m1, m2, m3, m4, m5, silent = silent)
  invisible(mc)
}