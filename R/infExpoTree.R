infExpoTree <- function(beta,mu,psi,rho,times,ttypes,survival=TRUE) 
{
  pars <- as.vector(c(beta,mu,psi,rho))
  f <- .Call("infTreeEval",parameters=pars,
      times=as.numeric(times),ttypes=as.integer(ttypes),survival=as.integer(survival))
  return(f)
}

