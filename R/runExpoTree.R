runExpoTree <- function(N,beta,mu,psi,rho,times,ttypes,survival=TRUE) 
{
  pars <- as.vector(c(N,beta,mu,psi,rho))
  p <- .Call("expoTreeEval",parameters=pars,
      times=as.numeric(times),ttypes=as.integer(ttypes))
  # p = (p0,p1,p2,...,pN)
  lik <- p[2]
  if (survival) {
    surv <- expoTreeSurvival(N,beta,mu,psi,rho,max(times))
    lik <- lik - surv
  }
  return(lik)
}


