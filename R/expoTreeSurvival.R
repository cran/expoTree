expoTreeSurvival <- function(N,beta,mu,psi,rho,torig) 
{
  pars <- as.vector(c(N,beta,mu,psi,rho))
  p <- .Call("expoTreeSurvival",parameters=pars,torig=as.numeric(torig))
  # p = (p0,p1,p2,...,pN)
  return(p[2])
}

