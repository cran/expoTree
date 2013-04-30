\name{expoTree.optim}
\alias{expoTree.optim}
\title{
Perform optimization to find a maximum likelihood estimate.
}
\description{
Perform optimization to find a maximum likelihood estimate.
}
\usage{
expoTree.optim(forest,fix=rep(5,F),fix.val=rep(5,0),
      pars=vector(length=sum(!fix))*1,lo=rep(sum(!fix),0),hi=rep(sum(!fix),0),
      survival=TRUE,vflag=0,method="Nelder-Mead",control=list(trace=0)) 
}
\arguments{
  \item{forest}{
    List of trees in two-column format. Column 1 are the event times
    (branching or sampling) and column 2 is the event type (0 = sampling, 1 =
    branching).
  }
  \item{fix}{
    Logical vector specifying which parameters to keep constant.
  }
  \item{fix.val}{
    Values for fixed parameters. Also specify dummy values for variable
    parameters. These are ignored.
  }
  \item{pars}{
    Starting values for the parameters.
  }
  \item{lo}{
    Lower bound for parameters.
  }
  \item{hi}{
    Upper bound for parameters.
  }
  \item{vflag}{
    Set verbosity level.
  }
  \item{survival}{
    Condition on the likelihood of observing the tree.
  }
  \item{method}{
    Choose optimization method. Passed on to 'optim'. See 'optim' for details.
  }
  \item{control}{
    Control parameters for the optimization method. See 'optim' for details.
  }
}
\details{
  Special case:
  It is possible to fix the ratio between mu and psi. If mu is fixed at a
  negtive value, it is seen as a fixed ratio.
  Example: 
  r <- 0.5
  fix <- c(F,F,T,F,T)
  fix.val <- c(0,0,-r,0,0)
  At every evaluation, mu is calculated accordingly:
  psi/(psi+mu) = r  ==> mu = psi*(1/r-1)
}
\value{
  Output from the optimization. List with the following entries:
    par        : parameter estimates
    value      : negative log-likelihood
  For the other returned values, see help(optim) respectively
}
\references{
  Leventhal, Guenthard, Bonhoeffer & Stadler, 2013
}
\author{
  Gabriel E. Leventhal
}
\seealso{
  \code{\link{expoTree}}
}
\examples{
  # simulate trees
  N <- 15
  beta <- 1
  mu <- 0.1
  psi <- 0.1
  rho <- 0
  nsamp <- 20
  epi <- sim.epi(N,beta,mu,psi,nsamp)
  tree <- cbind(epi$times,epi$ttypes)

  extant <- sum(2*tree[,2]-1)
  lineages <- sum(2*tree[,2]-1)+cumsum(1-2*tree[,2])
  max.lineages <- max(lineages)

  # calculate likelihood for the forest
  lik <- runExpoTree(pars=c(N,beta,mu,psi,rho),times=tree[,1],ttypes=tree[,2],
                     survival=TRUE,return.full=FALSE)
  cat("Likelihood = ",lik,"\n")

  if (! any(is.nan(lik))) {
    expoTree.optim(list(tree),pars=c(N,beta,psi),
                   lo=c(max(max.lineages),0,0),hi=c(50,3,1),
                   fix=c(FALSE,FALSE,TRUE,FALSE,TRUE),
                   fix.val=c(0,0,-psi/(psi+mu),0,0))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{likelihood}
\keyword{phylogenetics}
