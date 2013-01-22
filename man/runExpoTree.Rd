\name{runExpoTree}
\alias{runExpoTree}
\title{
Density-dependent likelihood caluculation
}
\description{
Calculates the density dependent likelihood of a phylogenetic tree. 
It takes branching and sampling times as an argument and integrates the
likelihoood function over the whole tree.
}
\usage{
runExpoTree(N, beta, mu, psi, rho, times, ttypes, survival = TRUE)
}
\arguments{
  \item{N}{
    Total population size. Must be at least the maximal number of extant lineages
    at any given time.
  }
  \item{beta}{
    Branching rate. For epidemic trees, this is the infection rate. For species
    trees, this is the speciation rate.
  }
  \item{mu}{
    Extinction rate. For epidemic trees, this is the death/recovery rate. For species
    trees, this is the extinction rate.
  }
  \item{psi}{
    Sampling rate.
  }
  \item{rho}{
    Initial sampling rate.
  }
  \item{times}{
    Vector of event times.
  }
  \item{ttypes}{
    Vector of event types.
  }
  \item{survival}{
    Condition on the likelihood of observing the tree.
  }
}
\value{
  Log-Likelihood
}
\references{
  Leventhal, Bonhoeffer, Guenthard & Stadler, 2012
}
\author{
  Gabriel E Leventhal
}
\seealso{
  \code{\link{expoTree}}
}
\examples{
times <- cumsum(runif(10))
ttypes <- rep(1,10)

N <- 15
beta <- 1
mu <- 0.1
psi <- 0
rho <- 1

lik <- runExpoTree(N,beta,mu,psi,rho,times,ttypes)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{likelihood}
\keyword{phylogenetics}
