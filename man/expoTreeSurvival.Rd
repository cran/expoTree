\name{expoTreeSurvival}
\alias{expoTreeSurvival}
\title{
  Calculate the probability that at least one lineage is sampled.
}
\description{
  Calculate the probability that at least one lineage is sampled.
}
\usage{
  expoTreeSurvival(N, beta, mu, psi, rho, torig) 
}
\arguments{
  \item{N}{
    Susceptible population size.
  }
  \item{beta}{
    Branching rate. For epidemic trees, this is the infection rate.
  }
  \item{mu}{
    Extinction rate. For epidemic trees, this is the death/recovery rate.
  }
  \item{psi}{
    Sampling rate.
  }
  \item{rho}{
    Initial sampling rate.
  }
  \item{torig}{
    Time of the origin of the tree (i.e. height of the tree).
  }
}
\value{
  Likelihood
}
\references{
  Tanja Stadler, Journal of Theoretical Biology, 2010
}
\author{
  Gabriel E Leventhal
}

\seealso{
\code{\link{expoTree}}
}
\examples{
  N <- 15
  beta <- 1
  mu <- 0.1
  psi <- 0
  rho <- 1
  torig <- 10.0

  lik <- expoTreeSurvival(N,beta,mu,psi,rho,torig)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{likelihood}
\keyword{phylogenetics}
