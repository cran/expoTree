\name{sim.epi}
\alias{sim.epi}
\title{
  Simulates epidemic outbreaks and sampled transmission trees.
}
\description{
  Simulate an instance of an epidemic outbreak along with a sampled
  transmission tree.
}
\usage{
  sim.epi(N,beta,mu,psi,max.samples,min.outbreak=min(10,max.samples),max.time=-1.0) 
}
\arguments{
  \item{N}{
    Population size (carrying capacity). Set N=0 for density-independent
    transmission.
  }
  \item{beta}{
    Transmission rate.
  }
  \item{mu}{
    Death rate.
  }
  \item{psi}{
    Sampling rate.
  }
  \item{max.samples}{
    Simulate epidemic until this many individuals are sampled. Set to zero for
    species phylogenies (see Details).
  }
  \item{min.outbreak}{
    Minimum number of infecteds that is considered an outbreak (avoids
    stochastic extinctions).
  }
  \item{max.time}{
    Maximum time to run epidemic, even if max.samples isn't reached.
  }
}
\details{
  Fixed-time sampling.
  The default behavior is to simulate sampled phylogenies for sequentially
  sampled data. It is possible to simulate fixed time sampling by setting
  psi=0. In this case the simulation is run until the total number of
  infecteds is equal to max.samples, when the simulation is stopped and all
  lineages are sampled.
}
\value{
  List with the following entries:

  - times   : event times of the sampled transmission tree.
  - ttypes  : event types of the sampled transmission tree (0 = sampling, 1 = branching).
  - itimes  : infection times of all (sampled and unsampled) individuals.
  - dtimes  : removal times of all individuals.
  - dtypes  : type of removal (0 = sampling, -1 = death, 1 = still infected at
              the present)
  - id      : identifier of the individual
  - parent  : identifier of the infector of an individual (can be used to
              recreate the transmission history)
}
\references{
  Leventhal, Bonhoeffer & Stadler, 2012
}
\author{
  Gabriel E Leventhal
}
\seealso{
  \code{\link{expoTree}}
}
\examples{
  N <- 100
  beta <- 1
  mu <- 0.1
  psi <- 0.4
  nsamp <- 100
  epi <- sim.epi(N,beta,mu,psi,nsamp)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{likelihood}
\keyword{phylogenetics}
