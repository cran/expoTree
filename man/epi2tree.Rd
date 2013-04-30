\name{epi2tree}
\alias{epi2tree}
\title{
  Convert an 'epi' object from 'sim.epi' to 'phylo' format (ape).
}
\description{
  Converts an 'epi' object from a simulation of 'sim.epi' to 'phylo' format (ape).
}
\usage{
  epi2tree(epi)
}
\arguments{
  \item{epi}{
    'epi' object.
  }
}
\author{
  Gabriel E Leventhal
}
\seealso{
  \code{\link{expoTree}}, \code{\link{ape}}
}
\examples{
  x <- sim.epi(100,1,0.1,0.4,100)
  tree <- epi2tree(x)
}
\keyword{phylogeny}
\keyword{phylo}
\keyword{ape}
