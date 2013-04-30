\name{plotLTT}
\alias{plotLTT}
\title{
Plot lineages through time.
}
\description{
Creates a (or adds to an existing) plot of the lineages-through-time of a
list of phylogenetic trees.
}
\usage{
plotLTT(trees,col=rgb(.5,.5,.5,.5),xlab="Time", ylab="Lineages",log="y",add=FALSE,
        plot=TRUE,plot.tree=NA,type="s",...) 

}
\arguments{
  \item{trees}{
    A list of trees in two-column format. Column 1 is the branching/sampling
    time and column two is the event type (1 = branching, 0 = sampling) 
  }
  \item{col}{
    Color of the lineages.
  }
  \item{xlab}{
    x-axis label.
  }
  \item{ylab}{
    y-axis label.
  }
  \item{log}{
    Use logarithmic axis (defaults to logarithmic y-axis).
  }
  \item{add}{
    Add the lineages to the current plotting device.
  }
  \item{plot}{
    If set to FALSE, only the number of lineages at each time point are
    returned.
  }
  \item{plot.tree}{
    Plot the phylogenetic tree in the background. Only works for one tree.
  }
  \item{type}{
    Line type. Defaults to steps 's'.
  }
  \item{...}{
    Other parameters to pass to 'plot()'.
  }
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
  N <- 100
  beta <- 1
  mu <- 0.1
  psi <- 0.4
  nsamp <- 100

  epis <- lapply(1:10,function(i) sim.epi(N,beta,mu,psi,nsamp))
  trees <- lapply(epis,function(x) cbind(x$times,x$ttypes))
  plotLTT(trees)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{likelihood}
\keyword{phylogenetics}
