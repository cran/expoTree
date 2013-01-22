plotLTT <- function(trees,col=rgb(.5,.5,.5,.5),xlab="Time",
                     ylab="Lineages",log="y",add=FALSE,...) 
{
  max.time <- sapply(trees,function(t) max(t[,1]))
  extant <- sapply(trees,function(t) sum(2*t[,2]-1))
  lineages <- lapply(trees,function(t) sum(2*t[,2]-1)+cumsum(1-2*t[,2]))
  if (! add) {
    plot(1,type="n",xlim=c(-max(max.time),0),ylim=c(1,max(unlist(lineages))),log=log,
         xlab=xlab,ylab=ylab,...)
  }
  for (i in 1:length(lineages)) {
    lines(-trees[[i]][,1],lineages[[i]],col=col,type="s")
  }
}

