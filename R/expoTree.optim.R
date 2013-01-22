expoTree.optim <- function(forest,lo=rep(5,0),hi=rep(5,0),
      fix=rep(5,F),fix.val=rep(5,0),pars=vector(length=sum(!fix))*1,
      survival=TRUE,method="pso",control=list(trace=0)) 
{
  if (! require(pso)) {
    stop("Particle Swarm Optimization package 'pso' required.")
  }

  if (! is.list(forest)) stop("Must supply a list of trees.")

  fn <- function(par,forest,fix,fix.val,survival) {
    x <- fix.val
    cnt <- 1
    for (i in 1:5) {
      if (! fix[i]) {
        x[i] <- par[cnt]
        cnt <- cnt + 1
      }
    }
    x[1] <- floor(x[1])
    if (x[3] < 0) x[3] <- x[4]*(1.0/(-x[3])-1.0)
    lik <- sapply(forest,function(tree) {
      l <- -Inf
      if (x[1] > 0) {
        l <- runExpoTree(N=x[1],beta=x[2],mu=x[3],psi=x[4],rho=x[5],
                  times=tree[,1],ttypes=tree[,2],survival=survival) 
      } else {
        l <- infExpoTree(beta=x[2],mu=x[3],psi=x[4],rho=x[5],
                  times=tree[,1],ttypes=tree[,2],survival=survival) 
      }
      return(l)
    })
    return(sum(lik))
  }

  control$fnscale <- -1
  opt <- c()
  if (method == "pso") {
    opt <- psoptim(par=pars,fn=fn,lower=lo,upper=hi,control=control,
                   forest=forest,fix=fix,fix.val=fix.val,survival=survival)
  } else if (method == "nelder-mead") {
    opt <- optim(par=pars,fn=fn,control=control,
                 forest=forest,fix=fix,fix.val=fix.val,survival=survival)
  }
  opt
}

