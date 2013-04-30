pso2gnuplot <- function(opt) {
  max.f <- -min(sapply(opt,min))
  do.call(rbind,lapply(opt,function(i) {
    t(rbind(opt$x[[i]],-opt[[i]]-max.f))
  }))
}

