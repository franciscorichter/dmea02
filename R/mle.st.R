mle.st <- function(st, init_par = c(2,1,60), impsam=F){
  n_trees = length(st)
  if(impsam==T){
    W = vector(mode = 'numeric',length = n_trees)
    D = vector(mode = 'numeric',length = n_trees)
    for(j in 1:n_trees){
      rec = st[[j]]
      W[j] = rec$weight
      D[j] = length(rec$wt)
    }
    m1=lm(W~D)
    p = subplex(par = init_par, fn = llik.st, setoftrees=st, impsam=impsam, correction=m1$coefficients[2])$par
  }else{
    p = subplex(par = init_par, fn = llik.st, setoftrees=st,impsam=impsam)$par
  }

  return(p)
}
