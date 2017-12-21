mle.st <- function(st, init_par = c(2,1,60), impsam=F, correction = TRUE){
  n_trees = length(st)
  if(impsam==T){
    W = vector(mode = 'numeric',length = n_trees)
    D = vector(mode = 'numeric',length = n_trees)
    for(j in 1:n_trees){
      rec = st[[j]]
      W[j] = rec$logweight
      D[j] = length(rec$wt)
    }
    if(correction==TRUE){
      m1=lm(W~D)
      m1=m1$coefficients[2]
    }
    else{
      m1=0
    }
    p = subplex(par = init_par, fn = llik.st, setoftrees=st, impsam=impsam, correction=m1)$par
  }else{
    p = subplex(par = init_par, fn = llik.st, setoftrees=st,impsam=impsam)$par
  }

  return(p)
}
