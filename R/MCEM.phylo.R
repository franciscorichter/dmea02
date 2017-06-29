MCEM.phylo <- function(tree, init_par, n_trees=10, n_it=30, printpar=TRUE, impsam=FALSE, tol=0.0001, parallel=F){
  n_pars = length(init_par)
  pars = init_par
  Pars = matrix(nrow=(n_it+1),ncol=n_pars)
  i = 1
  qt = c(0,0,0)
  q = vector(mode='numeric',length = n_it)
  Pars[i,] = pars
  while(i <= n_it){ #dist(rbind(x1, x2)) > tol &
    if(printpar){
      print(paste('iteration #',i,':'))
      print(pars)
    }
    p = proc.time()
    trees <- sim.srt(tree=tree, pars=pars, n_trees=n_trees, parallel = parallel)
    pars = mle.st(trees,init_par = init_par)
    Pars[(i+1),] = pars
    qt = proc.time()-p
    q[i] = qt[3]
    i = i+1
  }
  return(list(pars=Pars, times=q))
}

