MCEM.phylo <- function(tree, init_par, maxNtree=1000, printpar=TRUE, impsam=FALSE, parallel=F){
  print('MCEM computations initiated. This might take a while.')
  n_pars = length(init_par)
  pars = init_par
  Pars = matrix(nrow=(n_it+1),ncol=n_pars)
  i = 1
  qt = 0
  q = vector(mode='numeric',length = n_it)
  Pars[i,] = pars
  n_trees=10
  while(n_trees <= maxNtree){ #dist(rbind(x1, x2)) > tol &
    time = proc.time()
    if(printpar){
      print(paste('iteration #',i,':'))
      print(pars)
      print(paste('previous iteration took',qt,'secs'))
    }
    p = proc.time()
    trees <- sim.srt(tree=tree, pars=pars, n_trees=n_trees, parallel = parallel)
    pars = mle.st(trees,impsam=impsam)
    Pars[(i+1),] = pars
    qt = get.time(time)
    q[i] = qt
    i = i+1
    if(i%%20 == 0){
      n_trees = n_trees*10
    }
  }
  return(list(pars=Pars, times=q))
}

