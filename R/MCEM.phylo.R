MCEM.phylo <- function(tree, init_par, cutTime=10, printpar=TRUE, impsam=FALSE, parallel=F, ips=20){
  print('MCEM computations initiated. This might take a while.')
  TIME = proc.time()
  n_pars = length(init_par)
  pars = init_par
  Pars = matrix(nrow=(1000),ncol=n_pars)
  i = 1
  qt = 0
  q = vector(mode='numeric',length = 1000)
  Pars[i,] = pars
  n_trees=50
  while(get.time(TIME,mode='min') < cutTime){
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
    if(i%%ips == 0){
      n_trees = n_trees + 100
    }
  }
  Pars = Pars[1:i,]
  q = q[1:i]
  return(list(pars=Pars, times=q))
}

