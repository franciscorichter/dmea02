sim.est <- function(n_trees, ct=15, pars, init_par=c(1.8,0.13,60), impsam=FALSE, seed=0, parallel=F, method=1){ # simulate a tree, drop fossil, and estimate parameters back after bootstrap reconstruction
  if (seed != 0) set.seed(seed)
  st = sim.tree(lambda0 = pars[1], mu0 = pars[2], K = pars[3], ct=ct)
  p <- mle.tree(st$tree)
  sit = st$tree.extant
  trees = sim.srt(tree=sit, pars=p, parallel = parallel, n_trees = n_trees,method=method)
  pars = mle.st(st=trees,init_par)
  return(data.frame(real=p, est=pars))
}
