mle.tree <- function(tree,init_par = c(2,1,60)){
  p = subplex(par = init_par, fn = llik, tree=tree)$par
  return(p)
}
