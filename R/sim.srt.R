sim.srt <- function(tree, pars, parallel=F, n_trees, WT=FALSE){    # simulate set of reconstructed trees
  if(parallel){
    no_cores <- detectCores()- 1
    cl <- makeCluster(no_cores)
    registerDoParallel(cl)
    trees <- foreach(i = 1:n_trees, combine = list) %dopar% dmea02::rec.tree(tree=tree, pars=pars)
    stopCluster(cl)
  }
  else{
    trees = vector('list',length=n_trees)
    for (i in 1:n_trees){
      if(method == 1) rec = rec.tree(tree=tree, pars=pars)
      if(WT){trees[[i]] = rec$wt}
      else{trees[[i]] = rec}
    }
  }
  return(trees)
}
