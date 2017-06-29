sl = paste(letters[1],letters,":0",sep="")
for (i in 2:26){
  ll = paste(letters[i],letters,":0",sep="")
  sl = c(sl,ll)
}
SL = paste(LETTERS[1],LETTERS,":0",sep="")
for (i in 2:26){
  ll = paste(LETTERS[i],LETTERS,":0",sep="")
  SL = c(SL,ll)
}
S1 = paste(LETTERS[1],0:9,":0",sep="")
for (i in 2:26){
  ll = paste(LETTERS[i],1:10,":0",sep="")
  SL = c(SL,ll)
}
s1 = paste(letters[1],0:9,":0",sep="")
for (i in 2:26){
  ll = paste(letters[i],1:10,":0",sep="")
  SL = c(SL,ll)
}
sl=c(sl,SL,S1,s1)

#functions
compphyl <- function(newi,identf,ct){
  #set to extant species to the present time
  identf[,1] = as.character(identf[,1])
  identf[,2] = ct-identf[,2]
  for(i in 1:length(identf[,1])){
    ind = regexpr(identf[i,1],newi)[1] + 2
    newi = paste(substr(newi,1,ind),as.character(identf[i,2]),substring(newi,ind+2),sep="")
  }
  return(newi)
}

drop.fossil <- function (phy, tol = 1e-08)
{
  p = phylo2vectors(phy)
  ct = sum(p$wt)
  n <- Ntip(phy)
  x <- dist.nodes(phy)[n + 1, ][1:n]
  dphy = drop.tip(phy, root.edge = T , which(x < max(x) - tol))
  if(Ntip(dphy)>2){
    p2 = phylo2vectors(dphy)
    if(sum(p2$wt) != ct){
      rem_time = ct - sum(p2$wt)
      p2$wt[1] = p2$wt[1] + rem_time
      dphy = vectors2phylo(p2)
    }
  }
  return(dphy)
}

#experiments
sim.est <- function(n_trees, pars, init_par=c(1.8,0.13,60), seed=0, ct=15, parallel=F){ # simulate a tree, drop fossil, and estimate parameters back after bootstrap reconstruction
  if (seed != 0) set.seed(seed)
  s = sim.tree(lambda0 = pars[1], mu0 = pars[2], K = pars[3], ct=ct)
  st = s$tree
  p <- mle.tree(st)
  sit = s$tree.extant
  trees = sim.srt(tree=sit, pars=p, parallel = parallel, n_trees = n_trees)
  pars = mle.st(trees)
  return(data.frame(real=p, est=pars))
}

get.time <- function(time){
  dif = proc.time()-time
  return(dif[3])
}
