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

get.time <- function(time,mode='sec'){
  dif = proc.time()-time
  ti = as.numeric(dif[3])
  if(mode == 'min')  ti = ti/60
  if(mode == 'hou') ti = ti/3600
  return(ti)
}


inverse = function (f, lower = -100, upper = 100) {
  function (y) uniroot((function (x) f(x) - y), lower = lower, upper = upper)[1]
}

pms <- function(t,sumspec,mu,rem){
  pms = sumspec*exp(-sumspec*t)*(1-exp(-(rem-t)*mu))
  return(pms)
}

sumspec = 10
remt=5
mu = 0.1
pms2 <- function(t){
  return(pms(t,sumspec=sumspec,mu=mu,rem=remt))
}
missspec <- inverse(pms2, lower = 0, upper = 1)
missspec(0.2)



par_est_vis <- function(P,par,PR,true){
  # P is the recost estim values
  # par is the parameter you want to use
  # PR is the real estim
  if (par == 1){
    int = true
    parname = 'lambda'
    bin = 0.02
  }
  if (par==2){
    int= true
    parname = 'mu'
    bin = 0.02
  }
  if (par == 3){
    int = true
    parname = 'K'
    n = dim(P)[1]
    P1 = P[P[,3]<100 & PR[,3]<100,]
    PR = PR[P[,3]<100 & PR[,3]<100,]
    P = P1
    n2 = dim(P)[1]
    print(paste(1-n2/n,'proportion of data was excluded for vizualization purposes'))
    bin = 2
  }
  if (par ==4){
    int = 40
    parname = 'K'
    P = P[P[,4]<100,]
  }
  hist_top <- ggplot()+geom_histogram(aes(P[,par]),binwidth=bin) + geom_vline(xintercept=int) + xlab('MLE from incomplete tree')
  empty <- ggplot()+geom_point(aes(1,1), colour="white")+
    theme(axis.ticks=element_blank(),
          panel.background=element_blank(),
          axis.text.x=element_blank(), axis.text.y=element_blank(),
          axis.title.x=element_blank(), axis.title.y=element_blank())

  scatter <- ggplot()+geom_point(aes(P[,par], PR[,par]))+ geom_abline(intercept = 0, slope = 1) + ylab(TeX(paste('$\\hat{\\',parname,'}_{C}$',sep='')))+xlab(TeX(paste('$\\hat{\\',parname,'}_{I}$',sep='')))
  hist_right <- ggplot()+geom_histogram(aes(PR[,par]),binwidth=bin)+coord_flip()+ geom_vline(xintercept=int) +xlab('MLE of complete tree')

  grid.arrange(hist_top, empty, scatter, hist_right, ncol=2, nrow=2, widths=c(4, 1), heights=c(1, 4))
}


count.missing <- function(st){
  total = length(st$tree$wt)
  extant = length(st$tree.extant$wt)
  missing = (total - extant)/2
  return(missing)
}
