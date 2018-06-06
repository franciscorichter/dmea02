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


cuttree <- function(tree,dim){
  if(dim>length(tree$E)){
    stop('dimension greater than tree, so no need to cut')
  }
  wt = tree$wt[1:(dim+1)]
  E = tree$E[1:(dim)]
  n = tree$n[1:(dim+1)]
  newtree = list(wt=wt,E=E,n=n)
  return(newtree)
}

update.tree <- function(tree,t_spe,t_ext){
  wt = tree$wt
  #n = tree$n
  E = tree$E
  ct = sum(wt)
  if(t_ext > ct){
    stop('Extinction beyond present!')
  }
  if(t_ext < t_spe){
    stop('Speciation after extinction!')
  }
  # speciation
  K = length(wt)
  k = length(wt[cumsum(wt) < t_spe])
  #n = c(n[1:(k+1)],(n[(k+1):K]+1))
  if((k+1)<K){
    lastbit = E[(k+1):(K-1)]
  }
  else{
    lastbit = NULL
  }
  E = c(E[0:k],1,lastbit)
  # (this can be witten in one line, is convinent?)
  if(k==0){
    wt = c(t_spe, wt[1]-t_spe, wt[2:K])
  }
  if(k > 0 & k < (K-1)){
    wt = c(wt[1:k], t_spe - sum(wt[1:k]), wt[k+1]-(t_spe - sum(wt[1:k])), wt[(k+2):K])
  }
  if(k == (K-1)){
    wt = c(wt[1:(K-1)],t_spe-sum(wt[1:(K-1)]),ct-t_spe)
  }
  #extinction
  K = length(wt)
  k = length(wt[cumsum(wt) < t_ext])
  #n = c(n[1:(k+1)],(n[(k+1):K]-1))
  if((k+1)<K){
    lastbit = E[(k+1):(K-1)]
  }
  else{
    lastbit = NULL
  }
  E = c(E[0:k],0,lastbit)
  # (this can be witten in one line, is convinent?)
  if(k==0){
    wt = c(t_ext, wt[1]-t_ext, wt[2:K])
  }
  if(k > 0 & k < (K-1)){
    wt = c(wt[1:k], t_ext - sum(wt[1:k]), wt[k+1]-(t_ext - sum(wt[1:k])), wt[(k+2):K])
  }
  if(k == (K-1)){
    wt = c(wt[1:(K-1)],t_ext-sum(wt[1:(K-1)]),ct-t_ext)
  }
  #n = c(1,cumsum(E)+cumsum(E-1))
  tree = list(wt=wt,E=E)
  return(tree)
}
sampprob <- function(t,s,mu,r){  ## equation (*)
  term1 = s*(1-exp(-mu*(r-t)))
  c = exp(-(s/mu)*exp(-mu*r))
  if(c==0){
    term2 = 0
  }else{
    term2 = c^{-exp(mu*t)+1}
  }
  f = term1*exp(-s*t)*term2
  return(f)
}

# simulation of missing part [NEED TO CLEAN AND SIMPLIFY]
sim.extinct <- function(tree,pars,model='dd',seed=0, adjustment=FALSE){
  if(seed>0) set.seed(seed)
  wt = tree$wt
  ct = sum(wt)
  tree$E = rep(1,(length(wt)-1))
  lambda0=pars[1]
  mu0=pars[2]
  K=pars[3]
  limit = lambda0*K/(lambda0-mu0)
  rs = dim = length(wt)
  if (limit < (dim-1)){
    print('parameters do not make sense, observed tree implies negative rates')
  }
  ms = NULL # missing species, for now we just add time. When we consider topology we do it with species as well
  e.lims = NULL # limits on extinctions
  cbt = 0
  N = 1
  nm= 0 # number of missing species
  rprob = NULL # true probability of Missing|observed
  sprob = NULL # sampling probability of Missing|observed
  et = NULL # event type
  h = 1 # index to fill probabilities
  dif1 = vector(mode = 'numeric',length = dim)
  dif2 = vector(mode = 'numeric', length = dim)
  for(i in 1:dim){
    rs = rs-1 # this are the remaning speciations
    cwt = wt[i]
    cbt = sum(wt[0:(i-1)])
    key = 0
    last='nothing'
    while(key == 0){
      if(model == "dd"){  # diversity-dependence model
        lambda = max(0,lambda0 - (lambda0-mu0)*N/K)
        mu = mu0
        lambda = rep(lambda,N)
      }
      if(model == 'cr'){ # constant-rate model
        lambda = rep(lambda0,N)
        mu = rep(mu0,N)
      }
      s = sum(lambda)
      if (s == 0){
        t.spe = Inf
      }
      else{
        t.spe = rexp(1,s)
      }
      if(nm > 0){ # if there are missing species simulate extinction times
        t.ext = vector(mode = 'numeric',length = nm)
        for(j in 1:nm){
          t.ext[j] = truncdist::rtrunc(1,'exp',a=0,b=(e.lims[j]-cbt),rate=mu)
        }
        extinctedone = which(t.ext == min(t.ext))
        t.ext = min(t.ext)
      }
      else{
        t.ext = Inf
      }
      if(t.ext ==Inf & t.spe == Inf){
        print('try another K!')
      }
      mint = min(t.spe,t.ext)
      if(nm > 0){
        probs = vector(mode = 'numeric',length = nm)
        for(j in 1:nm){
          probs[j] = 1-truncdist::ptrunc(mint,'exp',a=0,b=(e.lims[j]-cbt),rate=mu)
        }
      }
      else{probs = 1}
      if(mint < cwt){
        if(mint == t.spe & mint != t.ext){#speciation
          acep = runif(1)
          thre = pexp(ct-cbt,mu)
          if(acep < thre){
            ms = c(ms,cbt+t.spe)
            rprob[h] = dexp(x = t.spe, rate = (s+nm*mu))*(s/(s+nm*mu))
            sp = sampprob(t = mint,s = s,mu = mu,r = ct-cbt)#/integrate(sampprob,lower = 0, upper = ct-cbt,s=s,mu=mu,r=ct-cbt)$value
            sprob[h] = prod(probs)*sp
            et[h] = 'speciation'
            h = h + 1
            nm = nm + 1 # number of missing species
            if((N + rs-1) < limit){
              e.lims = c(e.lims,ct)
            }
            else{
              e.lims = c(e.lims,sum(wt[1:(dim - (N+rs-floor(limit)-1))])) # cual es la interpretacion de esta formula?
            }
            N = N+1
            #print(paste('at branching time',cbt+t.spe,'a missing species arises, resulting on',N,'current species, happening before the current waiting time',cwt))
            last = 'speciation'
          }
          else{last = 'nothing'}
          cwt = cwt - t.spe
          cbt = cbt + t.spe
        }
        else{#extinction
          pickone = sample(1:nm,1)
          t_spe = ms[pickone]
          t_ext = cbt + t.ext
          tree = update.tree(tree,t_spe=t_spe,t_ext=t_ext)
          rprob[h] = dexp(x = mint, rate = (s+nm*mu))*(mu/(s+nm*mu))
          et[h] = 'extinction'
          probs = probs[-extinctedone]
          sprob[h] = prod(probs)*truncdist::dtrunc(mint,'exp',a=0,b=(e.lims[extinctedone]-cbt),rate=mu)*(1-integrate(sampprob,lower = 0, upper = mint,s=s,mu=mu,r=ct-cbt)$value)#/integrate(sampprob,lower = 0, upper = ct-cbt,s=s,mu=mu,r=ct-cbt)$value)
          ms = ms[-pickone]
          cwt = cwt - t.ext
          cbt = cbt + t.ext
          N = N-1
          h = h+1
          nm = nm - 1
          e.lims = e.lims[-extinctedone]
          last = 'extinction'
        }
      }
      else{
        key = 1
        rprob[h] = pexp(q = cwt, rate = (s+nm*mu0),lower.tail = FALSE)
        et[h] = 'nothing'
        sprob[h] = prod(probs)*(1 - integrate(Vectorize(sampprob),lower = 0, upper = cwt,s=s,mu=mu,r=ct-cbt)$value)
        h = h+1
        dif1[i] = cwt/wt[i]
        dif2[i] = (mint - cwt)/wt[i]
        if(adjustment & cwt<(mint - cwt) & last=='speciation'){ #Adjusting last speciation
          ms = ms[-length(ms)]
          e.lims = e.lims[-length(e.lims)]
          nm = nm - 1
          N = N-1
        }
      }
    }
    N= N+1
  }
  tree$rprob = rprob
  tree$sprob = sprob
  tree$et = et
  tree$weight = prod(rprob)/prod(sprob)
  logweight = log(rprob)-log(sprob)
  tree$logweight = sum(logweight)
  E = tree$E
  n = c(1,1+cumsum(E)+cumsum(E-1))
  tree$n = n
  return(tree)
}
#negative logLikelihood of a tree
nllik.tree = function(pars,tree){
  b = c(pars[1],(pars[1]-pars[2])/pars[3],pars[2])
  ldt = tree$wt[length(tree$wt)]
  dt = tree$wt[1:(length(tree$wt)-1)]
  lastn = tree$n[length(tree$n)]
  n = tree$n[1:(length(tree$n)-1)]
  E = tree$E
  sigma = n*(b[1]-b[2]*n + b[3]) #n-dimentional
  rho = pmax(b[1]*E-b[2]*n*E+b[3]*(1-E),0)
  l = -sum(-sigma*dt+log(rho)-lastn*ldt)
  if(min(b)<0){l = Inf}
  return(l)
}
# negative logLikelihood of a set of trees
nllik.st = function(pars, st){
  m = length(st$rec)
  l = vector(mode = 'numeric',length = m)
  w = vector(mode = 'numeric',length = m)
  for(i in 1:m){
    s = st$rec[[i]]
    w[i] = st$w[i]
    l[i] = nllik.tree(pars,tree=s)
  }
  L = sum(l*w)
  return(L)
}
# relative likelihood
rel.llik <- function(S1,p0,p1){
  m = length(S1)
  f1 = vector(mode='numeric',length = m)
  f2 = vector(mode='numeric',length = m)
  d = vector(mode='numeric',length = m)
  for(i in 1:m){
    s = S1[[i]]
    f1[i] = nllik.tree(pars=p1,tree=s)
    f2[i] = nllik.tree(pars=p0,tree=s)
    d[i] = length(s$tree$wt)
    if(is.na(f1[i])) print(s)
  }
  Delta = -log(sum(f1/f2)/m)
  return(Delta)
}
# MLE for a set of trees
mle.st <- function(S,init_par = c(0.5,0.5,100)){
  po = subplex(par = init_par, fn = nllik.st, st=S,hessian = TRUE)
  return(po)
}
