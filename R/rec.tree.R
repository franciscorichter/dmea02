rec.tree <- function(tree,pars,model='dd',seed=0){
  if(seed>0) set.seed(seed)
  wt = tree$wt
  ct = sum(wt)
  lambda0=pars[1]
  mu0=pars[2]
  K=pars[3]
  limit = lambda0*K/(lambda0-mu0)
  rs = dim = length(wt)
  if (limit < (dim-1)){
    print('parameters does not make sense, observed tree implies negative rates')
  }
  ms = NULL # missing species, for now we just add time. When we consider topology we do it with species as well
  e.lims = NULL # limits on extinctions
  cbt = 0
  N = 1
  nm= 0 # number of missing species
  prob = NULL
  gprob=0 # gost probability
  for(i in 1:dim){
    rs = rs-1
    cwt = wt[i]
    cbt = sum(wt[0:(i-1)])
    key = 0
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
        t.spe = rexp(1,s)  # waiting time of iteration i
      }
      if(nm > 0){ # if there are missing species
        t.ext = NULL
        for(j in 1:nm){
          t.ext[j] = rtrunc(1,'exp',a=0,b=(e.lims[j]-cbt),rate=mu)
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
      if(mint < cwt){
        if(mint == t.spe){#speciation
          acep = runif(1)
          thre = pexp(ct-cbt,mu)
          if(acep < thre){
            ms = c(ms,cbt+t.spe)
            nm = nm + 1 # number of missing species
            if((N + rs-1) < limit){
              e.lims = c(e.lims,ct)
            }
            else{
              e.lims = c(e.lims,sum(wt[1:(dim - (N+rs-floor(limit)-1))])) # cual es la interpretacion de esta formula?
            }
            N = N+1
            #print(paste('at branching time',cbt+t.spe,'a missing species arises, resulting on',N,'current species, happening before the current waiting time',cwt))
          }
          cwt = cwt - t.spe
          cbt = cbt + t.spe
        }
        else{#extinction
          pickone = sample(1:nm,1)
          t_spe = ms[pickone]
          t_ext = cbt + t.ext
          tree = update.tree(tree,t_spe=t_spe,t_ext=t_ext)
          ms = ms[-pickone]
          cwt = cwt - t.ext
          cbt = cbt + t.ext
          N = N-1
          nm = nm - 1
          e.lims = e.lims[-extinctedone]
         }
      }
      else{
         key = 1
      }
    }
    N= N+1
  }
  tree$prob = prob
  return(tree)
}



dists <- function(wt,lambda,mu=NULL,b=NULL,log=FALSE){
  nm = length(b) # number of missing species
  s = sum(lambda)
  if(nm==0){
    product = 1
  }
  else{
    product = prod((exp(-wt*mu)-exp(-mu*b))/(1-exp(-mu*b)))
  }
  prob = exp(-wt*s)*product
  if(prob>1) print('prob grater than 1???')
  if(log){
    prob = log(prob)
  }
  return(prob)
}
