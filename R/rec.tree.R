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
  rprob = NULL # true probability of Missing|observed
  sprob = NULL # sampling probability of Missing|observed
  h = 1
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
        # calculate F(Delta t)
        maxprob = pmiss(t=cwt,s = s,mu = mu,r = (ct-cbt))
        specp = runif(1)
#
        if(maxprob>specp){
#          t.spe = rexp(specp,rate=s)
          t.spe = uniroot(f = pmiss, s = s, mu = mu, r = (ct-cbt), shift = specp, lower = 0, upper = cwt)$root  # waiting time of iteration i
         # print(paste('speciation at',t.spe))
        }else{
          t.spe = 99999999
        #t.spe = rexp(1,s)
        #prob = pmiss(t=cwt,s = s,mu = mu,r = (ct-cbt))
        #unif = runif(1)
        #if(prob > unif){
        #  t.spe = 9999999999
        }
      }
      if(nm > 0){ # if there are missing species
        t.ext = vector(mode = 'numeric',length = nm)
        for(j in 1:nm){
          t.ext[j] = rtrunc(1,'exp',a=0,b=(e.lims[j]-cbt),rate=mu)
          #probs[j] = 1 - ptrunc()
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
          probs[j] = 1-ptrunc(mint,'exp',a=0,b=(e.lims[j]-cbt),rate=mu)
        }
      }
      else{probs = 1}
      if(mint < cwt){
        if(mint == t.spe & mint != t.ext){#speciation
#          acep = runif(1)
#          thre = pexp(ct-cbt,mu)
#          if(acep < thre){
            ms = c(ms,cbt+t.spe)
            rprob[h] = dexp(x = t.spe, rate = (s+nm*mu))*(lambda[1]/(s+nm*mu))
            sprob[h] = prod(probs)*dexp(x = t.spe, rate = (s+nm*mu))*(lambda[1]/(s+nm*mu))
            #print(paste('spec',sprob[h]))
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
          #}
          cwt = cwt - t.spe
          cbt = cbt + t.spe
        }
        else{#extinction
          pickone = sample(1:nm,1)
          t_spe = ms[pickone]
          t_ext = cbt + t.ext
          tree = update.tree(tree,t_spe=t_spe,t_ext=t_ext)
          rprob[h] = dexp(x = mint, rate = (s+nm*mu))*(mu0/(s+nm*mu0))
          probs = probs[-extinctedone]
          sprob[h] = prod(probs)*dtrunc(mint,'exp',a=0,b=(e.lims[extinctedone]-cbt),rate=mu)*pexp(mint,rate = s,lower.tail = FALSE)*(mu0/(s+nm*mu0))
          ms = ms[-pickone]
          cwt = cwt - t.ext
          cbt = cbt + t.ext
          N = N-1
          #print(paste('ext',sprob[h]))
          h = h+1
          nm = nm - 1
          e.lims = e.lims[-extinctedone]
         }
      }
      else{
         key = 1
         rprob[h] = pexp(q = cwt, rate = (s+nm*mu0),lower.tail = FALSE)
         sprob[h] = prod(probs)*pexp(cwt,rate = s,lower.tail = FALSE)
         #print(paste('nothing',sprob[h]))
         h = h+1
      }
    }
    N= N+1
  }
  tree$rprob = rprob
  tree$sprob = sprob
  tree$weight = prod(rprob)/prod(sprob)
  return(tree)
}

pmiss <- function(t,s,mu,r,shift = 0){
  P = -exp(-s*t)+1-(s*exp(-r*mu)/(mu-s))*(exp(t*(mu-s))-1)-shift
  return(P)
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
