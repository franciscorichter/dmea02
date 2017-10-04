rec.tree <- function(tree,pars,model='dd',seed=0, adjustment=FALSE){
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
  et = NULL # event type
  h = 1
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
            rprob[h] = dexp(x = t.spe, rate = (s+nm*mu))*(lambda[1]/(s+nm*mu))
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
          rprob[h] = dexp(x = mint, rate = (s+nm*mu))*(mu0/(s+nm*mu0))
          et[h] = 'extinction'
          probs = probs[-extinctedone]
          sprob[h] = prod(probs)*truncdist::dtrunc(mint,'exp',a=0,b=(e.lims[extinctedone]-cbt),rate=mu)*(1-integrate(sampprob,lower = 0, upper = mint,s=s,mu=mu,r=ct-cbt)$value)#/integrate(sampprob,lower = 0, upper = ct-cbt,s=s,mu=mu,r=ct-cbt)$value)
          ms = ms[-pickone]
          cwt = cwt - t.ext
          cbt = cbt + t.ext
          N = N-1
          #print(paste('ext',sprob[h]))
          #print(paste('probs',prod(probs),'trunc',dtrunc(mint,'exp',a=0,b=(e.lims[extinctedone]-cbt),rate=mu),'1-int',(1-integrate(sampprob,lower = 0, upper = mint,s=s,mu=mu,r=ct-cbt)$value/integrate(sampprob,lower = 0, upper = ct-cbt,s=s,mu=mu,r=ct-cbt)$value)))
          #print(paste('INTEGRATETRUNCA',integrate(dtrunc,lower = 0, upper = (e.lims[extinctedone]-cbt),spec='exp',a=0,b=(e.lims[extinctedone]-cbt),rate=mu)$value))
          #print(paste('b',(e.lims[extinctedone]-cbt)))
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
         sprob[h] = 1 - integrate(sampprob,lower = 0, upper = cwt,s=s,mu=mu,r=ct-cbt)$value#/integrate(sampprob,lower = 0, upper = ct-cbt,s=s,mu=mu,r=ct-cbt)$value
         #print(paste('nothing',sprob[h]))
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
  #tree$weight = prod(rprob)/prod(sprob)
  logweight = log(rprob)-log(sprob)
  tree$weight = sum(logweight)#+length(tree$wt)*log(sum(tree$wt)) #logweight
  tree$dif1 = dif1
  tree$dif2 = dif2
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


sampprob <- function(t,s,mu,r){
  term1 = s*(1-exp(-mu*(r-t)))
  c = exp(-(s/mu)*exp(-mu*r))
  f = term1*exp(-s*t)*c^{-exp(mu*t)}*c
  return(f)
}
