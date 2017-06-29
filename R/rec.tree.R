rec.tree <- function(tree,pars,model='dd'){
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
  prob = 0
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
        probs = dexp(x = t.spe,rate = s,log = TRUE)
      }
      if(nm > 0){ # if there are missing species
        t.ext = NULL
        for(j in 1:nm){
          t.ext[j] = rtrunc(1,'exp',a=0,b=(e.lims[j]-cbt),rate=mu)
        }
        extinctedone = which(t.ext == min(t.ext))
        t.ext = min(t.ext)
        probe = dtrunc(x=t.ext,'exp',a=0,b=(e.lims[extinctedone]-cbt),rate=mu,log=TRUE)
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
            prob = prob + probs + thre + dists(wt = cwt, lambda=lambda[-1],mu = rep(mu0,nm),b = e.lims,log = TRUE)
            nm = nm + 1 # number of missing species
            if((N + rs-1) < limit){
              e.lims = c(e.lims,ct)
            }
            else{
              e.lims = c(e.lims,sum(wt[1:(dim + (floor(limit)-(N+rs-1)))]))
            }
            N = N+1
          }
          else{
            prob = prob + probs + (1-thre) + dists(wt = cwt, lambda=lambda[-1],mu = rep(mu0,nm),b = e.lims,log = TRUE)
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
          prob = prob + probe  + dists(wt = cwt,lambda = lambda,mu = rep(mu0,nm),b = e.lims,log = TRUE)
        }
      }
      else{
        prob = prob + dists(wt = cwt, lambda=lambda,mu = rep(mu0,nm),b = e.lims,log = TRUE)
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
    product = prod((exp(-wt*mu)-exp(-wt*b))/(1-exp(-mu*b)))
  }
  prob = exp(-wt*s)*product
  if(log){
    prob = log(prob)
  }
  return(prob)
}
