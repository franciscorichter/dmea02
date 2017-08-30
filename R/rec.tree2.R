rec.tree2 <- function(tree,pars,model='dd',seed=0){
  if(seed>0) set.seed(seed)
  wt = tree$wt
  ct = sum(wt)
  lambda0=pars[1]
  mu0=pars[2]
  K=pars[3]
  limit = lambda0*K/(lambda0-mu0)
  rs = dim = length(wt)
  if (limit < (dim-1)){
    # print('parameters does not make sense, observed tree implies negative rates')
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
    # print(paste('at branching time',cbt,'an observed species arises, resulting on',N,'current species, happening before the current waiting time',cwt))
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
        print(maxprob)
        t.spe = rexp(1,s)  # waiting time of iteration i
        probs = dexp(x = t.spe,rate = s,log = TRUE)
        if(probs>0){
          print('probs dando logs positivos??')
          print(t.spe)
          print(s)
        }
      }
      if(nm > 0){ # if there are missing species
        t.ext = NULL
        for(j in 1:nm){
          t.ext[j] = rtrunc(1,'exp',a=ms[j],b=e.lims[j],rate=mu)
        }
        extinctedone = which(t.ext == min(t.ext))
        t.ext = min(t.ext)-cbt
        #probe = dtrunc(x=t.ext,'exp',a=nm[extinctedone],b=e.lims[extinctedone],rate=mu,log=TRUE)
        if(probe>0) print('probe dando logs positivos??')
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
            #prob = prob + probs + thre + dists(wt = cwt, lambda=lambda[-1],mu = rep(mu0,nm),b = e.lims,log = TRUE)
            prob = c(prob, probs + log(thre) + dists(wt = cwt, lambda=lambda[-1],mu = rep(mu0,nm),b = e.lims,log = TRUE)+gprob)
            gprob=0
            nm = nm + 1 # number of missing species
            if((N + rs-1) < limit){
              e.lims = c(e.lims,ct)
            }
            else{
              e.lims = c(e.lims,sum(wt[1:(dim - (N+rs-floor(limit)+1))])) # cual es la interpretacion de esta formula?
            }
            N = N+1
            #print(paste('at branching time',cbt+t.spe,'a missing species arises, resulting on',N,'current species, happening before the current waiting time',cwt))
          }
          else{ # nothing happens, move to next (gost) point
            #prob = prob + probs + (1-thre) + dists(wt = cwt, lambda=lambda[-1],mu = rep(mu0,nm),b = e.lims,log = TRUE)
            #prob = c(prob,probs + (1-thre) + dists(wt = cwt, lambda=lambda[-1],mu = rep(mu0,nm),b = e.lims,log = TRUE))
            gprob = gprob + probs + log(1-thre) + dists(wt = cwt, lambda=lambda[-1],mu = rep(mu0,nm),b = e.lims,log = TRUE)
            #print(paste('there was a non-acepted species at branching time',cbt+t.spe,'a non included species apears'))
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
          #prob = prob + probe  + dists(wt = cwt,lambda = lambda,mu = rep(mu0,nm),b = e.lims,log = TRUE)
          prob = c(prob,probe  + dists(wt = cwt,lambda = lambda,mu = rep(mu0,nm),b = e.lims,log = TRUE)+gprob)
          gprob = 0
          #print(paste('at branching time',cbt,'a missing species got exincted, resulting on',N,'current species, happening before the current waiting time',cwt))
        }
      }
      else{
        #prob = prob + dists(wt = cwt, lambda=lambda,mu = rep(mu0,nm),b = e.lims,log = TRUE)
        prob = c(prob,dists(wt = cwt, lambda=lambda,mu = rep(mu0,nm),b = e.lims,log = TRUE)+gprob)
        gprob = 0
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

