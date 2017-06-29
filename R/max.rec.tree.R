rec.tree <- function(wt,pars,model='dd'){
  ct = sum(wt)
  sumt = 0
  if(model == 'dd'){
    lambda0 = pars[1]
    mu0 = pars[2]
    K = pars[3]
  }
  #n
  i = 1
  ms = NULL # missing species, for now we just add time. When we consider topology we do it with species as well
  sumt = 0
  while(i < length(wt)){
    ctm = wt[i]
    if(model == "dd"){  # diversity-dependence model
      lambda = max(0,lambda0 - (lambda0-mu0)*N/K)
      mu = mu0
      lambda = rep(lambda,N)
      mu = rep(mu,N)
    }
    if(model == 'cr'){ # constant-rate model
      lambda = rep(lambda0,N)
      mu = rep(mu0,N)
    }

    s = sum(lambda)+sum(mu)
    if (s == 0){break}
    tm = rexp(1,s)  # waiting time of iteration i
    sumt = sumt + tm
    if(tm < ctm){
      if(length(ms) == 0){ # speciation
        ms = c(ms,sumt)
      }
      if(max(lambda)==0){#extinction
        if(length(ms)==0){
          print('Not extinction and speciation posible!!')
        }
        pickone = sample(1:length(ms),1)
        t_spe = ms[pickone]
        t_ext = sumt
        tree = update.tree(tree,t_spe,t_ext)
        ms[pickone] = NULL
      }
      else{ # extinction (in the case of max dimensions)
        # existe una probabilidad de que alguna de las missing species se extinga,
        # dado que todas deben extinguirse.
        # tabien hay probabilidad de que algunas de las que existen se hayan especiado,
        # dado que las que se especiaron se extingan.
        #E ...
        nm = 0
        #update.tree()

      }
    }
  }

}
