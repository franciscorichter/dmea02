llik.st = function(pars, setoftrees, impsam = F){
  m = length(setoftrees)
  l = vector(mode = 'numeric',length = m)
  for(i in 1:m){
    s = setoftrees[[i]]
    if(impsam){
      weight = s$weight
      l[i] = llik(pars=pars,tree=s)*weight
    }
    else{
      l[i] = llik(pars=pars,tree=s)
    }
  }
  L = sum(l)
  return(L)
}
