llik.st = function(pars, setoftrees, impsam = F, correction=0){
  m = length(setoftrees)
  l = vector(mode = 'numeric',length = m)
  w = vector(mode = 'numeric',length = m)
  for(i in 1:m){
    s = setoftrees[[i]]
    l[i] = llik(pars=pars,tree=s)+length(s$tree$wt)*correction
    w[i] = s$weight
  }
  if(impsam){
      weight = exp(w-correction*D)/max(exp(w-correction*D))
      l = l*weight
  }
  L = sum(l)
  return(L)
}
