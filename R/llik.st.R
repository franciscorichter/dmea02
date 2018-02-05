llik.st = function(pars, setoftrees, impsam = F, correction=0){
  m = length(setoftrees)
  l = vector(mode = 'numeric',length = m)
  w = vector(mode = 'numeric',length = m)
  #D = vector(mode = 'numeric',length = m)
  for(i in 1:m){
    s = setoftrees[[i]]
    l[i] = llik(pars=pars,tree=s)+length(s$tree$wt)*correction
    if(impsam) w[i] = s$logweight
    #D[i] = length(s$wt)
  }
  if(impsam){
      #weight = exp(w-correction*D)/sum(exp(w-correction*D))
      l = l*(exp(w)/sum(exp(w)))
  }
  L = sum(l)
  return(L)
}
