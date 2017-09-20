llik = function(pars,tree,sum=TRUE){
  b = c(pars[1],(pars[1]-pars[2])/pars[3],pars[2])
  t = tree$wt[1:(length(tree$wt)-1)]
  n = tree$n[1:(length(tree$n)-1)]
  E = tree$E
  sigma = n*(b[1]-b[2]*n + b[3]) #n-dimentional
  rho = pmax(b[1]*E-b[2]*n*E+b[3]*(1-E),0)
  if(sum) l = -(sum(-sigma*t+log(rho)))#+length(t)*log(sum(t)))
  else l = sigma*t-log(rho)
  if(min(b)<0){l = Inf}
  return(l)
}
