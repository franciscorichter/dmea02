llik = function(pars,tree,newlik=FALSE,newlik2=FALSE,newlik3=FALSE){
  b = c(pars[1],(pars[1]-pars[2])/pars[3],pars[2])
  ct=sum(tree$wt)
  t = tree$wt[1:(length(tree$wt)-1)]
  lastn=tree$n[length(tree$n)]
  n = tree$n[1:(length(tree$n)-1)]
  E = tree$E
  sigma = n*(b[1]-b[2]*n + b[3]) #n-dimentional
  rho = pmax(b[1]*E-b[2]*n*E+b[3]*(1-E),0)
  #if(sum)
  l = -(sum(-sigma*t+log(rho))+length(t)*log(sum(t)))
  if(newlik){
    l=-(sum(-sigma*t+log(rho)))#-length(t)*log(sum(t)))
  }
  if(newlik2){
    l=-(sum(-sigma*t+log(rho)-log(t))+ length(t)*log(sum(t)))
  }
  if(newlik3){
    rem=c(ct,ct-cumsum(t))
    lastt=rem[length(rem)]
    rem=rem[-length(rem)]
    l=-(sum(sigma*t+log(rho)-pexp(q = rem,rate = sigma,log.p = TRUE))+pexp(q=lastt,rate = lastn*(b[1]-b[2]*lastn + b[3]),lower.tail = TRUE,log.p = TRUE))
  }
  #else l = sigma*t-log(rho)
  if(min(b)<0){l = Inf}
  return(l)
}
