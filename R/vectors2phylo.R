vectors2phylo <- function(list){
  t=list$wt
  n=list$n
  E=list$E
  S=list$S
  ct=sum(t)
  newick = paste(sl[1],";",sep="")
  N=1
  identf = data.frame(Spec="aa",Time=0) # Labels of species
  for (i in 1:(length(t)-1)){
    # speciation
    sumt = sum(t[1:i])
    if( is.null(S)){
      BD = sample(1:N,1)
      species = as.character(identf[BD,1])
    }else{
      species = S[i]
    }
    if (E[i] == 1){
      ind = regexpr(species,newick)[1]-1
      atm = sumt-identf[which(identf[,1]==species),2]
      newick = paste(substr(newick,1,ind),"(",substr(newick,ind+1,ind+4),",",sl[i+1],"):",as.character(atm),substring(newick,ind+5),sep="")
      identf = rbind(identf,data.frame(Spec=substr(sl[i+1],1,2),Time=sumt))
      identf[identf$Spec == species,2] = sumt
      N = N+1
    }
    # extinction
    if (E[i]==0){
      ind = regexpr(species,newick)[1] + 2
      atm = sumt-identf[which(identf[,1]==species),2]
      identf = identf[!identf$Spec==species,]
      newick = paste(substr(newick,1,ind),as.character(atm),substring(newick,ind+2),sep="")
      N=N-1
    }
  }
  newick = compphyl(newi=newick,identf=identf,ct=ct)
  newick = read.tree(text=newick)
  return(newick)
}


