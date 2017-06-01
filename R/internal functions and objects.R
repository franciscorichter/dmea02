sl = paste(letters[1],letters,":0",sep="")
for (i in 2:26){
  ll = paste(letters[i],letters,":0",sep="")
  sl = c(sl,ll)
}
SL = paste(LETTERS[1],LETTERS,":0",sep="")
for (i in 2:26){
  ll = paste(LETTERS[i],LETTERS,":0",sep="")
  SL = c(SL,ll)
}
S1 = paste(LETTERS[1],0:9,":0",sep="")
for (i in 2:26){
  ll = paste(LETTERS[i],1:10,":0",sep="")
  SL = c(SL,ll)
}
s1 = paste(letters[1],0:9,":0",sep="")
for (i in 2:26){
  ll = paste(letters[i],1:10,":0",sep="")
  SL = c(SL,ll)
}
sl=c(sl,SL,S1,s1)
compphyl <- function(newi,identf,ct){
  #set to extant species to the present time
  identf[,1] = as.character(identf[,1])
  identf[,2] = ct-identf[,2]
  for(i in 1:length(identf[,1])){
    ind = regexpr(identf[i,1],newi)[1] + 2
    newi = paste(substr(newi,1,ind),as.character(identf[i,2]),substring(newi,ind+2),sep="")
  }
  return(newi)
}

