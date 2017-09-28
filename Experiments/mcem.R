s = sim.tree(seed=104,mu0=0.1)
cmle = mle.tree(s$tree)
ddmle = as.numeric(DDD::dd_ML(branching.times(s$phylo.extant),initparsopt = c(0.8,0.4,40)))

nem = 10
MCEM = vector(mode='list',length=nem)
for(i in 1:nem){
  MCEM[i] = MCEM.phylo(tree=s$tree.extant,init_par = c(0.5,0.05,90),impsam = TRUE,parallel = TRUE,ips = 100,cutTime = 10)
}

sss = c('s1','s2','s3','s4','s5','s6','s7','s8','s9','s10')
tmcem = rbind(cbind(1:100,as.data.frame(MCEM[1])[-101,],sss[1]))
colnames(tmcem) = c('it','lambda','mu','K','s')

for(i in 2:10){
  t2mcem = rbind(cbind(1:100,as.data.frame(MCEM[i])[-101,],sss[i]))
  colnames(t2mcem) = c('it','lambda','mu','K','s')
  tmcem = rbind(tmcem,t2mcem)
}

ggplot(data=tmcem, aes(x=it,y=lambda,col=s))+geom_line()+geom_hline(yintercept=cmle[1],col='blue')#+geom_hline(yintercept=ddmle[1],col='blue')
ggplot(data=tmcem, aes(x=it,y=mu,col=s))+geom_line()+geom_hline(yintercept=cmle[2],col='blue')#+geom_hline(yintercept=ddmle[2],col='blue')
ggplot(data=tmcem, aes(x=it,y=K,col=s))+geom_line()+geom_hline(yintercept=cmle[3],col='red')+geom_hline(yintercept=ddmle[3],col='blue')


dim(tmcem)

nem = 10
MCEM2 = vector(mode='list',length=nem)
for(i in 1:nem){
  MCEM2[i] = MCEM.phylo(tree=s$tree.extant,init_par = c(2.5,0.9,90),impsam = TRUE,parallel = TRUE,ips = 100,maxNtree = 10)
}


sss = c('s11','s12','s13','s14','s15','s16','s17','s18','s19','s20')
tmcem2 = rbind(cbind(1:100,as.data.frame(MCEM2[1])[-101,],sss[1]))
colnames(tmcem2) = c('it','lambda','mu','K','s')
tmcem = rbind(tmcem,tmcem2)
for(i in 2:10){
  t2mcem = rbind(cbind(1:100,as.data.frame(MCEM2[i])[-101,],sss[i]))
  colnames(t2mcem) = c('it','lambda','mu','K','s')
  tmcem = rbind(tmcem,t2mcem)
}

ggplot(data=tmcem, aes(x=it,y=lambda,col=s))+geom_line()+geom_hline(yintercept=cmle[1],col='red')+geom_hline(yintercept=ddmle[1],col='blue')
ggplot(data=tmcem, aes(x=it,y=mu,col=s))+geom_line()+geom_hline(yintercept=cmle[2],col='red')+geom_hline(yintercept=ddmle[2],col='blue')
ggplot(data=tmcem, aes(x=it,y=K,col=s))+geom_line()+geom_hline(yintercept=cmle[3],col='red')+geom_hline(yintercept=ddmle[3],col='blue')



