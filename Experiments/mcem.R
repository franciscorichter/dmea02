s = sim.tree(seed=104,mu0=0.4)
cmle = mle.tree(s$tree)
ddmle = as.numeric(DDD::dd_ML(branching.times(s$phylo.extant),initparsopt = c(0.8,0.4,40)))

mcem1=MCEM.phylo(tree=s$tree.extant,init_par = c(0.5,0.09,90),impsam = TRUE,parallel = TRUE)
mcem2=MCEM.phylo(tree=s$tree.extant,init_par = c(0.5,0.09,90),impsam = TRUE,parallel = TRUE)
mcem3=MCEM.phylo(tree=s$tree.extant,init_par = c(0.5,0.09,90),impsam = TRUE,parallel = TRUE)
mcem4=MCEM.phylo(tree=s$tree.extant,init_par = c(0.5,0.09,90),impsam = TRUE,parallel = TRUE)
mcem5=MCEM.phylo(tree=s$tree.extant,init_par = c(0.5,0.09,90),impsam = TRUE,parallel = TRUE)



mcemT = rbind(cbind(1:45,mcem1$pars[1:45,],'s1'),cbind(1:45,mcem2$pars[1:45,],'s2'),cbind(1:45,mcem3$pars[1:45,],'s3'),cbind(1:45,mcem4$pars[1:45,],'s4'),cbind(1:45,mcem5$pars[1:45,],'s5'),cbind(1:45,mcem6$pars[1:45,],'s6'),cbind(1:45,mcem7$pars[1:45,],'s7'),cbind(1:45,mcem8$pars[1:45,],'s8'),cbind(1:45,mcem9$pars[1:45,],'s9'))
colnames(mcemT) = c('it','lambda','mu','K','s')
mcemT=as.data.frame(mcemT)
qplot(x=as.numeric(mcemT$it),y=as.numeric(as.character(mcemT$lambda)),col=mcemT$s)+geom_line()
qplot(x=as.numeric(mcemT$it),y=as.numeric(as.character(mcemT$mu)),col=mcemT$s)+geom_line()
qplot(x=as.numeric(mcemT$it),y=as.numeric(as.character(mcemT$K)),col=mcemT$s)+geom_line()


mcem10=MCEM.phylo(tree=s$tree.extant,init_par = c(0.5,0.07,90),impsam = TRUE,parallel = TRUE)

mcem11=MCEM.phylo(tree=s$tree.extant,init_par = c(1.5,0.07,90),impsam = TRUE,parallel = TRUE)
mcem12=MCEM.phylo(tree=s$tree.extant,init_par = c(1.5,0.07,90),impsam = TRUE,parallel = TRUE)
mcem13=MCEM.phylo(tree=s$tree.extant,init_par = c(1.5,0.07,90),impsam = TRUE,parallel = TRUE)
mcem14=MCEM.phylo(tree=s$tree.extant,init_par = c(1.5,0.07,90),impsam = TRUE,parallel = TRUE)
mcem15=MCEM.phylo(tree=s$tree.extant,init_par = c(1.5,0.07,90),impsam = TRUE,parallel = TRUE)


cmle = mle.tree(s$tree)
ggplot(x=as.numeric(mcemT$it),y=as.numeric(as.character(mcemT$K)),col=mcemT$s)+geom_vline(intercept=cmle[3])
ggplot(data=mcemT,aes(x=as.numeric(it),y=as.numeric(as.character(lambda)),col=s))+geom_line()+geom_hline(yintercept=cmle[1])
ggplot(data=mcemT,aes(x=as.numeric(it),y=as.numeric(as.character(mu)),col=s))+geom_line()+geom_hline(yintercept=cmle[2])
ggplot(data=mcemT,aes(x=as.numeric(it),y=as.numeric(as.character(K)),col=s))+geom_line()+geom_hline(yintercept=cmle[3])

ddmle = dd_ML(branching.times(s$phylo.extant),initparsopt = c(0.8,0.4,40),idparsopt = 1:3,cond = 0,res = 200)

ggplot(data=mcemT,aes(x=as.numeric(it),y=as.numeric(as.character(lambda)),col=s))+geom_line()+geom_hline(yintercept=as.numeric(ddmle[1]))
ggplot(data=mcemT,aes(x=as.numeric(it),y=as.numeric(as.character(mu)),col=s))+geom_line()+geom_hline(yintercept=as.numeric(ddmle[2]))
ggplot(data=mcemT,aes(x=as.numeric(it),y=as.numeric(as.character(K)),col=s))+geom_line()+geom_hline(yintercept=as.numeric(ddmle[3]))
