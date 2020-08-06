### do a co-co analysis with John Wares' SNP dataset
rm(list=ls())
library(adegenet) # read.genepop
library(scales) # alpha colors

meta <- read.delim('data/meta.txt')
dat <- read.genepop('data/May19.snps.gen')
locusnames <- colnames(dat@tab)
freq.ind <- dat@tab
af <- c()
for(i in 1:dim(freq.ind)[2])
{
  tmp.count <- tapply(freq.ind[,i],dat@pop,sum,na.rm=T)
  tmp.n <- tapply(freq.ind[,i],dat@pop,length) - tapply(is.na(freq.ind[,i]),dat@pop,sum) ### to account for missing samples
  tmp.freq <- as.numeric(tmp.count/(2*tmp.n))
  af <- cbind(af,tmp.freq)
  colnames(af)[i] <- colnames(freq.ind)[i]
  
}
write.csv(af,'output/allelefreq.csv',row.names=names(tmp.n))

lats.ordered <- meta$lats[match(names(tmp.n),meta$popID)]
pdf('output/May19.snps-CoCo.pdf')
plot(af[,1]~lats.ordered,type="n",ylab="n",xlab="p")
apply(af[,1:dim(af)[2]],2,function(x){points(x~lats.ordered,type="l",col=alpha("grey",0.5))})


### fit broken stick model

source('R/broken.R')
bs.coef.all <- NULL
for (j in 1:dim(af)[2])
{
  plt <- data.frame(lats.ordered,af[,j])
  plt <- plt[complete.cases(plt),]
  bs.fit <- bs(plt[,1],plt[,2])
  flt.fit<- flt(plt[,1],plt[,2])
  sig.fit <- sig(plt[,1],plt[,2])
  cf <- bs.fit$par[1:4]
  clnpar <- extract.bs(bs.fit$par)
  lrt1 <- c(G1=2*abs(bs.fit$value-flt.fit$value),p.val1=(1-pchisq(2*abs(bs.fit$value-flt.fit$value),3)))
  lrt2 <- c(G2=2*abs(bs.fit$value-sig.fit$value),p.val2=(1-pchisq(2*abs(bs.fit$value-sig.fit$value),1)))
  bs.coef.all <- rbind(bs.coef.all,c(locus=names(plt)[2],clnpar,lrt1,lrt2,cf,sig.fit$par))
}

bs.coef.all <- as.data.frame(bs.coef.all)
names(bs.coef.all)[8:14] <- c("X1","Y1","X2","Y2","c","w","sig.sd")
for(i in 2:14) {bs.coef.all[,i] <- as.numeric(as.character(bs.coef.all[,i]))}
bs.coef.all$locus <- colnames(dat@tab)

write.table(file="output/May19.snps-CoCo.cline-parameters.csv",sep=',',row.names=F,bs.coef.all)


bs.coef.all2 <- bs.coef.all[bs.coef.all$p.val1<0.05,] # significant broken stick regressions (why p-value 1 and not p-value 2? )
af2 <- af[,colnames(af)%in%bs.coef.all2$locus] # significant broken stick regressions (why p-value 1 and not p-value 2? )

require(ks)
fhat <- with(bs.coef.all2,kde(x=cbind(mid,log(abs(slope))),H=Hscv(x=cbind(mid,log(abs(slope))))))
plot(fhat,drawpoints=F,cont=c(99,95,75),ylab="log(abs(slope))",main="All broken sticks significant")#,ylim=c(-4,9),xlim=c(14,19))
plot(fhat,drawpoints=F,cont=c(99,95,75),ylab="log(abs(slope))",main="All broken sticks significant")#,ylim=c(-4,9),xlim=c(14,19))

with(bs.coef.all2,text(x=mid,y=log(abs(slope)),locus,
                      col=1+(p.val1>0.05),cex=.5))


# plot clines - significant ones only 

for (k in 1:dim(af2)[2])
{
  plt <- cbind(lats.ordered,p=af2[,k])
  plt <- plt[complete.cases(plt),]
  plot(plt,type="b",lwd=2,main=colnames(af2)[k])
  bs.row <- which(bs.coef.all$locus==colnames(af2)[k])
  points(x=c(0,bs.coef.all$X1[bs.row],bs.coef.all$X2[bs.row],46),
         y=c(bs.coef.all$Y1[bs.row],bs.coef.all$Y1[bs.row],bs.coef.all$Y2[bs.row],bs.coef.all$Y2[bs.row]),type="l",col="red")
}




#for (l in 6:dim(kaw.ibd)[2])
#{
#  plt <- kaw.ibd[,c(2,l)]
#  plt <- plt[complete.cases(plt),]
#  
#  plot(plt,type="b",lwd=2,main=names(plt)[2])
#  bs.row <- which(bs.coef.all$locus==names(plt)[2])
#  points(x=c(0,bs.coef.all$X1[bs.row],bs.coef.all$X2[bs.row],46),
#         y=c(bs.coef.all$Y1[bs.row],bs.coef.all$Y1[bs.row],bs.coef.all$Y2[bs.row],bs.coef.all$Y2[bs.row]),type="l",col="red")

#  points(x=plt[,1],y=((1+tanh(2*(plt[,1]-bs.coef.all$c[bs.row])/bs.coef.all$w[bs.row]))/2),
#         type="l",col="green")

#}
dev.off()