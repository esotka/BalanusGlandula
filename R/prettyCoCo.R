### pretty co-co plot
rm(list=ls())
bs.coef.all <- read.csv(file="output/May19.snps-CoCo.cline-parameters.csv")

bs.coef.all2 <- bs.coef.all[!bs.coef.all$locus%in%c("percCOI.C","percEF1.C"),]
bs.coef.all2 <- bs.coef.all2[bs.coef.all2$p.val1<0.05,] # significant broken stick regressions (why p-value 1 and not p-value 2? )
#af2 <- af[,colnames(af)%in%bs.coef.all2$locus] # significant broken stick regressions (why p-value 1 and not p-value 2? )

require(ks)
pdf("output/prettyCoCo.pdf")
fhat <- with(bs.coef.all2,kde(x=cbind(-mid,log(abs(slope))),H=Hscv(x=cbind(-mid,log(abs(slope))))))
plot(fhat,drawpoints=F,cont=c(99,95,75),ylab="log(slope)",xlim=c(-55,-35),xaxt="none",xlab="midpoints (ºN)")
axis(side=1,labels = sort(seq(35,55,5),decreasing = T),at = seq(-55,-35,5))
with(bs.coef.all2,points(x=-mid,y=log(abs(slope)),pch=20,cex=0.5))
with(bs.coef.all[bs.coef.all$locus=="percCOI.C",],points(x=-mid,y=log(abs(slope)),pch=21,cex=1,col="red",lwd=2))
with(bs.coef.all[bs.coef.all$locus=="percEF1.C",],points(x=-mid,y=log(abs(slope)),pch=20,cex=1,col="red",lwd=2))
dev.off()

 
