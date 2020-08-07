
hexd <- data.frame(x = 1 + c(rep(-sqrt(3)/2, 2), 0, rep(sqrt(3)/2,2), 0), y = 1 + c(0.5, -0.5, -1, -0.5, 0.5, 1))
hexd <- rbind(hexd, hexd[1, ])

set.seed=123
x.val=seq(0.4,1.6,length.out=20)
y.val=c(1.17, 1.48, 0.7, 0.65, 1.36, 1.08, 0.94, 0.98, 0.99, 0.92,0.92, 1.03, 1.3, 1.02, 1.18, 1.08, 1.49, 0.86, 0.62, 1.45)
# y.val=runif(length(x.val),0.6,1.6)

cols=colorRampPalette(c("forestgreen","grey","white"))(length(y.val))
# png(filename="./hex/NADA2_hex.png",width=3.5,height=4,units="in",res=250,type="windows",bg=NA)
par(family="serif",oma=c(1,1,1,1),mar=c(0.1,0.1,0.1,0.1),xpd=F)

plot(hexd,axes=F,ann=F,type="n")
polygon(hexd,border="grey",lwd=4,col="white")
lines(x.val,y.val,lty=2,col="grey",lwd=1.5)
points(x.val,y.val,pch=21,bg=cols,cex=1.25,col=ifelse(y.val<1,"red","grey"),lwd=0.75)
lines(c(0.3,1.7),c(1,1),lty=2,col="Red",lwd=1)
text(1,0.6,"NADA2",col="red",pos=1,cex=1.5)
dev.off()
