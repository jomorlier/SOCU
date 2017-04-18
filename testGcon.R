TRL<-c(-1,-1)
TRU<-c(+1,+1)  
a1=20*pi/180      # gradient angle of 1st constraint
a2=40*pi/180      # gradient angle of 2nd constraint
gcon=NULL
gcon=cbind(gcon,c(cos(a1),sin(a1)))
gcon=cbind(gcon,c(cos(a2),sin(a2)))

x1<-seq(from=TRL[1],to=TRU[1],length=100)
x2<-seq(from=TRL[2],to=TRU[2],length=100)
xMat<-expand.grid(x=x1,y=x2)

f.value <- function(x) {
  #   if(!is.data.frame(x)) {
  #     newdata.num <- as.numeric(x)
  #     x <- data.frame(t(newdata.num))
  #   } 
  
  F=1
  for (con in c(1,2)) {
    F=F*pnorm(10*gcon[,con]%*%x)
  }
  return(F)
} 

SURFACE_PLOT=F
if (SURFACE_PLOT) {
  fnMat <- sapply(1:nrow(xMat),function(i)f.value(as.vector(t(xMat[i,]))))
  
  FIP<-ggplot(data=data.frame(xMat,z=fnMat),aes(x=x,y=y,z=z))
  FIP<-FIP+geom_tile(aes(fill=z))+stat_contour(color="black")
  FIP<-FIP+ scale_fill_gradientn(name="F(con)",colours=brewer.pal(6,"YlOrRd"))
  FIP<-FIP+coord_cartesian(xlim=c(TRL[1],TRU[1]),ylim=c(TRL[2],TRU[2]))
  
  plot(FIP)
}

delta=1
t=seq(-100,100)/100*delta

angle0=seq(20,180-20,20)
angle=angle0*pi/180
cgrad=c(1,0)
midp=c(0,0)
par(mfrow=c(1,1))
palette(rainbow(12))
for (i in 1:length(angle)) {
  ai = angle[i]
  Amat = matrix(data=c(cos(ai),sin(ai),-sin(ai),cos(ai)),nrow=2,ncol=2)
  cigrad = Amat %*% cgrad
  xMat2=data.frame(x=cigrad[1]*t+midp[1],y=cigrad[2]*t+midp[2])

  F2 <- sapply(1:nrow(xMat2),function(i)f.value(as.vector(t(xMat2[i,]))))

  if (i==1)  {
    plot(t,F2,type="l",col=i)
  } else {
    lines(t,F2,type="l",col=i)
  }
  
}
legend("topleft",legend=angle0,lty=0,col=1:length(angle),text.col=1:length(angle)) 
