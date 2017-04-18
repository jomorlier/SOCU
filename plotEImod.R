# How the maximum of EImod shifts in the 1D-case, when the plugin leads to a near-zero b.
# We assume that EImod can be approximated as a linear function near the boundary:
#       EImod(x) = pmax(a*x+b,0)
#
DEVICE=2  # 0: RWIN, 1: JPEG, 2: PDF
setDevice <- function (i) {
  par(mfrow=c(1,1))
  if (DEVICE == 0) {
    par(cex.axis=1.0, cex.lab=1.0, cex.main=1.0)
  } else if (DEVICE == 1) {
    jpeg(sprintf("../Results/EImod1D-%02d.jpg",i),width=700,height=700,units = "px")
    par(cex.axis=1.5, cex.lab=1.6, cex.main=1.0)
  } else if (DEVICE == 2) {
    pdf(sprintf("../Results/EImod1D-%02d.pdf",i), paper="special", width=6, height=6)
    par(cex.axis=2.0, cex.lab=1.9, cex.main=1.7, mar=c(5,5,4,2)+0.1)
  }
}

x <- 2.5*seq(-0.2,0.2,0.001)
gamma <- 2
sigma <- 1/5
FF <- pmin(gamma*pnorm(x/sigma),1)
a <- -20;
bSafe <- abs(a)*sqrt(2*pi)*sigma/2;     # for all b<bSafe the maximum of EImod shifts 
                                    # from 0 to negative values (!)
bvec <-  sort(c(0.02,0.2,2,20,200,bSafe),decr=T)
for (i in 1:length(bvec)) {
  b <- bvec[i]
  setDevice(i)
  EI <- pmax(a*x+b,0)
  #EI <- a*x+b
  max1=max(FF,EI*FF)
  #tit=paste0(sprintf("a=%5.1f,\n b=%5.2f, gamma=%5.2f",a,b,gamma))
  tit=paste0(sprintf("a=%6.1f, b=%5.2f",a,b))
  plot(x,EI*FF,type="l",lwd=3,ylim=c(max(-0.5,-b),1.1*max1),main=tit,ylab="EImod")
  par(cex.axis=1.5, cex.lab=1.6, cex.main=1.7)
  lines(x,FF,col="red")
  lines(x,EI,col="green")
  lines(x,2.5*max1*(pnorm(1000*x)-0.5),col="blue") # the boundary
  if (DEVICE > 0) {
    dev.off()
  }
}
par(mar=c(5, 4, 4, 2) + 0.1) # restore default

