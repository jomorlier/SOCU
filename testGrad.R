# this is just a small test routine to produce line plots for EImod(x) or F(x)

con.model=lst$con.model
objModel=lst$objModel
conGrad<-NULL
gGrad<-NULL
for(con in c(1:SORC$nConstraints)){
  conFunc <- function(x) {
    model<-con.model[[con]]
    x<-t(as.matrix(x))
    colnames(x) = colnames(model@X)
    predx <- predict(object = model, newdata = x, type = SORC$kType, checkNames = FALSE)
    if(is.null(SORC$constantSigma)){
      kriging.sd <- predx$sd
    }else{
      kriging.sd<-SORC$constantSigma    # /WK/ changed from 1 to SORC$constantSigma
    }
    prob <- predx$mean/kriging.sd
    return(prob)
  }
  probGrad<-grad(func=conFunc,x=SORC$infill,method="simple",method.args=list(eps=1e-6))
  gGrad<-c(gGrad,sqrt(sum(probGrad^2)))
  # gGrad<-c(gGrad,as.numeric(Grad%*%probGrad))
  conGrad<-rbind(conGrad,probGrad)
} # for (con)

delta=1e-4
t=seq(-10,10)/10*delta
# ngrad=Grad/GradSize
# ngrad=fgrad
# xMat=data.frame(x=ngrad[1]*t+SORC$solu[1],y=ngrad[2]*t+SORC$solu[2])
# E1<- -EImod(xMat,con.model,objModel,SORC)
# F1 <- CalculateF(xMat,con.model,SORC)$res
# plot(t,F1,type="l")

c1=conGrad[1,]/gGrad[1]
c2=conGrad[2,]/gGrad[2]
#cgrad=-(conGrad[1,] + conGrad[2,])
cgrad=+(c1+c2)
cgrad=cgrad/sqrt(sum(cgrad^2))
print(acos(sum(cgrad*c1))*180/pi)
print(acos(sum(cgrad*c2))*180/pi)
#cgrad=c(fgrad[2],-fgrad[1])
angle=c(0,-10,-5,+5,+10)*pi/180
for (i in 1:length(angle)) {
  ai = angle[i]
  Amat = matrix(data=c(cos(ai),sin(ai),-sin(ai),cos(ai)),nrow=2,ncol=2)
  cigrad = Amat %*% cgrad
  xMat2=data.frame(x=cigrad[1]*t+SORC$solu[1],y=cigrad[2]*t+SORC$solu[2])
  E2<- -EImod(xMat2,con.model,objModel,SORC)
  F2 <- CalculateF(xMat2,con.model,SORC)$res
  if (i==1)  {
    plot(t,F2,type="l",col="red")
  } else {
    lines(t,F2,type="l",col="black")
    
  }
  
}
  
Grad<-grad(func=myEI,x=SORC$infill,method="simple",method.args=list(eps=1e-6),
           objModel=objModel,plugin=SORC$plugin[length(SORC$plugin)],type=SORC$kType)
egrad<-Grad/sqrt(sum(Grad^2))

#print(1/sum(egrad*c2))
#print(1/sum(cgrad*c2))

