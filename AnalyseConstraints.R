AnalyseConstraints<-function(Gres){
  gGauss<-apply(Gres,2,FUN=function(x){
    mean=mean(x)
    sd=sd(x)
    #probability of each constraints to be an active constraint
    pActiveCon<-pnorm(-sd,mean=mean,sd=sd,lower.tail = FALSE)
    
    return(c(mean=mean,sd=sd,pActiveCon=pActiveCon))
  })
  colnames(gGauss)<-sprintf(paste("g",c(1:ncol(gGauss)),sep=""))
  if( ((nrow(Gres)-5*ncol(Gres))%%10==0))print(gGauss)
}