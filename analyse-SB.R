library(ggplot2)
library("scales")

problems<-c("sphere4angle-bigWFP")  


DEVICE=0  # 0: RWIN, 1: JPEG, 2: PDF
devname = c("RWIN","JPEG", "PDF") 

for(problem in problems){
  load(paste("Results/",problem,".RData",sep=""))
  lbig <- length(bigSORC)
  err = rep(NA,lbig)
  for (k in 1:lbig) {
    mySORC<-bigSORC[[k]]$SORC
    gamScheme <- mySORC$defGamma$gammaScheme
    if (mySORC$defGamma$active==F) gamScheme="constant"
    mtitle<-paste(mySORC$name,"problem, run ",k,
                  "\n fCalc=",mySORC$fCalc,", noiseVar=",mySORC$noiseVar,", TR=",mySORC$doTRIKE,
                  "\n gamma=",mySORC$defGamma$gammaInit,", gammaScheme =",gamScheme)
    error<-abs(mySORC$optim-mySORC$DF$bestFeas)
    if (DEVICE == 0) {
      par(cex.axis=1.0, cex.lab=1.0, cex.main=1.0)
    } else if (DEVICE == 1) {
      jpeg(paste("../Results/",problem,".jpeg",sep=""),width=700,height=700,units = "px")
      par(cex.axis=1.5, cex.lab=1.6, cex.main=1.0)
    } else if (DEVICE == 2) {
      pdf(paste("../Results/",problem,".pdf",sep=""), paper="special", width=6, height=6)
      par(cex.axis=1.5, cex.lab=1.6, cex.main=0.8)
    }
    par(mar=c(2, 4, 4, 2) + 0.1) 
    par(mfrow=c(2,2)) 
    
    
    plot(c(1:length(error)),abs(error),type="l",
         xlab="", #"iteration",
         ylab="optimization error",
         main=mtitle,
         log="y",ylim=c(1e-5,1e+3))
    
    
    ###plot feasibility rate
    myDF<-mySORC$DF
    numViol<-myDF$numViol
    numViol[numViol>1]<-1
    rate<-c()
    for(i in c(1:length(numViol))){
      #browser()
      rate<-c(rate,(i-sum(numViol[1:i]))/i)
      
    }
    plot(c(1:length(rate)),rate,type="l",
         xlab="", #"iteration",
         ylab="feasibility rate",
         main=mtitle)
    
    par(mar=c(5, 4, 2, 2) + 0.1) 
    
    ###plot objective approximation error
    initSize<-length(which(is.na(mySORC$DF$pred)))
    predF<-c(rep(NA,initSize),mySORC$predF)
    objectiveErr<-abs(predF-mySORC$Fres)
    
    plot(c(1:length(objectiveErr)),objectiveErr,type="l",
         xlab="iteration",
         ylab="approximation error",
         log="y",
         main="")
    
    
    
    ##plot the distance to the optimum over iterations
    plot(c(1:length(myDF$dist)),myDF$dist,type="l",
         xlab="iteration",
         ylab="distance to the optimum",
         log="y", ylim=c(0.0000001,1),
         main="")
    
    if (DEVICE > 0) {
      dev.off()
    }
    par(mar=c(5, 4, 4, 2) + 0.1) # restore default
    
    err[k] =  min(mySORC$DF$bestFeas - mySORC$optim,na.rm=T)
    cat(sprintf("%s, run %d: final error = %9.2e\n",
                mySORC$name,k,err[k]))
  } # for (k)
  if (lbig>1) cat(sprintf("%s:   min error=%9.2e, mean final error = %9.2e, median final error = %9.2e, max final error = %9.2e\n",
                          mySORC$name,min(abs(err)),mean(abs(err)),median(abs(err)),max(abs(err))))

  myFrame<-data.frame(x=c(1:length(mySORC$bCritic)),y=mySORC$bCritic,type="bCritic")
  myFrame<-rbind(myFrame,data.frame(x=c(1:length(mySORC$bCritic)),y=myDF$dist,type="dist"))
  p<-ggplot(data=myFrame,aes(x=x,y=y,color=type))
  p<-p+geom_line(size=1.2)
  p<-p+scale_y_log10(breaks=c(1e-7,1e-6,1e-5,1e-4,1e-3,1e-2,1e-1,1e+0,1e+1,1e+2,1e+3,1e+4,1e+5,1e+6))
  plot(p)
} # for (problem)