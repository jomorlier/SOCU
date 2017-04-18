# --- analyse-B.R ------------
# To run this script, the loaded .RData files have to be run with EIMOD_DEBUG set to TRUE
# in mainSORC.R (only then mySORC$DF2 will exist)
#
problems<-c("G06-14b-WK")  #"sphere4-13-WK" "G06-2-WK"
#problems<-c("sphere4-11-WK")  #"sphere4-3-WK" "G06-2-WK"
#problems<-c(paste0("sphere4-",c(4:5,9),"-WK")) 
#problems<-c(paste0("G06-",c("10","11"),"-WK")) 

DEVICE=0  # 0: RWIN, 1: JPEG, 2: PDF
devname = c("RWIN","JPEG", "PDF") 

for(problem in problems){
  load(paste("../Results/",problem,".RData",sep=""))  # bigSORC
  lbig <- length(bigSORC)
  err = rep(NA,lbig)
  for (k in 1:lbig) {
    mySORC<-bigSORC[[k]]
    mtitle<-paste(mySORC$name,"problem, run ",k,
                  "\n fCalc=",mySORC$fCalc,", noiseVar=",mySORC$noiseVar,", TR=",mySORC$doTRIKE,
                  "\n gamma=",mySORC$defGamma$gammaInit,", gammaScheme =",mySORC$defGamma$gammaScheme)
    error<-abs(mySORC$optim-mySORC$DF$bestFeas)
    iter<-c(1:length(error))
    if (DEVICE == 0) {
      par(cex.axis=1.0, cex.lab=1.0, cex.main=1.0)
    } else if (DEVICE == 1) {
      jpeg(paste("../Results/B-",problem,".jpeg",sep=""),width=700,height=700,units = "px")
      par(cex.axis=1.5, cex.lab=1.6, cex.main=1.0)
    } else if (DEVICE == 2) {
      pdf(paste("../Results/B-",problem,".pdf",sep=""), paper="special", width=6, height=6)
      par(cex.axis=1.5, cex.lab=1.6, cex.main=0.8)
    }
    par(mar=c(2, 4, 4, 2) + 0.1) 
    par(mfrow=c(2,2)) 
    
    if (is.null(mySORC$DF2)) stop(sprintf("No element DF2 in list mySORC for problem %s",problem))
    plot(c(1:length(error)),error,type="l",
         xlab="", #"iteration"
         ylab="optimization error",
         main=mtitle,
         log="y",ylim=c(1e-5,1e+3))
    
    
    ###plot EImod for solu and infill point based on each iteration's Kriging models
    DF2<-mySORC$DF2
    initSize<-length(which(is.na(mySORC$DF$pred)))
    data<-c(DF2$Eim_Solu,DF2$EimInfil)
    plot(iter,c(rep(NA,initSize),DF2$Eim_Solu),type="l",lwd=2,
         xlab="", #"iteration",
         ylab="EImod",
         log="y",ylim=c(min(data),max(data)),
         main=mtitle,col="red")
    lines(iter,c(rep(NA,initSize),DF2$EimInfil),col="blue",lwd=2)
    legend("topleft",legend=c("solu","infill"),lwd=c(2,2),col=c("red","blue")); 
    
    par(mar=c(5, 4, 2, 2) + 0.1) 
    
    ###plot gamma
    plot(iter,c(rep(NA,initSize),mySORC$DF2$gamma),type="l",
         xlab="iteration",
         ylab="mean gamma",
         ylim=c(0,max(mySORC$DF2$gamma)),
         main="")
    lines(iter,rep(2,length(error)),lty=2,col="blue")
    
    
    ###plot F (feas-prob) for solu and infill point based on each iteration's Kriging models
    plot(iter,c(rep(NA,initSize),DF2$F_Solu),type="l", lwd=2,lty=2,
         xlab="iteration",
         ylab="prob_feasible F",
         ylim=c(0,1),
         main="",col="red")
    lines(iter,c(rep(NA,initSize),DF2$FInfil),col="blue",lwd=2,lty=2)
    legend("topleft",legend=c("solu","infill"),lwd=c(2,2),lty=c(2,2),col=c("red","blue")); 
    
    
    if (DEVICE > 0) {
      dev.off()
    }
    par(mar=c(5, 4, 4, 2) + 0.1) # restore default
    
    err[k] =  min(mySORC$DF$bestFeas - mySORC$optim,na.rm=T)
    cat(sprintf("%s, run %d: final error = %9.2e\n",
                mySORC$name,k,err[k]))
  } # for (k)
  if (lbig>1) cat(sprintf("%s:   mean final error = %9.2e\n",
                          mySORC$name,mean(err)))
} # for (problem)