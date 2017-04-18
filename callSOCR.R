#'
#'
#'Samineh Bagheri
#'11.03.2016
#'Cologne University of Applied Sciences
#'
#'

source("sourceSOCR.R")
dimension<-20
#define your problem or call one of the predefined problems
problem<-"G08"
testName<-"test"
nrun=1
source(paste("testProblem/",problem,".R",sep=""))

#initialize the SORC
saveName<-paste(problem,"-",testName,sep="")
startSeed<-10
bigSORC<-list()
for(mySeed in c(startSeed:(startSeed+nrun-1))){
  cat(testName,"configuration is running for the seed",mySeed,"\n")
  mySORC<-initializeSORC(fn=fn,dLower=-1,dUpper=1,name=saveName,
                         lower=lower,upper=upper,dimension=dimension
                         ,tol=1e-7
                         ,nConstraints=m,kType="UK"
                         ,initSize=5*dimension,initMethod="LHS"
                         ,vis=F
                         ,visZoom=F
                         ,solu=solu
                         ,budget=200
                         ,seed=mySeed
                         ,fCalc="beta"    #alpha, beta, mix 
                         ,wtType="random" #only relevant for fCalc=mix
                         ,noiseVar = 0.00
                         ,handleCrash=T
                         ,PLUGINTYPE="fixed" #"fixed", "adaptive", original
                         ,maxTime=1
                         ,defEPS=list(active=F)
                         ,repairInfeas = F
                         ,constantSigma=NULL
                         ,useEI=T
                         ,verbose=0
                         ,useSBRBF=list(active=T,kernel="Cubic",param=2,ptail=TRUE,squares=TRUE))  
  
  mySORC<-mainSORC(mySORC)  
  bigSORC[[mySeed-startSeed+1]]<-mySORC
}

save(bigSORC,file=paste("Results/",saveName,".RData",sep=""))
