#'
#'
#'Samineh Bagheri
#'11.03.2016
#'Cologne University of Applied Sciences
#'
#'

#call all the necessary functions and libraries
#If you recieve by executing the following command
#Please consider installing the missing Packages
source("sourceSOCR.R")

#define your problem or call one of the predefined problems
problem<-"G06"
d<-dimension<-2
testName<-"test"
nrun=1
source(paste("testProblem/",problem,".R",sep=""))

#initialize the SORC
saveName<-paste(problem,"-",testName,sep="")
startSeed<-1
bigSORC<-list()
for(mySeed in c(startSeed:(startSeed+nrun-1))){
  cat(testName,"configuration is running for the seed",mySeed,"\n")
  mySORC<-initializeSORC(fn=fn,dLower=-1,dUpper=1,name=saveName,
                         lower=lower,upper=upper,dimension=dimension
                         ,tol=1e-5
                         ,nConstraints=m,kType="UK"
                         ,initSize=5*dimension,initMethod="LHS"
                         ,vis=F
                         ,visZoom=F
                         ,solu=solu
                         ,budget=100
                         ,seed=mySeed
                         ,fCalc="beta"  #alpha, beta, mix                         
                         ,noiseVar = 0.01
                         ,handleCrash=T
                         ,PLUGINTYPE="adaptive" #"fixed", "adaptive", original
                         ,maxTime=10
                         ,repairInfeas = F
                         ,constantSigma=NULL )  
  
  mySORC<-mainSORC(mySORC)  
  bigSORC[[mySeed-startSeed+1]]<-mySORC
}

save(bigSORC,file=paste("Results/",saveName,".RData",sep=""))
