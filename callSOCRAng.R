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
nrun=1
KS<-c(40,30,20,10,5,4,3,2,1)
for(k in KS){
  problem<-"sphere4angle"
  angle<-pi/k
  d<-dimension<-2
  testName<-paste("comp",which(k==KS),sep="")
  source(paste("testProblem/",problem,".R",sep=""))
  #initialize the SORC
  saveName<-paste(problem,"-",testName,sep="")
  startSeed<-3103543
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
                           ,budget=50
                           ,seed=mySeed
                           ,fCalc="beta"  #alpha, beta, mix                         only relevant for mix 
                           ,wtType="random"
                           ,noiseVar = 0.001
                           ,handleCrash=T
                           ,PLUGINTYPE="adaptive" #"fixed", "adaptive", original
                           ,maxTime=10
                           ,repairInfeas = T
                           ,defGamma=list(gammaScheme="adaptive",gammaThreshold=2)
                           ,constantSigma=NULL
                           ,ri=list(eps3=0.0,q=1,mmax=1000,eps2=0.0,eps1=1e-3) )  
    
    mySORC<-mainSORC(mySORC)  
    bigSORC[[mySeed-startSeed+1]]<-mySORC
  }
  
  save(bigSORC,file=paste("Results/",saveName,".RData",sep=""))
  
}




