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
problem<-"G06" # "sphere4" "G06" "G01" G04 G05
testName<-"14n-WK"
nrun=10
source(paste("testProblem/",problem,".R",sep=""))

#lower=solu-0.01 # 0.1 # 0.01  # **** debug only ****
#upper=solu+0.01 # 0.1 # 0.01

#initialize the SORC
saveName<-paste(problem,"-",testName,sep="")
startSeed<-120 # 120 121
bigSORC<-list()
for(mySeed in c(startSeed:(startSeed+nrun-1))){
  cat(testName,"configuration is running for the seed",mySeed,"\n")
  SORC<-initializeSORC(fn=fn,dLower=-1,dUpper=1,name=saveName,
                       lower=lower,upper=upper,dimension=dimension,
                       tol=1e-6, # 0 1e-6 1e-5
                       nConstraints=m,kType="UK",
                       initSize=max(3*dimension,6),
                       initMethod="LHS", 
                       vis=F, visZoom=T,
                       solu=solu,
                       budget=100, # 40 50 100
                       seed=mySeed,
                       doTRIKE=F,
                       #miniTR=T,
                       doNUGGET=F,
                       doConNUGGET=T,
                       doREPLACE=F,
                       noiseVar=0.05,     #NULL 0.1 0.05 0.05*70
                       conNoiseVar=NULL, # NULL 0.0005
                       wtType="iLin", #iLin , dLin, adaptive
                       fCalc="beta"  #alpha, beta, mix
                       ,smoothSA=F     # control parameter for GenSA (default: F)
                       ,maxTime=10     # 10, time in sec for GenSA (default: 30)
                       ,defEPS = list(active=F,
                                      kappa=1.0,
                                      EPSthreshold=3)
                       ,defGamma = list(active=F,
                                        gammaInit=2,
                                        gammaScheme= "adaptive" # constant linear adaptive
                       )
                       ,defMiniTR = list(active=F,
                                         scheme="adaptive", # "adaptive" "constant"
                                         mTRwidth=0.25)
                       #,gammaInit = 2
                       #,gammaScheme = "constant" # "constant", "linear"
                       #,nu=1 # 0.25 0.5 1.0
                       #,EPSCrate=5e-6 # 0 0.01
                       ,repairInfeas = T
                       #,ri=list(eps3=1e-10,eps2=1e-6,q=3,mmax=1000) 
                       ,NUmech=F
                       ,constantSigma=NULL # fixed value (0.2) or NULL (take sigma from Kriging model)
                       ,TRUEPLUGIN=T
                       ,PLUGINTYPE = "fixed" # fixed or adaptive
                       ,bcfac=10             # for OLD_BCRITIC only
  )  
  
  lst<-mainSORC(SORC)  
  SORC<-lst$SORC
  bigSORC[[mySeed-startSeed+1]]<-SORC
}

save(bigSORC,file=paste("../Results/",saveName,".RData",sep=""))

