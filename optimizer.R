
runOptimizer<-function(SORC,objModel,con.model){
#  myDF<-SORC$DF
#  defEPS<-SORC$defEPS
#  defGamma<-SORC$defGamma
###################################
# STEP 3: max EImod use           #
#    GenSA or COBYLA              #
###################################
saT<-proc.time()
###################################
# STEP 3a: optional trust region  #
#          and random start step  #
###################################
if(SORC$doTRIKE)SORC<-TRIKE(SORC)
currentU<-as.vector(as.numeric(SORC$TRU))
currentL<-as.vector(as.numeric(SORC$TRL))

xInit<-as.numeric(SORC$bestSolu) 

#deprecated
# print(c(currentL,currentU))
# randNum<-runif(1)
# if(randNum<0){
#   print("random start")
#   xInit<-runif(dimension,min=currentL,max=currentU)
# }else{
#   
#   
#   #calculating the closest solu
#   violatedInd<-which(SORC$DF$maxViol!=0.0 & SORC$DF$maxViol<1e-1 )
#   #close<-which(SORC$DF$maxViol[violatedInd]==min(SORC$DF$maxViol[violatedInd]))
#   close<-which(SORC$Fres[violatedInd]==min(SORC$Fres[violatedInd]))
#   
#   if(length(close)==0)
#     SORC$closestSOLU<-xInit
#   else
#     SORC$closestSOLU<-as.numeric(SORC$pop[violatedInd[close],])
#   
#   
#   if(SORC$closeSOLU){
#     xInit<-SORC$closestSOLU
#   }
# }

EIMOD_DEBUG <- SORC$EIMOD_DEBUG <- F
DOGENOUD=F          # /WK/ experimental
if (EIMOD_DEBUG) {
  # /WK/ --- only debug ---
  printEI <- function(x,txt) {
    Ei = myEI(x,objModel,plugin=SORC$plugin[length(SORC$plugin)],type=SORC$kType) # myEI in feasibilityMeasure.R
    EiMod = -EImod(x,con.model,objModel,SORC)
    conEI<-CalculateF(x,con.model,SORC)$res
    cat(sprintf(paste(SORC$optimizer,": %s = (%5.3f,%5.3f), EImod(%s) = %8.6f ( EI(%s) = %8.2e, F(%s)  = %8.2e) \n"), 
                txt,x[1],x[2], txt, EiMod, txt, Ei, txt, conEI)) 
    return(EiMod)
  }
  ei1<-printEI(xInit,"xInit")
  ei2<-printEI(SORC$solu," solu")
  SORC$eiList <- list(ei1=ei1,ei2=ei2)   # needed by updatingResult
  SORC$xInit <- xInit
  # /WK/ --- end  debug ---
}  

if (!DOGENOUD && SORC$optimizer=="GenSA") {
 # if(nrow(SORC$pop)>=31)browser()
  optimres<-GenSA(par=xInit,
                  fn=function(x,con.model,objModel,SORC){
                    # y<--log10(-EImod(x,con.model,objModel,SORC))
                    y<-EImod(x,con.model,objModel,SORC)
                    return(y)
                  },
                  lower=currentL, 
                  upper=currentU, 
                  control=list(max.time=SORC$maxTime,
                               verbose=F,smooth=SORC$smoothSA),   # /WK/ smooth=F seems very important here!!!
                  con.model=con.model,objModel=objModel,SORC=SORC) 
 # browser()
  
 
}
if(DOGENOUD){
  # /WK/ --- alternative optimizer, experimental ---
  optimres<-genoud(starting.values=xInit,
                   fn=function(x,con.model,objModel,SORC){
                     return(EImod(x,con.model,objModel,SORC))
                   },
                   nvars=length(xInit),
                   Domains=cbind(currentL,currentU),
                   boundary.enforcement=2,    # no child outside Domains   
                   print.level=0,
                   max=F,                     # do minimization
                   con.model=con.model,objModel=objModel,SORC=SORC) 
}

if(SORC$optimizer=="COBYLA"){
  optimres<-cobyla(x0=xInit,
                   fn=function(x,con.model,objModel,SORC){
                     # y<--log10(-EImod(x,con.model,objModel,SORC))
                     y<-EImod(x,con.model,objModel,SORC)
                     return(-y)
                   },
                   lower=currentL, 
                   upper=currentU, 
                   con.model=con.model,objModel=objModel,SORC=SORC)
}


if (SORC$defMiniTR$active) {
  #--- /WK/ optional: a fixed trust region step within the same iteration: ---
  if (SORC$defMiniTR$scheme=="adaptive") {
    SORC$mTRwidth <- SORC$lcurrent/2
  } else {   # "constant"
    SORC$mTRwidth <- SORC$defMiniTR$mTRwidth # 0.1 or 0.25 
  }
  val1 <- optimres$value
  optimres1<-optimres
  SORC$optimres1<-optimres1
  optimres<-GenSA(par=optimres$par,
                  fn=function(x,con.model,objModel,SORC){
                    y<-EImod(x,con.model,objModel,SORC)
                    return(y)
                  },
                  lower=pmax(optimres$par-SORC$mTRwidth,SORC$dLower), 
                  upper=pmin(optimres$par+SORC$mTRwidth,SORC$dUpper), 
                  control=list(max.time=10,verbose=F,smooth=SORC$smoothSA),
                  con.model=con.model,objModel=objModel,SORC=SORC) 
  if (val1 < optimres$value)         stop("Sanity check value in miniTR failed.")
  if (any(optimres$par<SORC$dLower)) stop("Sanity check lower in miniTR failed.")
  if (any(optimres$par>SORC$dUpper)) stop("Sanity check upper in miniTR failed.")
  optimres$counts <- optimres$counts+optimres1$counts  # total counts from two calls to GenSA 
  optimres2 <- optimres
}

#     GRID_DEBUG=FALSE
#     if (GRID_DEBUG) {
#       xMat <- gridSearch(lower=currentL,upper=currentU,control=list(length=10))
#       #xMat <- gridSearch(lower=SORC$solu-0.2,upper=SORC$solu+0.2,control=list(length=10))
#       ei <- -EImod(xMat, con.model, objModel, SORC)
#     }

if (EIMOD_DEBUG) {
  # /WK/ --- only debug ---
  ei3<-printEI(optimres$par,"infil")
  newdata <- data.frame(t(as.numeric(xInit)))
  newdata <- rbind(newdata,SORC$solu)
  newdata <- rbind(newdata,optimres$par)
  predx <- predict(object = objModel, newdata = newdata, type = SORC$kType, 
                   checkNames = FALSE)
  cat(sprintf(paste(SORC$optimizer,": obj(xInit)$mean,sd = %8.6f, %8.6f\n"), predx$mean[1], predx$sd[1]))
  cat(sprintf(paste(SORC$optimizer,": obj( solu)$mean,sd = %8.6f, %8.6f\n"), predx$mean[2], predx$sd[2]))
  cat(sprintf(paste(SORC$optimizer,": obj(infil)$mean,sd = %8.6f, %8.6f\n"), predx$mean[3], predx$sd[3]))
  ##EiMod = -EImod(xInit+0.001,con.model,objModel,SORC)         # a nearby point
  ##cat(sprintf("GenSA: EImod(xInit+0.001) = %8.6f \n",EiMod))
  # /WK/ --- end  debug ---
  #ei0<-printEI(xInit,"xInit")
  #ei7 = myEI(xInit,objModel,plugin=SORC$plugin,type=SORC$kType)
  #print(ei7)
  #if (SORC$iter>=7) browser()
  SORC$eiList = c(SORC$eiList,ei3=ei3)
}

##Analysis 
#if(all(SORC$solu<=currentU) && all(SORC$solu>=currentL))
 # print("TRUE Solution is somewhere in TR")
#else
#  print("TRUE Solution is *NOT* in TR")

#       optimres<-GenSA(par=xInit,
#                       fn=EImod,
#                       lower=currentL,
#                       upper=currentU,
#                       control=list(max.time=30,verbose=F),
#                       con.model=con.model,objModel=objModel,SORC=SORC) 
if(any(optimres$par==currentL) || any(optimres$par==currentU))
{
  print("SA returned a solu on the TR bounds")
  # SORC$lcurrent<-min(SORC$lcurrent*2,4) #temporary
}

#  }else{
#       optimres<-cobyla(x0=as.numeric(SORC$bestSolu),
#                        fn=EImod,lower=rep(SORC$dLower,SORC$dimension),
#                        upper=rep(SORC$dUpper,SORC$dimension),
#                        con.model=con.model,objModel=objModel,SORC=SORC) 
#  }


#     
#     GenSAres<-GenSA(par=NULL,
#                     fn=EImod,lower=rep(SORC$dLower,SORC$dimension),
#                     upper=rep(SORC$dUpper,SORC$dimension),
#                     control=list(max.time=60,verbose=F),
#                     con.model=con.model,objModel=objModel,SORC=SORC)

# browser()
if((SORC$verbose>0) && ((SORC$iter-SORC$initSize)%%SORC$printIter==0))print(paste("optimizing EImod took:",round((proc.time()-saT)[3],digits=3), ", counts: ",optimres$counts))


diff<-SORC$pop-outer(rep(1,nrow(SORC$pop)),optimres$par)
diff<-abs(diff)

if(any(apply(diff,1,FUN=function(x){sum(x)})==0)){
  #  browser()
    optimres$par<-runif(dimension,min=as.numeric(SORC$bestSolu+(2e-5)),max=as.numeric(SORC$bestSolu+(3e-5)))
}

SORC$infill<-optimres$par




return(SORC)

}
