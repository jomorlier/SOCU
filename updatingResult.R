
setPluginValue<-function(SORC,bIndex){

    switch(SORC$PLUGINTYPE,
           #SB:only for dignostic issues and tests on G06,
           #the plugin value is not appropriate for all problems
           "original"={    cplugin<-SORC$Fres[bIndex]},
           "fixed"={cplugin<- 1e+5},
           #"fixed"={cplugin<- max(-6000,SORC$DF$bestFeas[nrow(SORC$DF)])}, 
           #"fixed"={cplugin<-max(SORC$DF$bestFeas[-which(is.infinite(SORC$DF$bestFeas))])},
           "adaptive-old"={SORC$DF$bestFeas[nrow(SORC$DF)]+0.5*abs(SORC$DF$bestFeas[nrow(SORC$DF)])},
           "adaptive"={cplugin<-SORC$DF$bestFeas[nrow(SORC$DF)]+SORC$bCritic[length(SORC$bCritic)]}
    )

  
  return(cplugin)
}


updatingResult<-function(SORC,con.model,objModel){
  myDF<-SORC$DF
  defEPS<-SORC$defEPS
  defGamma<-SORC$defGamma
  
  
  
  
  
  
  
  
  
  
  
  
###################################
# STEP 4: Update the results      #
###################################
  if((SORC$verbose>3) && ((SORC$iter-SORC$initSize)%%SORC$printIter==0))print(paste("log(EI(solu))=",round(log10(-EImod(x=SORC$solu,con.model=con.model,objModel=objModel,SORC=SORC)),digits=2),
            "log(EI-SA)=",round(log10(-EImod(x=SORC$infill,con.model=con.model,objModel=objModel,SORC=SORC)),digits=2))       )
SORC$EISOLU<-c(SORC$EISOLU,-EImod(x=SORC$solu,con.model=con.model,objModel=objModel,SORC=SORC))
SORC$EISA<-c(SORC$EISA,-EImod(x=SORC$infill,con.model=con.model,objModel=objModel,SORC=SORC))



SORC$EI<--EImod(x=SORC$infill,con.model=con.model,objModel=objModel,SORC=SORC)
##prediction of kriging objective on the infill point
infillPred<-predict(object = objModel, newdata = t(SORC$infill), type = SORC$kType, 
                    checkNames = FALSE)

##prediction of kriging objective on the optimum
optimPred<-predict(object = objModel, newdata = t(SORC$solu), type = SORC$kType, 
                   checkNames = FALSE)

#mean of the objective model at the infill point
SORC$predF<-rbind(SORC$predF,infillPred$mean)
#sd of the objective model on the infill point
SORC$sd.i<-rbind(SORC$sd.i,infillPred$sd)

#mean of the objective model at the optimum
SORC$predSolu<-rbind(SORC$predSolu,optimPred$mean)
#sd of the objective model on the optimum
SORC$sd.solu<-rbind(SORC$sd.solu,optimPred$sd) 


predg<-c()
predgsd<-c()
for(con in c(1:SORC$nConstraints)){
  model<-con.model[[con]]
  newdata <- t(as.numeric(SORC$infill))
  infillGPred<-predict(object = model, newdata = newdata, type = SORC$kType, 
                       checkNames = FALSE)
  predg<-c(predg,infillGPred$mean)
  predgsd<-c(predgsd,infillGPred$sd)
}
SORC$predC<-rbind(SORC$predC,predg)
SORC$sd.con<-rbind(SORC$sd.con,predgsd)
SORC$pop<-rbind(SORC$pop,SORC$infill)

####################################
# Calculating bCritic              #
####################################
# Violation<-apply(SORC$Gres,2,FUN=function(x){
#   ret<-quantile(abs(x))[2]
#   return(ret)
# })
#Violation<-abs(SORC$Gres[nrow(SORC$Gres),])
#browser()
#cat("violation:",Violation,"\n")
#activeConstraints<-which(Violation<1e-1)
# if(length(activeConstraints)==0)
#   print("no active constraint is yet detected")
# else
#  cat("active constraint(s):",activeConstraints,"\n")
# 

if(SORC$PLUGINTYPE=="adaptive" ){
 # browser()
  #calculating gradient of EI 
  Grad<-grad(func=myEI,x=SORC$infill,method="simple",method.args=list(eps=1e-6),
             objModel=objModel,
             plugin=SORC$plugin[length(SORC$plugin)],
             type=SORC$kType)
  #browser()
  GradSize<-sqrt(sum(Grad^2))
  probGrad<-1
  gGrad<-c()
  conGrad<-NULL
  gTrueGrad<-c()
  conTrueGrad<-NULL
  conInfill<-NULL
  sigma<-Inf

  for(con in c(1:SORC$nConstraints)){
   # browser()
    conFunc <- function(x) {
      model<-con.model[[con]]
      x<-t(as.matrix(x))
      if(!SORC$useSBRBF$active) colnames(x) = colnames(model@X)
      predx <- predict(object = model, newdata = x, type = SORC$kType, checkNames = FALSE)
      kriging.mean <- predx$mean
      if(is.null(SORC$constantSigma)){
         kriging.sd <- predx$sd
      }else{
        kriging.sd<-SORC$constantSigma    # /WK/ changed from 1 to SORC$constantSigma
      }
      prob<- kriging.mean/kriging.sd      # this is the function -z(x)
      sigma<<-(min(sigma,kriging.sd))
      return(prob)
    }
    probGrad<-grad(func=conFunc,x=SORC$infill,method="simple",method.args=list(eps=1e-6))
    # gGrad<-c(gGrad,as.numeric(Grad%*%probGrad))
    conGrad<-rbind(conGrad,probGrad)
    gGrad<-c(gGrad,sqrt(sum(probGrad^2)))
    
    predx <- predict(object=con.model[[con]], newdata=t(as.matrix(SORC$infill)), type=SORC$kType, checkNames=FALSE)
    conInfill <- c(conInfill,predx$mean)
    
    # only debug /WK/
    probTrueGrad<-grad(func=function(x){(SORC$fn(x))[con+1]},x=SORC$infill,method="simple",method.args=list(eps=1e-6))
    conTrueGrad<-rbind(conTrueGrad,probTrueGrad)
    gTrueGrad<-c(gTrueGrad,sqrt(sum(probTrueGrad^2)))
    #
    GradMod<-grad(func=function(x){-EImod(x,con.model,objModel,SORC)},x=SORC$infill,method="simple",method.args=list(eps=1e-6))
    mgrad <- GradMod/sqrt(sum(GradMod^2))
    GradF<-grad(func=function(x){CalculateF(x,con.model,SORC)$res},x=SORC$infill,method="simple",method.args=list(eps=1e-6))
    fgrad <- GradF/sqrt(sum(GradF^2))
    
    
  } # for (con)
  
  minAbs <- max(min(abs(conInfill)),1e-4)     # NEW! /WK/ 
  index <- which(abs(conInfill)<10*minAbs)    # this is an estimate of all active constraints
  index <- 1:SORC$nConstraints   # this is a simpler heuristic: take all constraints
#  index<-c(1:length(activeConstraints)) #SB: A very simpel definition for active constraint is used in line 86 before bCritic is calculated
  cgrad <- 0
  for (con in index) {
    cgrad <- cgrad + conGrad[con,]/gGrad[con]    
  }
  cgrad=cgrad/sqrt(sum(cgrad^2))    # direction of slowest descent is assumed to be the direction of all
                                    # active constraint gradients. This is true at least for G06.

  # /WK/ ?? Are the next four lines needed anywhere ??
#   predx <- predict(object = model, newdata = t(as.matrix(SORC$infill)), type = SORC$kType, checkNames = FALSE)
#   kriging.mean <- predx$mean
#   kriging.sd <- predx$sd  
#   sigma<-kriging.sd
#   
  cosa<-(conGrad%*%Grad)/(gGrad*GradSize)
  cosa2<-(conGrad%*%cgrad)/(gGrad*1)     # NEW! /WK/
  cosb<-(Grad%*%cgrad)/(GradSize*1)     # NEW! /WK/
  if (GradSize==0) cosb<-1   # avoid cosb=NaN in this pathological case

#   if (!is.null(SORC$DF)) {
#     if (SORC$DF$bestFeas[nrow(SORC$DF)]< -6960 & sqrt(sum((SORC$infill-SORC$solu)^2))<1e-4) {
#       browser()
#     }
#   }
  
  index<-which(cosa>0)      # NEW! /WK/ cosa>0 (not <0): constraint-gradient in the half-space of grad(EI)
                            # Why '>0'? - Because conGrad has the function -z(x). For z(x) it would be '<0'  
  #index<-which(cosa2>0)    # Alternative (not recommended): include all constraints whose 
                            # constraint-gradient is in in half-space of cgrad 
  #print(index)
  cosa<-abs(cosa)
  cosa2<-abs(cosa2)
  cosb<-abs(cosb)   # just for safety, should be positive anyway
  
OLD_BCRITIC=FALSE

if(length(index)!=0) { #} && gGrad !=0){    # && GradSize!=0 
  bCritics<-c()
  for(ind in index){                                                        
    bCritics<-c(bCritics,(GradSize*cosb*sqrt(2*pi)*0.5/(cosa2[ind]*gGrad[ind]))) 
  }
  if (OLD_BCRITIC) {
    print("!!! OLD_BCRITIC !!!")
    cosa<-(conGrad%*%Grad)/(gGrad*GradSize)
    index<-which(cosa<0)       
    cosa<-abs(cosa)
    for(ind in index){                                                    # old version: SORC$bcfac=10
      bCritics<-c(bCritics,(GradSize*sqrt(2*pi)*0.5/(cosa[ind]*(gGrad[ind]/SORC$bcfac)))) 
    }
  }
  bCritic<-max(bCritics)[1] 
  usedInd<-which(bCritics==bCritic)[1]
  
}else{
  bCritic<-0 
  usedInd<-1
}
if (is.nan(bCritic) || is.infinite(bCritic)) {
 # browser()    # just for safety, should normally not happen
  #SB: 11.01.2017 It is observed that this case happens for G04 Problem with SOCU120 setting
  bCritic<-0
}

########################################
# /WK/ Two improvements seem necessary:
#   a) handle the case when index is empty: should bCritic be zero or should we take one of the cosa>=0 inidices?
#   b) bCritic should be calculated for all constraints which are in index and the max value should be taken
########################################
SORC$bCritic<-c(SORC$bCritic,abs(bCritic))
print(sprintf("(a,k,kTrue)=(%10.3f,%10.3f,%10.3f), sigma is %9.5f, bCritic is %9.5f"
              ,GradSize,gGrad[usedInd],gTrueGrad[usedInd],sigma,bCritic)) 

}
#else{ # i.e. if SORC$PLUGINTYPE != "adaptive"
 # SORC$bCritic<-c(SORC$bCritic,0)
  #SORC$bCritic<-0.1*abs(SORC$DF$bestFeas[nrow(SORC$DF)])
 # SORC$bCritic<-10       # /WK/ ???? why 10 ????
  
  ### and why is this else-branch needed at all? SORC$bCritic will by ONLY used if PLUGINTYPE=="adaptive"
#}
#browser()



###################################
# STEP 5: Replacement add-on      #
###################################
if(SORC$doREPLACE){
  closeInd<-which(dist(rbind(SORC$infill,SORC$fpop))<1e-5)
  # if(length(closeInd)>0)browser()
  closeInd<-closeInd[closeInd<=nrow(SORC$fpop)]
  if(length(closeInd)>0){
    print("deleting a point")
    delInd<-closeInd
    SORC$fpop<-SORC$fpop[-delInd,]
    SORC$fGres<-as.matrix(SORC$fGres[-delInd,])
    SORC$fFres<-SORC$fFres[-delInd]
  }  
}

SORC$fpop<-rbind(SORC$fpop,SORC$infill)
# cat(nrow(SORC$fpop),nrow(SORC$fGres),length(SORC$fFres),"\n")



#evaluate the new solution
eval<-SORC$fn(SORC$infill)
obj<-eval[1]
gs<-eval[-1]

maxViol<-gs
maxViol[maxViol<SORC$tol]<-0
maxViol<-max(maxViol)
numViol<-sum(maxViol>0)

switch(SORC$wtType,
       "iLin"={SORC$wt=(SORC$iter-SORC$initSize)/(SORC$budget-SORC$initSize)},
       "dLin"={SORC$wt=1-((SORC$iter-SORC$initSize)/(SORC$budget-SORC$initSize))},
       "adaptive"={SORC$wt<-1-(length(which(myDF$numViol==0))/SORC$iter)},
       "random"={rand<-runif(1);if(rand>0.4){SORC$wt<-1}else{SORC$wt<-0}},
       stop("invalid wt scheme"))

# browser()
bfeas<-myDF$bestFeas[nrow(myDF)]

if(numViol==0 ){ 
  if(is.na(bfeas)==T ||( obj<myDF$bestFeas[nrow(myDF)]))
    bfeas<-obj
}
if(length(SORC$gamma)==1){
  gamma<-c(myDF$gamma,SORC$gamma)
}else{
  gamma<-rbind(myDF$gamma,SORC$gamma)
}
#browser()
myDF<-data.frame(iter=c(1:nrow(SORC$pop)),
                 # tau=c(myDF$tau,round(SORC$tau,2)),
                 gamma=gamma,
                 wt=c(myDF$wt,round(SORC$wt,2)),
                 obj=c(myDF$obj,obj),
                 pred=c(rep(NA,SORC$initSize),SORC$predF),
                 frate=c(myDF$frate,length(which(myDF$numViol==0))/length(myDF$numViol)),
                 plugin=SORC$plugin,
                 bestFeas=c(myDF$bestFeas,bfeas),
                 type=c(as.character(myDF$type),SORC$type),
                 maxViol=c(myDF$maxViol,maxViol),
                 numViol=c(myDF$numViol,numViol),
                 solu=SORC$pop,
                 dist=c(myDF$dist,as.numeric(dist(rbind(SORC$solu,SORC$infill)))),
                 predSolu=c(rep(NA,SORC$initSize),SORC$predSolu),
                 sd.i=c(rep(NA,SORC$initSize),SORC$sd.i),
                 sd.solu=c(rep(NA,SORC$initSize),SORC$sd.solu)
)
SORC$DF<-myDF
SORC$Fres<-c(SORC$Fres,obj)
SORC$Gres<-rbind(SORC$Gres,gs)

AnalyseConstraints(SORC$Gres)
#browser()

###################################
# STEP 6: Adaptive mechanism      #
#         (epsilon or gamma)      #
###################################
#infeasgs<-gs
violated<-which(gs>0)
nviolated<-which(gs<=0) 
#infeasgs[violated]<-1
#infeasgs[nviolated]<-0

SORC$infeasCount[violated]<- SORC$infeasCount[violated]+1 # infeasgs
SORC$infeasCount[nviolated]<-0

SORC$feasCount[violated]<-0
SORC$feasCount[nviolated]<-SORC$feasCount[nviolated]+1


SORC$infeasCountG[violated]<- SORC$infeasCountG[violated]+1 # infeasgs
SORC$infeasCountG[nviolated]<-0

SORC$feasCountG[violated]<-0
SORC$feasCountG[nviolated]<-SORC$feasCountG[nviolated]+1


if(defEPS$active==F){
  SORC$EPS<-rep(0,SORC$nConstraints)
}else{
  
  #changing the EPS
  feasTrigger <- which(SORC$feasCount==defEPS$EPSthreshold)
  infeasTrigger <- which(SORC$infeasCount==defEPS$EPSthreshold)
  triggered<-union(infeasTrigger,feasTrigger)
  SORC$EPS[triggered]  <-SORC$EPS[triggered]  + defEPS$kappa*gs[triggered]
  
  SORC$EPS <- pmin(SORC$EPS,1)
  SORC$EPS <- pmax(SORC$EPS,0)
  
  SORC$feasCount[feasTrigger]<-0
  SORC$infeasCount[infeasTrigger]<-0
  
}

#if(defGamma$active==F){
 # SORC$gamma<-rep(2,nConstraints) #SB: we have to talk about this logic, if defGamma$active==F is offline how to calculate beta(x)
  #WK: I think that defGamma$active should be always TRUE
#}else{
  feasTriggerG<-c()
  infeasTriggerG<-c()
  ##SB:finding feasTrigger and infeasTrigger for gamma scheme (later we can only use one united parameter for epsilon and gamma scheme both)
  
  feasTriggerG <- which(SORC$feasCountG==defGamma$gammaThreshold)
  infeasTriggerG <- which(SORC$infeasCountG==defGamma$gammaThreshold)
  # browser()
  switch(defGamma$gammaScheme,
         "constant"={SORC$gamma<-rep(defGamma$gammaInit,SORC$nConstraints)},
         "linear"={gamma<-defGamma$gammaInit-((SORC$iter+1-SORC$initSize)/(SORC$budget-SORC$initSize))*(defGamma$gammaInit-1);SORC$gamma<-rep(gamma,SORC$nConstraints)},
         "adaptive"={
           changeFactor<-rep(0.5,SORC$nConstraints);SD<-rep(1,SORC$nConstraints);
           triggered<-union(infeasTriggerG,feasTriggerG)
           for(con in triggered){
             model<-con.model[[con]]
             newdata<-SORC$infill
             newdata <- data.frame(t(as.numeric(newdata)))
             #browser()
             if(!SORC$useSBRBF$active) colnames(newdata) = colnames(model@X)
             predx <- predict(object = model, newdata = newdata, type = SORC$kType, 
                              checkNames = FALSE)
             SD[con] <- predx$sd
           };
           changeFactor[triggered]<- pnorm((-gs[triggered]*defGamma$kappa)/SD[triggered]);
           # if(length(triggered)>0)browser()
           SORC$gamma<-2*SORC$gamma*changeFactor;
         }
  ) # switch
  
  #       if(length(triggered)>0){print("changing gamma") 
  #         browser()}
  
  SORC$gamma <- pmin(SORC$gamma,2)
  SORC$gamma <- pmax(SORC$gamma,1)  # /WK/ was 1 before, why 0.75?
  
  SORC$feasCountG[feasTriggerG]<-0
  SORC$infeasCountG[infeasTriggerG]<-0
  
#}

  if((SORC$verbose>3) && ((SORC$iter-SORC$initSize)%%SORC$printIter==0))cat(sprintf("SORC$gamma : %s\n",paste(SORC$gamma,collapse=", ")))

#only debugging
  if((SORC$verbose>3) && ((SORC$iter-SORC$initSize)%%SORC$printIter==0)){
    cat(sprintf("SORC$feasCount   : %s\n",paste(SORC$feasCountG,collapse=" ")))
    cat(sprintf("SORC$infeasCount : %s\n",paste(SORC$infeasCountG,collapse=" ")))  
    cat(sprintf("prod(gamma) : %s\n",paste(prod(SORC$gamma),collapse=" ")))
    
  }
  

  if(defEPS$active){
    if((SORC$verbose>3) && ((SORC$iter-SORC$initSize)%%SORC$printIter==0)){
    cat(sprintf("log10(SORC$EPS) : %s\n",paste(log10(SORC$EPS),collapse=", ")))
    cat(sprintf("SORC$feasCount   : %s\n",paste(SORC$feasCount,collapse=" ")))
    cat(sprintf("SORC$infeasCount : %s\n",paste(SORC$infeasCount,collapse=" ")))  
    }
    # only debug info:
    nC <- SORC$nConstraints
    newrow <- SORC$FGEPSres[1,]   # to copy the names
    newrow[1] <- obj
    newrow[2:(nC+1)] <- gs
    newrow[(nC+2):(2*nC+1)] <- SORC$EPS
    SORC$FGEPSres <- rbind(SORC$FGEPSres,newrow)
    rownames(SORC$FGEPSres) <-  SORC$DF$iter
  }
  



###################################
# STEP 7: more bookkeeping        #
###################################
SORC$fFres<-c(SORC$fFres,obj)
SORC$fGres<-rbind(SORC$fGres,gs) 
SORC$iter<-nrow(SORC$pop)
#SORC$tau<--(tauInit/SORC$budget)*SORC$iter+tauInit
#SORC$EISA<-EISA
#SORC$EISOLU<-EISOLU

if (SORC$EIMOD_DEBUG) {
  F_2<-CalculateF(SORC$solu,con.model,SORC)$res
  F_3<-CalculateF(SORC$infill,con.model,SORC)$res
  C_1<-CalculateCon(SORC$infill,con.model,SORC)$conVec
  CTr<-SORC$fn(SORC$infill)[-1]
  C_S<-CalculateCon(SORC$solu,con.model,SORC)$conVec
  # browser()
  myDF2 <- data.frame(iter=nrow(SORC$pop)
                      ,gamma=mean(SORC$gamma)     # NOTE: with 'adaptive', different values are possible in SORC$gamma
                      ,EimXinit=SORC$eiList$ei1
                      ,Eim_Solu=SORC$eiList$ei2
                      ,EimInfil=SORC$eiList$ei3
                      ,F_Solu=F_2
                      ,FInfil=F_3
                      ,xInit.1=SORC$xInit[1]
                      ,xInit.2=SORC$xInit[2]
                      ,CInfil.1=C_1[1]   # constraint 1 Kriging model value at infill point 
                      ,CInfil.2=C_1[2]
                      ,CT_inf.1=CTr[1]   # constraint 1 true value at infill point
                      ,CT_inf.2=CTr[2]
                      ,CSolu.1=C_S[1]    # constraint 1 Kriging model value at true solution 
                      ,CSolu.2=C_S[2]     
  )
  SORC$DF2<-rbind(SORC$DF2,myDF2)
}
if((SORC$verbose>-1) && ((SORC$iter-SORC$initSize)%%SORC$printIter==0)) print(myDF[nrow(myDF),c(1,c((2+SORC$nConstraint):(9+SORC$nConstraint)))])
# print(c(myDF[nrow(myDF),1],prod(myDF[nrow(myDF),c(2:(SORC$nConstraints+2))]),myDF[nrow(myDF),c(SORC$nConstraints+2:6)]))
#if(SORC$vis && SORC$dimension==2){
  # print(paste("EI(solu)=",round(-EImod(x=SORC$solu,con.model=con.model,objModel=objModel,SORC=SORC),digits=2)
  #             ,"max-EI=",round(SORC$EImax,digits=2)))
#}

#deprecated
#if(SORC$doREPLACE!=T && SORC$handleCrash!=T){
 # SORC$fpop<-SORC$pop
 # SORC$fGres<-SORC$Gres
 # SORC$fFres<-SORC$Fres
#}
save(SORC,file=paste("../Results/",SORC$name,SORC$SORCseed,".RData",sep=""))

###################################
# STEP 8: optional visualization  #
###################################
visT<-proc.time()
if(SORC$vis && SORC$dimension==1){
  visualizeOne(objModel=objModel,fn=SORC$fn,SORC=SORC,con.model=con.model)
}
#browser()
if(SORC$vis && SORC$dimension!=1){
  if (SORC$defMiniTR$active) {
    TRL<-SORC$optimres1$par-SORC$lcurrent/2
    TRL<-pmax(TRL,SORC$dLower)
    TRU<-SORC$optimres1$par+SORC$lcurrent/2
    TRU<-pmin(TRU,SORC$dUpper)
  } else {
    TRL<-SORC$TRL
    TRU<-SORC$TRU
  }

  
  SORC<-visualizeTwo(objModel=objModel,fn=SORC$fn,SORC=SORC,con.model=con.model,TRL=TRL,TRU=TRU)
  if((SORC$verbose>0) && ((SORC$iter-SORC$initSize)%%SORC$printIter==0))print(paste("visualisation took:",round((proc.time()-visT)[3],digits=3)))
  
}





###################################
#         Find the current        #
#          Best Solution          #
###################################
#If Feasible Solu are exitisng in the popolation:
# The best solu is the best feasible solution with the smallest objective value

#If no feasible solution is provided:
# The best solution is the one with the lowest max. Viol

#Is there any feasible solu?
if(any(myDF$numViol==0)){
  fIndex<-myDF$iter[which(myDF$numViol==0)]
  fDF<-myDF[fIndex,]
  bIndex<-fDF$iter[min(which(fDF$obj==min(fDF$obj)))]
}else{
  #No feasible solution available
  bIndex<-which(myDF$maxViol==min(myDF$maxViol)) 
  # plugin<-10^10
}
if(SORC$iter>SORC$initSize) if(bIndex==SORC$bIndex){
  SORC$AI<-0
}else{
  SORC$AI<-SORC$Fres[SORC$bIndex]-SORC$Fres[bIndex]
}
bIndex <- bIndex[1]       # /WK/ to avoid a crash later
SORC$bIndex<-bIndex
bestI<-SORC$pop[bIndex,]
# browser()
cplugin<-setPluginValue(SORC,bIndex)
if(is.infinite(cplugin))cplugin<-10^10    # /WK/ shouldn't this be larger, e.g. 10^10 as in the comment above?
SORC$cplugin<-cplugin
SORC$bestSolu<-bestI
SORC$plugin<-c(SORC$plugin,SORC$cplugin)


return(SORC)


}