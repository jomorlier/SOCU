#'
#'
#'Samineh Bagheri
#'mainSORC.R
#'
#'03.03.2016
#'
#'SAMCO Workshop Constraint Group
#'















mainSORC<-function(SORC){
  
  #adapting predictRBF to be compatible for this framework
  #SB: I know this is not the nicest way of switching between km.predict and SBRBFpredict 
  #but for the moment it works and later we should work more effecient with environments
  
  if(SORC$useSBRBF$active){
    source("SBRBF.R")
    predict<-function(object, newdata , type ,checkNames ){
      if(!is.matrix(newdata))
        newdata<-as.matrix(newdata)
      output<-predictSBRBF(newdata=newdata,model=object)
      return(list(mean=output$mu,sd=output$sigma))
      
    }
    assign("predict",predict, envir=globalenv())
  }else{

    assign("predict",DiceKriging::predict, envir=globalenv())
  }
  
  
  maxViol<-apply(SORC$Gres,1,FUN=function(x){
    viol<-x
    viol[viol<0]<-0
    return(max(viol))
  })
  numViol<-apply(SORC$Gres,1,FUN=function(x){
    return(sum(x>0))
  })

  if(length(which(numViol==0))!=0)
    bfeas<-min(SORC$Fres[which(numViol==0)])
  else
    bfeas=Inf
  SORC$plugin<-rep(NA,nrow(SORC$pop))
  SORC$bCritic<-rep(0,nrow(SORC$pop))
  myDF<-data.frame(iter=c(1:nrow(SORC$pop)),
                   #tau=rep(NA,nrow(SORC$pop)),
                   plugin=SORC$plugin,
                   wt=rep(0,nrow(SORC$pop)),
                   obj=SORC$Fres,
                   frate=length(which(numViol==0)/length(numViol)),
                   bestSolu=rep(NA,nrow(SORC$pop)),
                   bestFeas=rep(bfeas,nrow(SORC$pop)),
                   type=as.character(rep("init",nrow(SORC$pop))),
                   maxViol=maxViol,
                   numViol=numViol,
                   solu=SORC$pop,
                   dist=rep(NA,SORC$initSize))
  SORC$DF<-myDF
  #tauInit<-max(median(maxViol),SORC$tol)
 # browser()
  #SORC$tau<-tauInit
  SORC$nuCoef<-1
  SORC$TRU<-SORC$dUpper
  SORC$TRL<-SORC$dLower
  SORC$AI<-0
  SORC$EI<-Inf
  #SORC$lcurrent<-4
  SORC$bIndex<-Inf
  SORC$TRU<-rep(SORC$dUpper,SORC$dimension)
  SORC$TRL<-rep(SORC$dLower,SORC$dimension)
  SORC$fpop<-SORC$pop
  SORC$fGres<-SORC$Gres
  SORC$fFres<-SORC$Fres
  SORC$FGEPSres <- cbind(SORC$Fres,SORC$Gres,NA*SORC$Gres) # /WK/ debug only
  defEPS<-SORC$defEPS
  defGamma<-SORC$defGamma
  SORC$EPS<-rep(0,SORC$nConstraints)
  SORC$gamma<-rep(defGamma$gammaInit,SORC$nConstraints)
  SORC$EISOLU<-c()
  SORC$EISA<-c()
  SORC$NU<-rep(1,SORC$nConstraints)
  SORC$DF<-myDF
  fn<-SORC$fn
  
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
  cplugin<-setPluginValue(SORC,bIndex)
  if(is.infinite(cplugin))cplugin<-10^10    # /WK/ shouldn't this be larger, e.g. 10^10 as in the comment above?
  SORC$cplugin<-cplugin
  SORC$bestSolu<-bestI
  SORC$plugin<-c(SORC$plugin,SORC$cplugin)
  
  
  while(SORC$iter < SORC$budget){
    iTime<-proc.time()
    myDF<-SORC$DF
   
    ###################################
    # STEP 2: Train the models        #
    ###################################

    #.Train the objective and constraint functions with RTM function 
    # which stands for(ROBUST TRAINED MODELS)
    RTM<-robustModeling(SORC)
    SORC<-RTM$SORC
    objModel<-RTM$objModel
    con.model<-RTM$con.model
    
    ###################################
    # STEP 3: Optimize                #
    #   the Surrogate problem         #
    ###################################
    
    SORC<-runOptimizer(SORC,objModel,con.model)
    SORC$type<-"SA"
    
    ###################################
    # STEP 4: Update the data.frames  #
    #                                 #
    ###################################
    SORC<-updatingResult(SORC,con.model,objModel)
    
    
   
    ###################################
    #    nu - mechanism               #
    #      deprectaed                 #
    ###################################
    
    if(SORC$NUmech && SORC$DF$numViol[nrow(SORC$DF)]!=0 && SORC$DF$maxViol[nrow(SORC$DF)]<1e-1){
      ##calculating the slope and the yintercept of the EI on the infill point for nu-mechanism
      yintercept<-myEI(SORC$infill,objModel,plugin=SORC$plugin[length(SORC$plugin)],type=SORC$kType)
      Grad<-grad(func=myEI,x=SORC$infill,method="simple",
                 objModel=objModel,
                 plugin=SORC$plugin[length(SORC$plugin)],
                 type=SORC$kType)
      probGrad<-1
      gGrad<-c()
      for(con in c(1:SORC$nConstraints)){
        
        
        probGrad<-grad(func=function(x){
          model<-con.model[[con]]
          x<-t(as.matrix(x))
          colnames(x) = colnames(model@X)
          predx <- predict(object = model, newdata = x, type = SORC$kType, checkNames = FALSE)
          kriging.mean <- predx$mean
          kriging.sd <- predx$sd
          return(kriging.mean/kriging.sd)},x=SORC$infill,method="simple")
       # gGrad<-c(gGrad,as.numeric(Grad%*%probGrad))
        gGrad<-c(gGrad,sqrt(sum(Grad^2)))
        
      }
      
      #gGrad<-1
      slope<-sqrt(sum(Grad^2))
     # SORC$nuCoef<- slope/yintercept
      SORC$nuCoef<-slope/myEI(SORC$infill,objModel=objModel,
                              plugin=SORC$plugin[length(SORC$plugin)],
                              type=SORC$kType)
      SORC$NU<-rep(1,SORC$nConstraints)
      predg<-SORC$predC[nrow(SORC$predC),]
      predgsd<-SORC$sd.con[nrow(SORC$sd.con),]
      gs<-SORC$Gres[nrow(SORC$Gres),]
      prob<-SORC$infill-predg[which(gs>0)]/predgsd[which(gs>0)]
     # SORC$NU[which(gs>0)]<-max(1,SORC$nuCoef*dnorm(prob))
    #  print(paste("NUoriginal",SORC$nuCoef*pnorm(prob)/(dnorm(prob)*abs(gGrad))))
     # SORC$NU[which(gs>0)]<-pmax(1,SORC$nuCoef*pnorm(prob)/(dnorm(prob)*abs(gGrad)))
      SORC$NU[which(gs>0)]<-pmax(1,SORC$nuCoef*pnorm(0)/(dnorm(0)*abs(gGrad)))
      
      if(any(is.na(SORC$NU)))browser()
      if(SORC$verbose>3)print(paste("NU:", SORC$NU)) #debugging message

      if(any(SORC$NU!=1)){
        print("running the optimization process on the repaired nu")
        SORC<-runOptimizer(SORC,objModel,con.model)
        SORC$type<-"NU"
        SORC<-updatingResult(SORC,con.model,objModel)
         }#
    
    
    }# End of the mu mechanism
 
    ###################################
    #    repair - mechanism           #
    #       Not working optimal yet   #
    #                                 #
    ###################################       
    
    if(SORC$repairInfeas && SORC$DF$numViol[nrow(SORC$DF)]!=0 && SORC$DF$maxViol[nrow(SORC$DF)]<SORC$ri$repairMargin ){
      #
      xNewRepaired<-repairInfeasRI2(SORC$infill,SORC$Gres[nrow(SORC$Gres),], 
                                     con.model,SORC,modelType="kriging")
      SORC$infill<-xNewRepaired
      SORC$type<-"REP"
      SORC<-updatingResult(SORC,con.model,objModel)
    }#End of the repair mechanism
    
    
    if((SORC$verbose>0) && ((SORC$iter-SORC$initSize)%%SORC$printIter==0))print(paste("Computation in iteration",SORC$iter," took:",round((proc.time()-iTime)[3],digits=3)))
    
  }# End of the while(SORC$iter < SORC$budget)
  
  return(list(SORC=SORC,objModel=objModel,con.model=con.model))
}# End of the main function



#-------------------Helper functions------------------------------#

setOpts<-function(opts,defaultOpt){
  
  setting<-defaultOpt
  
  if(methods::hasArg(opts)){
    matching<-intersect(names(opts),names(defaultOpt))
    setting[matching] <- opts[matching]
    
    notMatching <- setdiff(names(opts),
                           names(defaultOpt))
    if(length(notMatching)!=0) warning(paste("The
                                             following arguments are ignored: ", notMatching))
    
  }
  
  
  return(setting)
}