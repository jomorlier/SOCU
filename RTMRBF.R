# TODO: add documentation
#
# 
#
RTMRBF<-function(SORC){
  
  if(SORC$handleCrash ){
    objT<-proc.time()
   # browser()
    objModel<-tryCatch(
      objModel <- train( X = as.matrix(data.frame(SORC$fpop)), Y = as.matrix(SORC$fFres),delta=SORC$noiseVar,
                         kernel=SORC$useSBRBF$kernel,param=SORC$useSBRBF$param,DISTM=SORC$useSBRBF$DISTM,ptail=SORC$useSBRBF$ptail,squares=SORC$useSBRBF$squares),
      error=function(e){cat("WhooOooO0oOps!! ill conditioned model \n We must delete a few points from the population\n");
        closeInd<-which(dist(rbind(SORC$infill,SORC$fpop))<1e-6)
        closeInd<-closeInd[closeInd<=nrow(SORC$fpop)]
        
        if(length(closeInd)>0){
          print(paste("deleting", length(closeInd) ,"point(s)"))
          delInd<-closeInd
          SORC$fpop<<-SORC$fpop[-delInd,]
          SORC$fGres<<-as.matrix(SORC$fGres[-delInd,])
          SORC$fFres<<-SORC$fFres[-delInd]
        };
        objModel <- train( X = as.matrix(data.frame(SORC$fpop)), Y = as.matrix(SORC$fFres),delta=SORC$noiseVar,
                           kernel=SORC$useSBRBF$kernel,param=SORC$useSBRBF$param,DISTM=SORC$useSBRBF$DISTM,ptail=SORC$useSBRBF$ptail,squares=SORC$useSBRBF$squares)
      }
    )
    
    
    if((SORC$verbose>0) && ((SORC$iter-SORC$initSize)%%SORC$printIter==0))print(paste("training objective function took:",round((proc.time()-objT)[3],digits=3)))
  }else{
    objModel <- train( X = as.matrix(data.frame(SORC$fpop)), Y = as.matrix(SORC$fFres),delta=SORC$noiseVar,
                       kernel=SORC$useSBRBF$kernel,param=SORC$useSBRBF$param,DISTM=SORC$useSBRBF$DISTM,ptail=SORC$useSBRBF$ptail,squares=SORC$useSBRBF$squares)
    
  }
  
  # TODO: take for the m constraint functions only those points from the population with maxViol<vEps
  # This should make the constraint models more precise near the best solution.
  #
  # ---- does not work yet, needs more testing ---
  #
  SELECT_CON_POP=FALSE
  if (SELECT_CON_POP) {
    if (SORC$iter > 0.2*SORC$budget) {
      # take those with maxViol<vEps.  
      vEps=1e-3
      ind <- which(SORC$DF$maxViol<vEps)
    } else {
      ind <- 1:nrow(SORC$fpop)
    }
    cat(sprintf("SELECT_CON_POP: Population size: %5d\n",length(ind)))
  } else {
    ind <- 1:nrow(SORC$fpop)
  }
  
  #.Train all m constraint functions
  if(SORC$handleCrash ){
    con.model<-list()
    conT<-proc.time()
    for(con in c(1:SORC$nConstraints)){
      con.model[[con]]<- tryCatch(
        {
          conpop <- SORC$fpop[ind,]
          contarget <- SORC$fGres[ind,con]
          
          con.model[[con]] <- train( X = as.matrix(data.frame(SORC$fpop)), Y = as.matrix(as.numeric(contarget),nrow=nrow(conpop)),delta=0.0,
                                   kernel=SORC$useSBRBF$kernel,param=SORC$useSBRBF$param,DISTM=SORC$useSBRBF$DISTM,ptail=SORC$useSBRBF$ptail,squares=SORC$useSBRBF$squares)
          #browser()
        },
        
        
        error=function(e){cat("WhooOooO0oOps!! ill conditioned con.model \n We must delete a few points from the population\n");
          closeInd<-which(dist(rbind(SORC$infill,SORC$fpop))<1e-6)
          closeInd<-closeInd[closeInd<=nrow(SORC$fpop)]
          #print(paste("population size is:",(nrow(SORC$fpop)-nrow(SORC$pop))))
          
          if(length(closeInd)>0){
            print(paste("deleting", length(closeInd) ,"point(s)"))
            delInd<-closeInd
            SORC$fpop<<-SORC$fpop[-delInd,]
            SORC$fGres<<-as.matrix(SORC$fGres[-delInd,])
            SORC$fFres<<-SORC$fFres[-delInd]
          };
          # browser()
          #print(paste("population size is:",(nrow(SORC$fpop)-nrow(SORC$pop))))
          con.model[[con]] <- train( X = as.matrix(data.frame(SORC$fpop)), Y = as.matrix(as.numeric(contarget),nrow=nrow(conpop)),delta=0.0,
                                   kernel=SORC$useSBRBF$kernel,param=SORC$useSBRBF$param,DISTM=SORC$useSBRBF$DISTM,ptail=SORC$useSBRBF$ptail,squares=SORC$useSBRBF$squares)
          
        }
      )
      
    }
    if((SORC$verbose>0) && ((SORC$iter-SORC$initSize)%%SORC$printIter==0)) print(paste("training constraint function(s) took:",round((proc.time()-conT)[3],digits=3)))
    
  }else{
    con.model<-list()
    conT<-proc.time()
    for(con in c(1:SORC$nConstraints)){
      con.model[[con]]<-train( X = as.matrix(data.frame(SORC$fpop)), Y =matrix(as.numeric(SORC$fGres[,con]),nrow=nrow(SORC$fpop)),delta=0.0,
             kernel=SORC$useSBRBF$kernel,param=SORC$useSBRBF$param,DISTM=SORC$useSBRBF$DISTM,ptail=SORC$useSBRBF$ptail,squares=SORC$useSBRBF$squares)
    }
    print(paste("training constraint function(s) took:",round((proc.time()-conT)[3],digits=3)))
  }
  
  return(RTM=list(SORC=SORC,objModel=objModel,con.model=con.model)) 
}