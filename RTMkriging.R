# TODO: add documentation
#
# 
#
RTMKriging<-function(SORC){
  
  if(SORC$handleCrash ){
    objT<-proc.time()
    
    objModel<-tryCatch(
      
      objModel <- km( design = data.frame(SORC$fpop), response = SORC$fFres,
                      covtype = "matern3_2",control=list(trace=F),nugget.estim=SORC$doNUGGET, 
                      noise.var=rep(SORC$noiseVar,nrow(SORC$fpop))),
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
        objModel <- km( design = data.frame(SORC$fpop), response = SORC$fFres,
                        covtype = "matern3_2",control=list(trace=F),nugget.estim=SORC$doNUGGET, 
                        noise.var=rep(SORC$noiseVar,nrow(SORC$fpop)))
      }
    )
    
    
    if((SORC$verbose>0) && ((SORC$iter-SORC$initSize)%%SORC$printIter==0))print(paste("training objective function took:",round((proc.time()-objT)[3],digits=3)))
  }else{
    objModel <- km( design = data.frame(SORC$fpop), response = SORC$fFres,
                    covtype = "matern3_2",control=list(trace=F),nugget.estim=SORC$doNUGGET, 
                    noise.var=rep(SORC$noiseVar,nrow(SORC$fpop)))
    
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
      con.model[con]<- tryCatch(
        #         con.model[con] <- km(design = data.frame(SORC$fpop)
        #                              , response = matrix(as.numeric(SORC$fGres[,con]),nrow=nrow(SORC$fpop))
        #                              , covtype = "matern3_2",control=list(trace=F),nugget.estim=SORC$doConNUGGET
        #                              , noise.var=rep(SORC$conNoiseVar,nrow(SORC$fpop))
        #         )
        {
          conpop <- SORC$fpop[ind,]
          contarget <- SORC$fGres[ind,con]
          con.model[con] <- km(design = data.frame(conpop)
                               , response = matrix(as.numeric(contarget),nrow=nrow(conpop))
                               , covtype = "matern3_2",control=list(trace=F),nugget.estim=SORC$doConNUGGET
                               , noise.var=rep(SORC$conNoiseVar,nrow(conpop))
          )
        },
        
        
        error=function(e){cat("WhooOooO0oOps!! ill conditioned con.model \n We must delete a few points from the population\n");
          closeInd<-which(dist(rbind(SORC$infill,SORC$fpop))<1e-6)
          closeInd<-closeInd[closeInd<=nrow(SORC$fpop)]
          browser()
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
          con.model[con] <- km(   design = data.frame(SORC$fpop)
                                  , response = matrix(as.numeric(SORC$fGres[,con]),nrow=nrow(SORC$fpop))
                                  , covtype = "matern3_2",control=list(trace=F),nugget.estim=SORC$doConNUGGET
                                  , noise.var=rep(SORC$conNoiseVar,nrow(SORC$fpop))
          )
          
        }
      )
      
    }
    if((SORC$verbose>0) && ((SORC$iter-SORC$initSize)%%SORC$printIter==0))print(paste("training constraint function(s) took:",round((proc.time()-conT)[3],digits=3)))
    
  }else{
    con.model<-list()
    conT<-proc.time()
    for(con in c(1:SORC$nConstraints)){
      con.model[con]<-km(   design = data.frame(SORC$fpop)
                            , response = matrix(as.numeric(SORC$fGres[,con]),nrow=nrow(SORC$fpop))
                            , covtype = "matern3_2",control=list(trace=F),nugget.estim=SORC$doConNUGGET
                            , noise.var=rep(SORC$conNoiseVar,nrow(SORC$fpop))
      )
    }
    print(paste("training constraint function(s) took:",round((proc.time()-conT)[3],digits=3)))
  }
  return(RTM=list(SORC=SORC,objModel=objModel,con.model=con.model)) 
}