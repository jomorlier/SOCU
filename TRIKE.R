#TRIKE.R


TRIKE<-function(SORC){
  lmin<-1e-6 # 1e-7 1e-6 1e-4 0.02 # 
  lmax<-2*(SORC$dUpper-SORC$dLower)
  cf<-1-(1/SORC$dimension) #contraction factor, used to be 0.5 for the test results of G06
  ef<-2                    #  expansion factor, used to be 2 for the results of G06
  eta<-1
  AI<-SORC$AI       # actual improvement
  EI<-SORC$EI       # EImod, modified expected improvement
  lcurrent<-SORC$lcurrent
# if(SORC$iter==11) 
 # browser()
  if(is.na(EI*AI)==F)if((EI)!=0){
 
  if(AI/EI >= eta){
    #expand
    print("expand")
    lcurrent<-min(lcurrent*ef,lmax)
  }else if(AI==0){
    #contract
    #browser()
    print("contract")
    lcurrent<-max(lcurrent*cf,lmin)
    factor<-4/lcurrent
    #SORC$fn<-rescaleWrapper(fn=fn,lower=factor*lower,upper=factor*upper,dimension=dimension,newlower=dLower,newupper=dUpper)
    
  }
    
  }else{
    lcurrent<-min(lcurrent*ef,lmax)
    
  }
  TRU<-SORC$bestSolu+lcurrent/2
  TRU[TRU>SORC$dUpper]<-SORC$dUpper
  
  TRL<-SORC$bestSolu-lcurrent/2
  TRL[TRL<SORC$dLower]<-SORC$dLower
  
  
  SORC$TRU<-as.vector(as.numeric(TRU))
  SORC$TRL<-as.vector(as.numeric(TRL))
  SORC$lcurrent<-lcurrent
  print(paste("lcurrent:",SORC$lcurrent))
 # browser()
  return(SORC)
}