#TRIKEmini.R


TRIKEmini<-function(p1,p2,SORC){
  lmin<-1e-6 # 1e-7 1e-6 1e-4 0.02 # 
  lmax<-2*(SORC$dUpper-SORC$dLower)
  tauc<-0.4
  cf<-1-(1/SORC$dimension) #contraction factor, used to be 0.5 for the test results of G06
  ef<-2                    #  expansion factor, used to be 2 for the results of G06
  lcurrent<-SORC$lcurrent
  delta<-dist(rbind(p1,p2))
  state<-"stay";
  if (!is.null(SORC$p2prev)) {
    eta<-dist(rbind(p2,SORC$p2prev))
    if (eta>lcurrent) state<-"expandEta";
  }
  if (any(p2==SORC$TRL)) state<-"expand";
  if (any(p2==SORC$TRU)) state<-"expand";
  if (delta/lcurrent<tauc) state<-"contract";
  
  switch(state,
         contract = {lcurrent<-max(lcurrent*cf,lmin);},
         expand   = {lcurrent<-min(lcurrent*ef,lmax);},
         expandEta= {lcurrent<-eta;},
         stay     = {}
  )
  
  SORC$lcurrent<-lcurrent
  SORC$p2prev<-p2 # for next iteration
  cat(sprintf("state: %s,  lcurrent: %g\n",state,SORC$lcurrent))
 # browser()
  return(SORC)
}