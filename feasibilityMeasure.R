#'
#'
#'
#'
#'Samineh Bagheri
#'Cologne University of Applied Sciences
#'
#'

# /WK/ not used anywhere - delete? 
plog<-function(y,pShift=0.0){
  if(y-pShift>=0){
    ret<- log(1+y-pShift)
  }else{
    ret<- -log(1-(y-pShift))
  }
  return(ret)
}


# /WK/ own EI (expected improvement function), since the function EI (package DiceOptim) 
# sometimes returns exactly zero (!) when sd is smaller than a certain value, although the 
# calculation in myEI shows that the formula gives an EI>0. 
#
# (For larger sd, both myEI and EI return exactly the same results)
#
myEI <- function(x,objModel,plugin=NULL,type="UK") {
  if (!is.data.frame(x)) x <- t(x)
  predx = predict(object = objModel, newdata = x, type = type, checkNames = FALSE)
  yhat<-predx$mean
  sd  <-predx$sd
  arg<-(plugin-yhat)/sd
  ei<-(plugin-yhat)*pnorm(arg) + sd*dnorm(arg)
  #ei[sd==0] <- 0.0
  
  #SB: only for testing 
  
 # minEI<-1e-09
 # ei[ei < minEI] <- ei[ei < minEI]+minEI
  #return(-predx$mean)  # just a debug switch to test the hypothesis whether objective 
                        # Kriging model is as good as EI for our purposes /WK/
  return(ei)
}

# /WK/ new version (old one in CalculateF_OLD), where x can be either a vector (single point)
# or a data.frame (multiple points in rows).
# Returns a list(res,alpha) instead of the former concatenated vector
#
CalculateF<-function(x,con.model,SORC){
  
  if(!is.data.frame(x)) {
    newdata.num <- as.numeric(x)
    newdata <- data.frame(t(newdata.num))
  } else {
    newdata <- x
  }

  tol<-SORC$tol
  alg<-SORC$fCalc
  wt<-SORC$wt
  EPS<-SORC$EPS
  alpha<-1
  beta<-1
  mean.sd <- NULL
  
  #SB:nu mechanism
 # browser()
#   yintercept<-myEI(newdata,objModel,plugin=SORC$plugin,type=SORC$kType)
#   Grad<-grad(func=myEI,x=newdata,method="simple",
#              objModel=objModel,
#              plugin=SORC$plugin,
#              type=SORC$kType)
#   slope<-sqrt(sum(Grad^2))
#   nuCoef<- slope/yintercept
  #browser()
  
  for(con in c(1:SORC$nConstraints)){
    
    model<-con.model[[con]]
    if(!SORC$useSBRBF$active)colnames(newdata) = colnames(model@X)
   # browser()
    predx <- predict(object = model, newdata = newdata, type = SORC$kType, 
                     checkNames = FALSE)
    kriging.mean <- predx$mean
    
    kriging.sd <- predx$sd
   # xcrAlpha <- (tol - kriging.mean)/kriging.sd
   # xcrAlpha <- (-EPS[con] - kriging.mean)/kriging.sd #Do we need EPS[con] here???
    xcrAlpha <- ( - kriging.mean)/kriging.sd 
    #browser()
    xcrA.prob <- pnorm(xcrAlpha)
    #xcrA.prob[kriging.sd==0] <- 1       # /WK/ to fix the case sd==0 (assumed to be surely feasible)
    alpha<-alpha*xcrA.prob
    
    #xcrBeta <- (tol - kriging.mean)/kriging.sd
    if(is.null(SORC$constantSigma)){
      xcrBeta <- (-EPS[con] - kriging.mean)/kriging.sd
      
    }else{
      xcrBeta <- (-EPS[con] - kriging.mean)/SORC$constantSigma
      
    }

    xcrB.prob <- pnorm(xcrBeta)

    beta<-beta*(pmin(1,SORC$gamma[con]*xcrB.prob))^SORC$NU[con]  # why pmin? - if x is a data.frame, xcrB.prob is a vector
  
    mean.sd=c(mean.sd,mean(kriging.sd))
  }
  switch(alg,
         "alpha"={res<-alpha},
         "beta"={res<-beta},
         "mix"={res<-wt*alpha+(1-wt)*beta},
         stop("Invalid fCalc algorithm"))
  return(list(res=res,alpha=alpha,mean.sd=mean.sd)) 
}

# this is only a function for debug purposes: How would the ideal F look like? 
# Instead of calling the constraint Kriging models, we call the true constraint functions.
# We need to pass in an assumed standard deviation ideal.sd (either a number or a vector 
# with length equal to number of constraints)
#
CalculateF_ideal<-function(x,ideal.sd=1,SORC){
  
  if(!is.data.frame(x)) {
    newdata.num <- as.numeric(x)
    newdata <- data.frame(t(newdata.num))
  } else {
    newdata <- x
  }
  ideal.sd <- mean(ideal.sd) # /WK/ experimental
  if (length(ideal.sd)==1) ideal.sd <- rep(ideal.sd,SORC$nConstraints)
  alg<-SORC$fCalc
  wt<-SORC$wt
  alpha<-1
  beta<-1
  
  fnMat <- sapply(1:nrow(newdata),function(i)SORC$fn(as.vector(t(newdata[i,]))))
  fnMat <- fnMat[-1,]   # strip the row with objective, retain only g1..gn

  for(con in c(1:SORC$nConstraints)){
    ideal.mean <- fnMat[con,]
    
    xcrAlpha <- ( - ideal.mean)/ideal.sd[con]
    xcrA.prob <- pnorm(xcrAlpha)
    alpha<-alpha*xcrA.prob
    
    xcrBeta <- ( - ideal.mean)/ideal.sd[con]
    
    xcrB.prob <- pnorm(xcrBeta)
    beta<-beta*pmin(1,SORC$gamma*xcrB.prob)   # why pmin? - if x is a data.frame, xcrB.prob is a vector
  }
  
  switch(alg,
         "alpha"={res<-alpha},
         "beta"={res<-beta},
         "mix"={res<-wt*alpha+(1-wt)*beta},
         stop("Invalid fCalc algorithm"))
  return(list(res=res,alpha=alpha)) 
}

# /WK/ debug only: calculate the constraint model response for x
#
CalculateCon<-function(x,con.model,SORC)
{
  
  if(!is.data.frame(x)) {
    newdata.num <- as.numeric(x)
    newdata <- data.frame(t(newdata.num))
  } else {
    newdata <- x
  }
  
  tol<-SORC$tol
  mean.sd <- NULL
  conVec <- NULL
  conSD <- NULL
  
  for(con in c(1:SORC$nConstraints)){
    
    model<-con.model[[con]]
    if(!SORC$useSBRBF$active)colnames(newdata) = colnames(model@X)
    predx <- predict(object = model, newdata = newdata, type = SORC$kType, 
                     checkNames = FALSE)
    kriging.mean <- predx$mean
    kriging.sd <- predx$sd
    
    conVec <- c(conVec,kriging.mean)
    conSD <- c(conSD,kriging.sd)
  }
  return(list(conVec=conVec,conSD=conSD)) 
}



# /WK/ new version (old one in EImod_OLD), where CalculateF (x can be either vector or data
# frame) instead of CalculateF_OLD is called and where the optional new parameter SORC$nu
# balances the relative weight between EI and feasibility probability. 
#
EImod<-function(x,con.model,objModel,SORC){

  if(any(SORC$DF$numViol==0) ){
    if(SORC$useEI){
     # browser()
      objEI<-myEI(x,objModel,plugin=SORC$plugin[length(SORC$plugin)],type=SORC$kType)  # /WK/ new, to cure the buggy EI
      #objEI<-EI(x,objModel,plugin=SORC$plugin,type=SORC$kType)
      obj<-objEI 
    }else{
      if(class(x)=="numeric"){
        newdata <- data.frame(t(x))
        
      }else{
        newdata <- data.frame(x)
        
      }
      if(!SORC$useSBRBF$active)colnames(newdata) = colnames(objModel@X)
      obj<-predict(object = objModel, newdata =newdata , type = SORC$kType, 
                   checkNames = FALSE)$mean
    }
    
    conEI<-CalculateF(x,con.model,SORC)$res
    
   # conEI<-CalculateF(x,con.model,SORC)
    
  }else{
   # print("maximizing the feasibility probability")
    obj<-1
    conEI<-CalculateF(x,con.model,SORC)$alpha
    
  }
  myEImod<-(obj^SORC$nu)*(conEI)

  return(-myEImod)  
}

# /WK/ CalculateF_OLD and EImod_OLD are the old versions of CalculateF and EImod, 
# which make due to their apply-interface in 2DVis.R (and due to the fact that x can 
# only be a vector, not a data frame) the 2D-visualization slower by a factor of 130..150
#
CalculateF_OLD<-function(x,con.model,SORC){
  
  newdata.num <- as.numeric(x)
  newdata <- data.frame(t(newdata.num))
  colnames(newdata) = paste("x",c(1:SORC$dimension),sep="")
  tol<-SORC$tol
  alg<-SORC$fCalc
  wt<-SORC$wt
  alpha<-1
  beta<-1
  for(con in c(1:SORC$nConstraints)){
    model<-con.model[[con]]
    newdata.num <- as.numeric(x)
    newdata <- data.frame(t(newdata.num))
    if(!SORC$useSBRBF$active)colnames(newdata) = colnames(model@X)
    predx <- predict(object = model, newdata = newdata, type = SORC$kType, 
                     checkNames = FALSE)
    kriging.mean <- predx$mean
    
    kriging.sd <- predx$sd
    xcrAlpha <- (tol - kriging.mean)/kriging.sd
    xcrA.prob <- pnorm(xcrAlpha)
    #xcrA.prob[kriging.sd==0] <- 1       # /WK/ to fix the case sd==0 (assumed to be surely feasible)
    xcrA.dens <- dnorm(xcrAlpha)
    alpha<-alpha*xcrA.prob
    
    xcrBeta <- (tol - kriging.mean)/kriging.sd
    xcrB.prob <- pnorm(xcrBeta)
    #xcrB.prob[kriging.sd==0] <- 1       # /WK/ to fix the case sd==0 (assumed to be surely feasible)
    xcrB.dens <- dnorm(xcrBeta)
    beta<-beta*min(1,SORC$gamma*xcrB.prob)  
    #browser()
    
  }
  
  switch(alg,
         "alpha"={res<-alpha},
         "beta"={res<-beta},
         "mix"={res<-wt*alpha+(1-wt)*beta},
         stop("Invalid fCalc algorithm"))
  #browser()
  return(c(res=res,alpha=alpha)) 
}

EImod_OLD<-function(x,con.model,objModel,SORC){
  
  if(any(SORC$DF$numViol==0) ){
    objEI<-myEI(x,objModel,plugin=SORC$plugin[length(SORC$plugin)],type=SORC$kType)  # /WK/ new, to cure the buggy EI
    #objEI<-EI(x,objModel,plugin=SORC$plugin,type=SORC$kType)
    conEI<-CalculateF_OLD(x,con.model,SORC)[1]
    
    # conEI<-CalculateF_OLD(x,con.model,SORC)
  }else{
    objEI<-1
    conEI<-CalculateF_OLD(x,con.model,SORC)[2]
    
  }
  myEImod<-objEI*conEI
  
  return(-myEImod)  
}

# # /WK/ very experimental function
# gridSearch <- function(#fn=function(x,con.model,objModel,SORC){
# #                             return(EImod(x,con.model,objModel,SORC))
# #                           },
#                         lower,
#                         upper,
#                         control=list(delta=0.1,length=10),
#                         con.model=con.model,objModel=objModel,SORC=SORC)
# {
#   dim=length(lower)
#   if (dim>4) stop("gridSearch works currently only for length(lower)<=4")
#   x1 <- x2 <- x3 <- x4 <- 1
#   x = list()
#   for (i in 1:4) {
#     if (i<=dim) {
#       x[[i]] <- seq(from=lower[i],to=upper[i],length=control$length)
#     } else {
#       x[[i]] <- 1
#     }
#   }
#   xMat<-expand.grid(x[[1]],x[[2]],x[[3]],x[[4]])
#   #ei <- -EImod(xMat,con.model,objModel,SORC)
#   return(xMat)
# }


