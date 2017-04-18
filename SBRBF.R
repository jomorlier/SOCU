#'
#'
#'
#'SBRBF.R script provides a non-linear regression 
#'This script build an RBF Ridge model for a multivariate X element of R^d space 
#'according to least square method
#'This regressor has several parameter \code{delta}, \code{param, \code{kernel} and \code{centroid}
#'
#'15.02.2017
#'Samineh Bagheri
#'Cologne University of Applied Sciences
#'
#'The code contains two functions \code{trainSBRBF()} and \code{predictSBRBF()}
#'
#'
#'

#'Input:
#' @param X: is a matrix of size n x d 
#' where n is the size of the population (size of the sample points) 
#' and d is the dimension of the input space
#' @param Y: is often a matrix of size n x 1 represents the labels or the output of X
#' @param delta: is a scalar value often used for regulatization, it can be set to zero if the data is noise-free 
#' and the training sample is not very large
#' @param kernel[Gauss]: 
#' @param param:
#' @param delta:
#' @param DISTM[euclidean]: method of calculating distance, it is recommended to use \code{"manhattan"} distance for dimensions higher than 3
#' Output:
#' @param  SBRIDGE.model: which constains the model coeffecients theta 
#' theta is a matrix of size (d+1) X 1

train<-trainSBRBF<-function(X,Y,delta,kernel="Cubic",param=1,centroid=(1:nrow(X)),DISTM="euclidean",ptail=TRUE,squares=TRUE){
  standardize<-function(x){
    testit::assert("x is a vector",length(x)>1)
    y<-(x-mean(x))/sd(x)
    return(y)
  }
  switch(kernel,
         "Gauss" ={  phiFunc<-function(DIST){
           apply(DIST, c(1,2),function(x){exp(-param*x^2)})
         }},
         "Cubic"={  phiFunc<-function(DIST){
           DIST*DIST*DIST
         }}
  )
  
  DIST<-as.matrix(stats::dist(X,upper=T,method = DISTM))[,centroid]

  K<-phiFunc(DIST)
  
  P<-NULL

  if(ptail)P<-cbind(P,1,X)
  if(squares)P<-cbind(P,X^2)
  if(squares||ptail){
    zeros<-matrix(0,nrow=ncol(P),ncol=ncol(P))
    K<-cbind(K,P)
    K<-rbind(K,cbind(t(P),zeros)) 
    Y<-rbind(Y,matrix(0,nrow=ncol(P),ncol=1))
  }
  
  n<-ncol(K)
  KK<-t(K)%*%K + delta*diag(x = 1, nrow=n, ncol=n)
  s = svd(KK)
  invD <- 1/s$d
  invD[abs(s$d/s$d[1])<1e-14] <- 0
  InversedTerm <- s$v %*% diag(invD) %*% t(s$u)
  #InversedTerm<-MASS::ginv(KK)
  rhs<-t(K) %*% Y
  theta<-InversedTerm %*% rhs
  RN<-paste("theta",c(1:ncol(K)),sep="")
  rownames(theta)<-RN
  SBRBF.model<-list(theta=theta,
                    phi=K,
                    InversedTerm=InversedTerm,
                    delta=delta,
                    name="SBRBF",
                    kernel=kernel,
                    param=param,
                    input=X,
                    DISTM=DISTM,
                    ptail=ptail,
                    squares=squares,
                    xc=matrix(X[centroid,],nrow=length(centroid),ncol=ncol(X)))
  return(SBRBF.model)
}

#newx can be a matrix of size m*d where m is the number of new points 
#and d is the dimension of x

predictSBRBF<-function(newdata,model,type = NULL, 
                       checkNames = FALSE){
  
  testit::assert("You are calling SBRBF predict function but your model is from some other type",model$name=="SBRBF")
  ptail<-model$ptail
  squares<-model$squares
  
  newx<-newdata
  if(!is.matrix(newx))newx<-as.matrix(newx,ncol=length(newx))

  standardize<-function(x){
    testit::assert("x is a vector",length(x)>1)
    y<-(x-mean(x))/sd(x)
    return(y)
  }
  switch(model$kernel,
         "Gauss" ={  phiFunc<-function(DIST){
           apply(DIST, c(1,2),function(x){exp(-model$param*x^2)})
         }},
         "Cubic"={  phiFunc<-function(DIST){
           DIST*DIST*DIST
         }}
  )

    DISTnew<-t(apply(newx,1,function(x){as.matrix(dist(rbind(x,model$xc),method=model$DISTM))[-1,1]}))
    phis<-matrix(phiFunc(DISTnew),nrow=nrow(DISTnew))
    
    P<-NULL
    if(ptail)P<-cbind(P,1,newx)
    if(squares)P<-cbind(P,newx^2)
    if(squares||ptail){
      zeros<-matrix(0,nrow=ncol(P),ncol=ncol(P))
      phis<-cbind(phis,P)
     # phis<-rbind(phis,cbind(t(P),zeros)) 
      #Y<-rbind(Y,matrix(0,nrow=ncol(P),ncol=1))
    }
  # browser()
    
    n<-nrow(model$input)
    m<-nrow(newx)
    d<-ncol(newx)
    phiss<-matrix(phiFunc(matrix(0,nrow=m,ncol=1)),nrow=m)
    
    mu<-phis%*%model$theta
    sigma<-phiss-(phis%*%model$InversedTerm%*%t(model$phi)*phis)%*%matrix(rep(1,ncol(phis)),nrow=ncol(phis))
   # sigma[which(abs(sigma)<1e-6)]<-1e-6
    sigma[which(sigma<0)]<-1e-6
    if(any(sigma<0))browser()
    output<-list(mu=mu,sigma=sqrt(sigma))
  return(output)
}
predict<-predictSBRBF
  