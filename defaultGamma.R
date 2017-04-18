defaultGamma<-function(){
  defGamma<-list(
    #active=F,  
    gammaInit=2,
    gammaFinal=1,
    gammaScheme="constant",  #"adaptive", "linear", "constant",
    kappa=1,
    gammaThreshold=3
  )
}
