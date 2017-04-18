
#default values of SBRBF function

defaultSBRBF<-function(){
  defSBRBF<-list(
               active=F
               ,kernel="Cubic" #Gauss, Cubic
               ,param=1
               ,delta=0.0
               ,DISTM="euclidean" # "manhattan"
               ,ptail=TRUE
               ,squares=TRUE
  )
  return(defSBRBF)
}