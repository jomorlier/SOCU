# TODO: add documentation
#
# /WK/: It is possible an error that SORC is a return parameter of robustModeling, but at the same 
# time this function modifies SORC at the higher level with <<- operator. 
#
robustModeling<-function(SORC){
  if(!SORC$useSBRBF$active)
    RTM<-RTMKriging(SORC)
  else
    RTM<-RTMRBF(SORC)
 return(RTM)
}