


makexMat<-function(dd, TRL, TRU,SORC,xType="optim"){
  
 
  
  x1<-seq(from=TRL[dd[1]],to=TRU[dd[1]],length=100)
  x2<-seq(from=TRL[dd[2]],to=TRU[dd[2]],length=100)
  xMat0<-expand.grid(x=x1,y=x2) 
  switch(xType,
         "optim"={xMat<-matrix(SORC$solu,nrow=nrow(xMat0),ncol=SORC$dimension,byrow=TRUE)}
        # "zero"={xMat<-matrix(0,nrow=nrow(xMat0),ncol=SORC$dimension)},
        # "infill"={xMat<-matrix(SORC$infill,nrow=nrow(xMat0),ncol=SORC$dimension,byrow=TRUE)}
         )
#   xMat<-matrix(0,nrow=nrow(xMat0),ncol=SORC$dimension)
#   xMat[,dd[1]]<-xMat0[,1]
#   xMat[,dd[2]]<-xMat0[,2]
  xMat[,dd[1]]<-xMat0[,1];
  xMat[,dd[2]]<-xMat0[,2]
  # xMat<-cbind(xMat,X)
  #browser()
  xMat<-data.frame(xMat)
  names(xMat)<-paste("x",c(1:SORC$dimension),sep="")
  return(xMat)
}

makeDataFrame<-function(dd, TRL, TRU,SORC,xType,func,objModel,con.model){
  xMat<-makexMat(dd, TRL, TRU,SORC,xType)
  switch(func,
         "EI"={  z<-myEI(xMat,objModel,type=SORC$kType,plugin=SORC$plugin[length(SORC$plugin)])},
         "EImod"={z<--EImod(xMat,con.model,objModel,SORC)},
         "FF"={z<-CalculateF(xMat,con.model,SORC)$res },
         "Fideal"={ mean.sd <- CalculateF(makexMat(dd, TRL, TRU,SORC,xType),con.model,SORC)$mean.sd;     
         cat(sprintf("mean.sd : %s\n",paste(mean.sd,collapse=" ")));
         z <- CalculateF_ideal(makexMat(dd, TRL, TRU,SORC,xType),ideal.sd=mean.sd,SORC)$res})
 # browser()
  xMat<-data.frame(xMat[,dd],z)
  names(xMat)<-c("x","y","z")
  return(xMat)
  
}

visualizenD<-function(objModel,fn,SORC,con.model,TRL,TRU,dd){
  
  THREE_PLOTS=T         # if TRUE, arrange three plots: EI, F, EImod
  # if FALSE, arrange two plots:  EI, EImod
  OBJ_CON1_NEEDED=F     # if TRUE, calculate objP, con1P, con2P (they seem not needed anymore) 
  DEBUG_OLD_VERSION=F   # if TRUE, run also the old version for EI.value and EImod.value 
  # (slower by factor >100) and compare
  IDEAL_F=F             # debug plot of true F, can be very slow (additional 3 sec)
  
  num<-nrow(SORC$pop)
  dimension<-ncol(SORC$pop)
  
  if(SORC$visZoom){
#     dev<-SORC$iter
#     pow<-1.3 # 1.8 2 
#     delta<-max(0.0002,0.1/dev^pow)
#     TRL<-SORC$solu-delta          # /WK/ limit the zoom-in
#     TRU<-SORC$solu+delta  
#     # 
#     TRL<-pmin(TRL,SORC$bestSolu)
#     TRU<-pmax(TRU,SORC$bestSolu)
  }

  DD<-NULL
  
  for(i in c(1:(dimension-1))){
    for(j in c((i+1):dimension)){
      DD<-rbind(DD,c(i,j))
    }
  }
  
  
 # bestsolu<-data.frame(SORC$pop[SORC$bIndex,],z=NA)
 # names(bestsolu)<-c("x","y","z")
  #browser()
  manipulate({
    PLOTMAT<-as.list(matrix(0,nrow=dimension,ncol=dimension))
    for(i in c(1:nrow(DD))){
      dd<-as.numeric(DD[i,])
      #browser()
    
       bestsolu<-data.frame(SORC$pop[SORC$bIndex,dd],z=NA)
       closesolu<-data.frame(t(c(SORC$closestSOLU[dd],z=NA)))
       #browser()
       names(bestsolu)<-c("x","y","z")
       names(closesolu)<-c("x","y","z")
       
    if(!is.null(SORC$solu)){
      if(is.matrix(SORC$solu)) {soluMat<-data.frame(SORC$solu[,dd[1]],SORC$solu[,dd[2]],SORC$optim)
      }else{
        soluMat<-data.frame(SORC$solu[dd[1]],SORC$solu[dd[2]],SORC$optim)  
      }
      # browser()
      names(soluMat)<-c("x","y","z")
    }  
    
  EIP<-ggplot(data=makeDataFrame(dd, TRL, TRU,SORC,xType,func="EI",objModel,con.model),aes(x=x,y=y,z=z))
  EIP<-EIP+geom_tile(aes(fill=z))+stat_contour(color="black")
 # browser()
  EIP<-EIP+ scale_fill_gradientn(name="EI(obj)",colours=brewer.pal(7,"Greens"))
  EIP<-EIP+ geom_point(data=bestsolu,aes(x=x,y=y),color="blue",size=4)
  EIP<-EIP+ geom_point(data=closesolu,aes(x=x,y=y),color="pink",size=4)
  
  datapoint<-SORC$pop[,dd]
  names(datapoint)<-c("x1","x2")
  EIP<-EIP+geom_point(data=data.frame(datapoint,z=NA),aes(x=x1,y=x2))
  EIP<-EIP+coord_cartesian(xlim=c(TRL[1],TRU[1]),ylim=c(TRL[2],TRU[2]))
  PLOTMAT[[(dd[1]-1)*dimension+dd[2]]]<-EIP
#   FIP<-ggplot(data=makeDataFrame(dd, TRL, TRU,SORC,xType,func="FF",objModel,con.model),aes(x=x,y=y,z=z))
#   FIP<-FIP+geom_tile(aes(fill=z))+stat_contour(color="black")
#   FIP<-FIP+ scale_fill_gradientn(name="F(con)",colours=brewer.pal(6,"YlOrRd"))
#   FIP<-FIP+ geom_point(data=soluMat,aes(x=x,y=y),fill="yellow",size=4,shape=22)
#   FIP<-FIP+ geom_point(data=data.frame(bestsolu,z=NA),aes(x=x,y=y),color="blue",size=4)
#   FIP<-FIP+geom_point(data=data.frame(SORC$pop,z=NA),aes(x=x1,y=x2))
#   FIP<-FIP+coord_cartesian(xlim=c(TRL[1],TRU[1]),ylim=c(TRL[2],TRU[2]))
#   
#   if (IDEAL_F) {
#     F2P<-ggplot(data=makeDataFrame(dd, TRL, TRU,SORC,xType,func="Fideal",objModel,con.model),aes(x=x,y=y,z=z))
#     F2P<-F2P+geom_tile(aes(fill=z))+stat_contour(color="black")
#     F2P<-F2P+ scale_fill_gradientn(name="F(ideal)",colours=brewer.pal(6,"YlOrRd"))
#     F2P<-F2P+ geom_point(data=soluMat,aes(x=x,y=y),fill="yellow",size=4,shape=22)
#    # F2P<-F2P+ geom_point(data=bestsolu,aes(x=x,y=y),color="red",size=4)
#     F2P<-F2P+geom_point(data=data.frame(SORC$pop,z=NA),aes(x=x1,y=x2))
#     F2P<-F2P+coord_cartesian(xlim=c(TRL[1],TRU[1]),ylim=c(TRL[2],TRU[2]))
#   } else {
#     F2P<-ggplot(data=makeDataFrame(dd, TRL, TRU,SORC,xType,func="FF",objModel,con.model),aes(x=x,y=y,z=z))
#     F2P<-F2P+geom_tile(aes(fill=z))+stat_contour(color="black")
#     F2P<-F2P+ scale_fill_gradientn(name="F(con)",colours=brewer.pal(6,"YlOrRd"))
#     F2P<-F2P+ geom_point(data=soluMat,aes(x=x,y=y),fill="yellow",size=4,shape=22)
#    # F2P<-F2P+ geom_point(data=bestsolu,aes(x=x,y=y),color="red",size=4)
#     F2P<-F2P+geom_point(data=data.frame(SORC$pop,z=NA),aes(x=x1,y=x2))
#     F2P<-F2P+coord_cartesian(xlim=c(TRL[1],TRU[1]),ylim=c(TRL[2],TRU[2]))
#   }
  
  
 # SORC$EImax<-max(-EI.mod)
  EImodP<-ggplot(data=makeDataFrame(dd, TRL, TRU,SORC,xType,func="EImod",objModel,con.model),aes(x=x,y=y,z=z))
  EImodP<-EImodP+geom_tile(aes(fill=z))+stat_contour(color="black") 
  EImodP<-EImodP+ scale_fill_gradientn(name="EImod",colours=brewer.pal(6,"YlOrRd"))
  EImodP<-EImodP+ geom_point(data=soluMat,aes(x=x,y=y),fill="yellow",size=4,shape=22)
  EImodP<-EImodP+ geom_point(data=bestsolu,aes(x=x,y=y),color="blue",size=4)
  EImodP<-EImodP+ geom_point(data=closesolu,aes(x=x,y=y),color="pink",size=4)
  
  EImodP<-EImodP+geom_point(data=data.frame(datapoint,z=NA),aes(x=x1,y=x2))
  EImodP<-EImodP+coord_cartesian(xlim=c(TRL [1],TRU[1]),ylim=c(TRL[2],TRU[2]))
 
  #PLOTMAT[rev(dd)]<-EImodP
  PLOTMAT[[(dd[2]-1)*dimension+dd[1]]]<-EImodP
  
  #browser()

    }
    for(i in seq(1,dimension^2,(dimension+1))){
        p<-ggplot(data=data.frame(x=c(1:10),y=c(1:10),z=c(1:10)),aes(x=x,y=y))+geom_line()
        annText<-paste("Best Sol [",SORC$iter,"]: ", round(SORC$Fres[SORC$bIndex],digits=2)," | ",round(SORC$DF$maxViol[SORC$bIndex],digits=3),sep="")
        
        annText0<-paste("                     obj","     ", "  maxViol \n",annText)
        if(any(SORC$NU!=1)){  
          annText0<-paste("    adjusting nu     obj","     ", "  maxViol \n",annText)}
        
        #browser()
        
        lly<-TRU[2]-TRL[2]
        llx<-TRU[1]-TRL[1]
        
        textFrame<-data.frame(x=5,y=5,z=annText0)
        
        p<-p+geom_label(data=textFrame,aes(x=x,y=y,label=z),alpha=0.5,size=5)
        
        PLOTMAT[[i]]<-p
          
      
    }
   # browser()
    
   # newp<-grid.arrange(PLOTMAT[c(1:dimension^2)], ncol = dimension)
    newp<-do.call(grid.arrange, PLOTMAT )
    ggsave(paste("../Figures/",SORC$name,"-",SORC$fCalc,"-EI",num,".pdf",sep=""),newp,width=25,height=25)
    
  },
  xType=picker("optim","infill","zero")
  )
 # browser()
  return(SORC)
}