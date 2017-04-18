readndPlotPer <- function()
{ 
  n <- readline(prompt="Do you want to plot a high dimensional problem? yes [y], No [n] ")
  n <- as.character(n)
  if(n=="y"){
    output<-T
  }else{
    output<-F
  }
  return(output)
}
askDesiredDimension<-function(){
  dd1<-readline(prompt = "please enter the first desired dimensions for plotting: ")
  dd<-as.numeric(dd1)
  dd2<-readline(prompt = "please enter the second desired dimensions for plotting: ")
  dd<-c(dd,as.numeric(dd2))
  return(dd)
}



visualizeTwo<-function(objModel,fn,SORC,con.model,TRL,TRU){
  dd<-c(1,2)
  if(SORC$dimension >2){
    #SORC$vis<-readndPlotPer()
    SORC$vis<-T
    if(SORC$vis==F){ return(SORC)
    }else{
      SORC<-visualizenD(objModel,fn,SORC,con.model,TRL,TRU,dd)
      return(SORC)
    }
  }
  
  
  
  THREE_PLOTS=T         # if TRUE, arrange three plots: EI, F, EImod
                        # if FALSE, arrange two plots:  EI, EImod
  OBJ_CON1_NEEDED=F     # if TRUE, calculate objP, con1P, con2P (they seem not needed anymore) 
  DEBUG_OLD_VERSION=F   # if TRUE, run also the old version for EI.value and EImod.value 
                        # (slower by factor >100) and compare
  IDEAL_F=T             # debug plot of true F, can be very slow (additional 3 sec)
  SHOW_GRAD=F           # Show several gradient arrows
  PLOT_CON=F            # plot constraint borders for G6  problem 
  PLOT_SD=F             # plot sigma function on the search space
  
  
  num<-nrow(SORC$pop)
  #visZoom=F
  if(SORC$visZoom){
    dev<-SORC$iter
    pow<-2.1 # 1.8 2 
    delta<-max(0.000002,0.1/dev^pow)
    TRL<-SORC$solu-delta          # /WK/ limit the zoom-in
    TRU<-SORC$solu+delta  

     TRL<-pmin(TRL,SORC$bestSolu)
     TRU<-pmax(TRU,SORC$bestSolu)
  }
  
  
   x1<-seq(from=TRL[1],to=TRU[1],length=100)
   x2<-seq(from=TRL[2],to=TRU[2],length=100)
   xMat<-expand.grid(x=x1,y=x2)
   
  #deprectaeds 
  # if(is.null(SORC$probP)){
  #   z<-apply(xMat,1,FUN=function(x){
  #     y<-SORC$fn(x)
  #     return(y[1])
  #   })
  #   probMat<-data.frame(xMat,z=z)
  #   
  #   probP<-ggplot(data=probMat,aes(x=x,y=y,z=z))
  #   probP<-probP + geom_tile(aes(fill = z)) + stat_contour(color="black")
  #   if(!is.null(SORC$solu)){
  #     if(is.matrix(SORC$solu)) {soluMat<-data.frame(SORC$solu[,1],SORC$solu[,2],SORC$optim)
  #     }else{
  #       soluMat<-data.frame(SORC$solu[1],SORC$solu[2],SORC$optim)  
  #     }
  #    # browser()
  #     names(soluMat)<-c("x","y","z")
  #   }
  #   probP<-probP +geom_point(data=soluMat,aes(x=x,y=y),color="yellow",size=4)
  # } ##End of probP
  
   xMatP<-xMat
   names(xMatP)<-c("x1","x2")
 
   if (OBJ_CON1_NEEDED)   {
     obj<-predict(objModel,newdata=xMatP,type=SORC$kType)
     con1<-predict(con.model[[1]],newdata=xMatP,type=SORC$kType)
     if(SORC$nConstraints>1) con2<-predict(con.model[[2]],newdata=xMatP,type=SORC$kType)
     
     
     objD<-data.frame(xMat,av=obj$mean,sdm=obj$sd)
     con1D<-data.frame(xMat,av=con1$mean,sdm=con1$sd)
     if(SORC$nConstraints>1)  con2D<-data.frame(xMat,av=con2$mean,sdm=con2$sd)
     
     objP<-ggplot(data=objD,aes(x=x,y=y,z=av))
     objP<-objP + geom_tile(aes(fill = av)) + stat_contour(color="black")
   #  objP<-objP +coord_cartesian(xlim=c(SORC$visL,SORC$visU),ylim=c(SORC$visL,SORC$visU))
     objP<-objP +coord_cartesian(xlim=c(SORC$TRL[1],SORC$TRU[1]),ylim=c(SORC$TRL[2],SORC$TRU[2]))
     
     # objP<-objP + geom_point(data=soluMat,aes(x=x,y=y),color="yellow",size=4)
     
     
     con1P<-ggplot(data=con1D,aes(x=x,y=y,z=av))
     con1P<-con1P + geom_tile(aes(fill = av)) + stat_contour(color="black")+ scale_fill_gradientn(colours=brewer.pal(6,"YlOrRd"))
     # con1P<-con1P + geom_point(data=soluMat,aes(x=x,y=y),color="yellow",size=4)
     con1P<-con1P +coord_cartesian(xlim=c(SORC$TRL[1],SORC$TRU[1]),ylim=c(SORC$TRL[2],SORC$TRU[2]))
     
     
     if(SORC$nConstraints>1){
       con2P<-ggplot(data=con2D,aes(x=x,y=y,z=av))
       con2P<-con2P + geom_tile(aes(fill = av)) + stat_contour(color="black") + scale_fill_gradientn(colours=brewer.pal(6,"YlOrRd"))
       #con2P<-con2P + geom_point(data=soluMat,aes(x=x,y=y),color="yellow",size=4)
       con2P<-con2P +coord_cartesian(xlim=c(SORC$TRL[1],SORC$TRU[1]),ylim=c(SORC$TRL[2],SORC$TRU[2]))
     }
     
   } # if (OBJ_CON1_NEEDED)
   
   eivT<-proc.time()
   EI.value<-myEI(xMat,objModel,type=SORC$kType,plugin=SORC$plugin[length(SORC$plugin)])
   if (DEBUG_OLD_VERSION) {
     if((SORC$verbose>3) && ((SORC$iter-SORC$initSize)%%SORC$printIter==0))print(paste("EI.value took:",round((proc.time()-eivT)[3],digits=3)))
     eivT<-proc.time()
                                           #/WK/
     EI.value_OLD<-apply(as.matrix(xMat),1,myEI,objModel,type=SORC$kType,plugin=SORC$plugin[length(SORC$plugin)])
     if((SORC$verbose>3) && ((SORC$iter-SORC$initSize)%%SORC$printIter==0))print(paste("EI.value_OLD took:",round((proc.time()-eivT)[3],digits=3)))
     if (any(EI.value!=EI.value_OLD)) {
       browser()    # we should normally not get here 
     }
   }
  # FI.value<-apply(as.matrix(xMat),1,CalculateF,con.model,SORC)
   eivT<-proc.time()
   EI.mod<-EImod(xMat,con.model,objModel,SORC)
   if (DEBUG_OLD_VERSION) {
     if((SORC$verbose>3) && ((SORC$iter-SORC$initSize)%%SORC$printIter==0)) print(paste("EI.mod took:",round((proc.time()-eivT)[3],digits=3)))
     eivT<-proc.time()
     EI.mod_OLD<-apply(as.matrix(xMat),1,EImod,con.model,objModel,SORC)
     if((SORC$verbose>3) && ((SORC$iter-SORC$initSize)%%SORC$printIter==0))   print(paste("EI.mod_OLD took:",round((proc.time()-eivT)[3],digits=3)))
     if (any(EI.mod!=EI.mod_OLD)) {
       browser()    # we should normally not get here 
     }
   }
   

   F.value <- CalculateF(xMat,con.model,SORC)$res 
   if (IDEAL_F) {
     # this can take about 3 sec computation time
     mean.sd <- CalculateF(xMat,con.model,SORC)$mean.sd     
     if((SORC$verbose>3) && ((SORC$iter-SORC$initSize)%%SORC$printIter==0))cat(sprintf("mean.sd : %s\n",paste(mean.sd,collapse=" ")))
     F.ideal <- CalculateF_ideal(xMat,ideal.sd=mean.sd,SORC)$res
   }
     
   bestsolu<-data.frame(SORC$pop[SORC$bIndex,],z=NA)
   names(bestsolu)<-c("x","y","z")
   
   EIP<-ggplot(data=data.frame(xMat,z=EI.value),aes(x=x,y=y,z=z))
   EIP<-EIP+geom_tile(aes(fill=z))+stat_contour(color="black")
   EIP<-EIP+ scale_fill_gradientn(name="EI(obj)",colours=brewer.pal(7,"Greens"))
   EIP<-EIP+geom_point(data=data.frame(SORC$pop,z=NA),aes(x=x1,y=x2))
   EIP<-EIP+coord_cartesian(xlim=c(TRL[1],TRU[1]),ylim=c(TRL[2],TRU[2]))
   if (SHOW_GRAD) {
     arrowLength <- sqrt(sum((TRU-TRL)^2))/5
     df <- NULL
     for (f in c(0,1)) {
       xif=SORC$infill+f*arrowLength*c(-1,1)/2
       GradE<-grad(func=myEI,x=xif,method="simple",method.args=list(eps=1e-6),
                   objModel=objModel,plugin=SORC$plugin[length(SORC$plugin)],type=SORC$kType)
       egrad <- GradE/sqrt(sum(GradE^2))
       df <- rbind(df,data.frame(x=xif[1],xend=xif[1]+egrad[1]*arrowLength,
                                 y=xif[2],yend=xif[2]+egrad[2]*arrowLength,z=NA))
     }
     EIP<-EIP+geom_segment(data=df,aes(x=x,y=y,xend=xend,yend=yend,z=z),
                           color="blue",arrow = arrow(length = unit(0.08, "npc")))
   }
   
   FIP<-ggplot(data=data.frame(xMat,z=F.value),aes(x=x,y=y,z=z))
   FIP<-FIP+geom_tile(aes(fill=z))+stat_contour(color="black")
   FIP<-FIP+ scale_fill_gradientn(name="F(con)",colours=brewer.pal(6,"YlOrRd"))
   FIP<-FIP+ geom_point(data=soluMat,aes(x=x,y=y),fill="yellow",size=4,shape=22)
   FIP<-FIP+ geom_point(data=bestsolu,aes(x=x,y=y),color="blue",size=4)
   FIP<-FIP+geom_point(data=data.frame(SORC$pop,z=NA),aes(x=x1,y=x2))
   FIP<-FIP+coord_cartesian(xlim=c(TRL[1],TRU[1]),ylim=c(TRL[2],TRU[2]))
   if (SHOW_GRAD) {
     arrowLength <- sqrt(sum((TRU-TRL)^2))/5
     df <- NULL
     for (f in c(0,0.05,1)) {
       xif=SORC$infill+f*arrowLength*c(-1,1)/2
       GradF<-grad(func=function(x){CalculateF(x,con.model,SORC)$res},x=xif,method="simple",
                   method.args=list(eps=1e-6))
       fgrad <- GradF/sqrt(sum(GradF^2))
       df <- rbind(df,data.frame(x=xif[1],xend=xif[1]+fgrad[1]*arrowLength,
                                 y=xif[2],yend=xif[2]+fgrad[2]*arrowLength,z=NA))
     }
     FIP<-FIP+geom_segment(data=df,aes(x=x,y=y,xend=xend,yend=yend,z=z),
                           color="blue",arrow = arrow(length = unit(0.08, "npc")))
   }
   
   if (IDEAL_F) {
     F2P<-ggplot(data=data.frame(xMat,z=F.ideal),aes(x=x,y=y,z=z))
     F2P<-F2P+geom_tile(aes(fill=z))+stat_contour(color="black")
     F2P<-F2P+ scale_fill_gradientn(name="F(ideal)",colours=brewer.pal(6,"YlOrRd"))
     F2P<-F2P+ geom_point(data=soluMat,aes(x=x,y=y),fill="yellow",size=4,shape=22)
     F2P<-F2P+ geom_point(data=bestsolu,aes(x=x,y=y),color="blue",size=4)
     F2P<-F2P+geom_point(data=data.frame(SORC$pop,z=NA),aes(x=x1,y=x2))
     F2P<-F2P+coord_cartesian(xlim=c(TRL[1],TRU[1]),ylim=c(TRL[2],TRU[2]))
   } else {
     F2P<-ggplot(data=data.frame(xMat,z=F.value),aes(x=x,y=y,z=z))
     F2P<-F2P+geom_tile(aes(fill=z))+stat_contour(color="black")
     F2P<-F2P+ scale_fill_gradientn(name="F(con)",colours=brewer.pal(6,"YlOrRd"))
     F2P<-F2P+ geom_point(data=soluMat,aes(x=x,y=y),fill="yellow",size=4,shape=22)
     F2P<-F2P+ geom_point(data=bestsolu,aes(x=x,y=y),color="blue",size=4)
     F2P<-F2P+geom_point(data=data.frame(SORC$pop,z=NA),aes(x=x1,y=x2))
     F2P<-F2P+coord_cartesian(xlim=c(TRL[1],TRU[1]),ylim=c(TRL[2],TRU[2]))
   }
   
   #    FIP<-ggplot(data=data.frame(xMat,z=FI.value),aes(x=x,y=y,z=z))
#    FIP<-FIP+geom_tile(aes(fill=z))+stat_contour(color="black") + scale_fill_gradientn(colours=brewer.pal(6,"YlOrRd"))
#    
   
   SORC$EImax<-max(-EI.mod)
   EImodP<-ggplot(data=data.frame(xMat,z=-EI.mod),aes(x=x,y=y,z=z))
   EImodP<-EImodP+geom_tile(aes(fill=z))+stat_contour(color="black") 
   EImodP<-EImodP+ scale_fill_gradientn(name=expression(EI[mod]),colours=brewer.pal(6,"YlOrRd"))
   EImodP<-EImodP+ geom_point(data=soluMat,aes(x=x,y=y),fill="yellow",size=4,shape=22)
   EImodP<-EImodP+ geom_point(data=bestsolu,aes(x=x,y=y),color="blue",size=4)
   EImodP<-EImodP+coord_cartesian(xlim=c(TRL [1],TRU[1]),ylim=c(TRL[2],TRU[2]))
   
  if(!PLOT_CON) EImodP<-EImodP+geom_point(data=data.frame(SORC$pop,z=NA),aes(x=x1,y=x2))
   #EImodP<-EImodP+coord_cartesian(xlim=c(SORC$visL[1],SORC$visU[1]),ylim=c(SORC$visL[2],SORC$visU[2]))
   if(PLOT_CON){
     myx<-x1[which(x1>=SORC$solu[1])]
     rx1<-rescale(myx,to=c(13,100),from=c(-1,1))
     con1<--sqrt(100-(rx1-5)^2)+5
     con2<--sqrt(82.81-(rx1-6)^2)+5
     con1<-rescale(con1,to=c(-1,1),from=c(0,100))
     con2<-rescale(con2,to=c(-1,1),from=c(0,100))
     
     FRAME1<-data.frame(x=myx,y=con1,z=NA)
     FRAME2<-data.frame(x=myx,y=con2,z=NA)
     #browser()
     
     EImodP<-EImodP+geom_line(data=FRAME1,aes(x=x,y=y),color="black",size=1)
     EImodP<-EImodP+geom_line(data=FRAME2,aes(x=x,y=y),color="black",size=1)
     EImodP<-EImodP+ geom_point(data=soluMat,aes(x=x,y=y),fill="blue",size=6,shape=22)
     EImodP<-EImodP+ geom_point(data=bestsolu,aes(x=x,y=y),color="black",fill="yellow",size=3.5,shape=21)
     EImodP<-EImodP + xlab(expression(x[1]))
     EImodP<-EImodP + ylab(expression(x[2]))
     
   }
   
   if (SHOW_GRAD) {
     arrowLength <- sqrt(sum((TRU-TRL)^2))/5
     df <- NULL
     for (f in c(0,0.01,1)) {
       xif=SORC$infill+f*arrowLength*c(-1,1)/2
       GradEM<- grad(func=function(x){-EImod(x,con.model,objModel,SORC)},x=xif,
                     method="simple",method.args=list(eps=1e-6))
       emgrad <- GradEM/sqrt(sum(GradEM^2))
# --- only an intermediate test: (c1+c2) points along the red rim ---       
#        conGrad <- SORC$conGrad
#        c1=conGrad[1,]/SORC$gGrad[1]
#        c2=conGrad[2,]/SORC$gGrad[2]
#        cgrad=+(c1+c2)
#        emgrad=cgrad/sqrt(sum(cgrad^2))
       df <- rbind(df,data.frame(x=xif[1],xend=xif[1]+emgrad[1]*arrowLength,
                                 y=xif[2],yend=xif[2]+emgrad[2]*arrowLength,z=NA))
     }
     EImodP<-EImodP+geom_segment(data=df,aes(x=x,y=y,xend=xend,yend=yend,z=z),
                                 color="blue",arrow = arrow(length = unit(0.08, "npc")))
   }
     
   annText<-paste("Best Sol [",SORC$iter,"]: ", round(SORC$Fres[SORC$bIndex],digits=2)," | ",round(SORC$DF$maxViol[SORC$bIndex],digits=3),sep="")
   
   annText0<-paste("                     obj","     ", "  maxViol \n",annText)
   if(any(SORC$NU!=1)){  
     #browser()
     annText0<-paste("    adjusting nu     obj","     ", "  maxViol \n",annText)}

   #browser()
   
   lly<-TRU[2]-TRL[2]
   llx<-TRU[1]-TRL[1]
   
   textFrame<-data.frame(x=TRL[1]+1.5*llx/3,y=TRU[2]-(lly/13),z=annText0)
   EImodP<-EImodP+geom_label(data=textFrame,aes(x=x,y=y,label=z),alpha=0.5)
   #browser()
  # EImodP<-EImodP+geom_rect(xmin=TRL[1],xmax=TRU[1],ymin=TRL[2],ymax=TRU[2],alpha=0.01)
   
   if (THREE_PLOTS) {
     EIFIP<-arrangeGrob(EIP,FIP,nrow=2)         # arrangeGrob does not draw on current device
     #if (IDEAL_F) {
       EIFI2 <- arrangeGrob(EImodP,F2P,nrow=2) 
     #} else {
     #   EIFI2 <- EImodP
     #}
     newp<-grid.arrange(EIFIP,EIFI2, ncol = 2) # grid.arrange draws on current device
   } else {
     newp<-grid.arrange(EIP,EImodP, ncol = 2)
   }
   
   indMax = which(-EI.mod==max(-EI.mod))
   s <- "\n"
   if (length(indMax)>1) s <- sprintf(" and at %d other points\n",length(indMax));
   if((SORC$verbose>3) && ((SORC$iter-SORC$initSize)%%SORC$printIter==0)) cat(sprintf(paste0("EImod at infill point (%5.3f,%5.3f): %8.6f,\n",
                      "max(EImod) on grid at (%5.3f,%5.3f): %8.6f %s"),
               SORC$infill[1],SORC$infill[2],SORC$EI, 
               xMat[indMax[1],1],xMat[indMax[1],2],max(-EI.mod),s));
   #browser()
   
   
   if(PLOT_SD){
     obj<-predict(objModel,newdata=xMatP,type=SORC$kType)
     mySD<-obj$sd
     SDPLOT<-ggplot(data=data.frame(xMat,z=mySD),aes(x=x,y=y,z=z))
     SDPLOT<-SDPLOT+geom_tile(aes(fill=z))+stat_contour(color="black")
     SDPLOT<-SDPLOT+ scale_fill_gradientn(name="SD(obj)",colours=brewer.pal(7,"Blues"))
     SDPLOT<-SDPLOT+geom_point(data=data.frame(SORC$pop,z=NA),aes(x=x1,y=x2))
     SDPLOT<-SDPLOT+coord_cartesian(xlim=c(TRL[1],TRU[1]),ylim=c(TRL[2],TRU[2]))
     
   }
   
   #ggsave(paste("../Figures/",SORC$name,"-",SORC$fCalc,"-EI",num,".pdf",sep=""),newp,width=12,height=6)
   ggsave(paste("../Figures/",SORC$name,"-EI",num,".pdf",sep=""),newp,width=8.5,height=6)
   if(PLOT_CON) ggsave(paste("../Figures/",SORC$name,"-EImodCon",num,".pdf",sep=""),EImodP,width=5,height=4)
   if(PLOT_SD)  ggsave(paste("../Figures/",SORC$name,"-SD",num,".pdf",sep=""),SDPLOT,width=5,height=4)
   # width=8.5 for a fixed axis ratio (only then the angles between height lines and gradient arrows 
      # are o.k. /WK/)
   
#   plotFrame<-data.frame(x=t,obj=obj$mean,
#                         con=con$mean,
#                         EI=EI.value,
#                         FI=FI.value,
#                         EImod=EI.mod,
#                         objl=obj$lower95,
#                         obju=obj$upper95,
#                         conl=con$lower95,
#                         conu=con$upper95)
#   pointFrame<-data.frame(x=as.numeric(SORC$pop[,1]),
#                          obj=SORC$Fres,
#                          con=SORC$Gres[,1],
#                          EI=NA,
#                          FI=NA,
#                          EImod=NA,
#                          objl=NA,
#                          obju=NA,
#                          conl=NA,
#                          conu=NA)
# #  plot(t,obj$mean,main="obj model",type="l")
# #  points(as.numeric(SORC$pop[,1]),SORC$Fres)
# #  plot(t,con$mean,main="con model",type="l")
# #  plot(t,EI.value,main="EI",type="l")
# #  plot(t,FI.value,main="FI",type="l")
# #  plot(t,-EI.mod,main="EI-mod",type="l")
#   
#   p<-ggplot(data=plotFrame)
#   p<-p+geom_line(size=2,color="black",aes(x=t,y=obj))
#   p<-p+geom_line(size=1,color="black",linetype="dashed",aes(x=t,y=objl))
#   p<-p+geom_line(size=1,color="black",linetype="dashed",aes(x=t,y=obju))
#   p<-p+geom_point(data=pointFrame,aes(x=x,y=obj),size=5,color="black")
#   p<-p+geom_point(data=pointFrame,aes(x=x,y=obj),size=3,color="white")
#   p<-p+labs(title="Objective Model",x="x")
#   
#   
#   ggsave(paste("../Figures/objm",num,".pdf",sep=""), width = 8, height = 6)
#   p1<-p
#   
#   p<-ggplot(data=plotFrame)
#   p<-p+geom_line(size=2,color="blue",aes(x=t,y=con))
#   p<-p+geom_line(size=1,color="blue",linetype="dashed",aes(x=t,y=conl))
#   p<-p+geom_line(size=1,color="blue",linetype="dashed",aes(x=t,y=conu))
#   
#   p<-p+geom_point(data=pointFrame,aes(x=x,y=con),size=5,color="blue")
#   p<-p+geom_point(data=pointFrame,aes(x=x,y=con),size=3,color="white")
#   p<-p+labs(title="Constraint Model",x="x")
#   
#   p<-p+geom_hline(yintercept=0)
#   ggsave(paste("../Figures/conm",num,".pdf",sep=""), width = 8, height = 6)
#   p2<-p
#   
#   
#   p<-ggplot(data=plotFrame)
#   p<-p+geom_line(aes(x=t,y=EI),size=2)
#   #p<-p+geom_vline(xintercept = pointFrame$x,linetype="dashed",color="red",size=0.5)
#   p<-p+labs(title="EI(x) (Objetcive function)",x="x")
#   ggsave(paste("../Figures/ei",SORC$name,num,".pdf",sep=""), width = 8, height = 6)
#   p3<-p
#   
#   p<-ggplot(data=plotFrame)
#   p<-p+geom_line(aes(x=t,y=FI),size=2)
#   #p<-p+geom_vline(xintercept = pointFrame$x,linetype="dashed",color="red",size=0.5)
#   p<-p+labs(title="F(x) (Feasibility Measure)",x="x")
#   ggsave(paste("../Figures/fi",SORC$name,num,".pdf",sep=""), width = 8, height = 6)
#   p4<-p
#   
#   p<-ggplot(data=plotFrame)
#   p<-p+geom_line(aes(x=t,y=-EImod),size=2)
#  # p<-p+geom_vline(xintercept = pointFrame$x,linetype="dashed",color="red",size=0.5)
#   p<-p+labs(title="EImod(x) (Constrained EI)",x="x")
#   p5<-p
#   ggsave(paste("../Figures/eimod",SORC$name,num,".pdf",sep=""), width = 8, height = 6)
#   
#   newp<-grid.arrange(p1, p3,p2,p5, ncol = 2)
#  # browser()
#   ggsave(paste(SORC$name,num,".pdf",sep=""),newp)
#   
 #browser()
  
  
return(SORC)
}