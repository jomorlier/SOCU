visualizeOne<-function(objModel,fn,SORC,con.model){
  num<-nrow(SORC$pop)
  t<-seq(from=SORC$dLower,to=SORC$dUpper,length=300)
  y<-sapply(t,SORC$fn)
  names(t)<-rep("x1",300)
  obj<-predict(objModel,newdata=t,type=SORC$kType)
  con<-predict(con.model[[1]],newdata=t,type=SORC$kType)
 # EI.value<-apply(as.matrix(t),1,EI,objModel,type=SORC$kType)
  
  EI.value<-apply(as.matrix(t),1,EI,objModel,type=SORC$kType,plugin=SORC$plugin)
  #browser()
  FI.value<-apply(as.matrix(t),1,CalculateF,con.model,SORC)
  EI.mod<-apply(as.matrix(t),1,EImod,con.model,objModel,SORC)
  plotFrame<-data.frame(x=t,obj=obj$mean,
                        con=con$mean,
                        EI=EI.value,
                        FI=FI.value,
                        EImod=EI.mod,
                        objl=obj$lower95,
                        obju=obj$upper95,
                        conl=con$lower95,
                        conu=con$upper95)
  pointFrame<-data.frame(x=as.numeric(SORC$pop[,1]),
                         obj=SORC$Fres,
                         con=SORC$Gres[,1],
                         EI=NA,
                         FI=NA,
                         EImod=NA,
                         objl=NA,
                         obju=NA,
                         conl=NA,
                         conu=NA)
#  plot(t,obj$mean,main="obj model",type="l")
#  points(as.numeric(SORC$pop[,1]),SORC$Fres)
#  plot(t,con$mean,main="con model",type="l")
#  plot(t,EI.value,main="EI",type="l")
#  plot(t,FI.value,main="FI",type="l")
#  plot(t,-EI.mod,main="EI-mod",type="l")
  
  p<-ggplot(data=plotFrame)
  p<-p+geom_line(size=2,color="black",aes(x=t,y=obj))
  p<-p+geom_line(size=1,color="black",linetype="dashed",aes(x=t,y=objl))
  p<-p+geom_line(size=1,color="black",linetype="dashed",aes(x=t,y=obju))
  p<-p+geom_point(data=pointFrame,aes(x=x,y=obj),size=5,color="black")
  p<-p+geom_point(data=pointFrame,aes(x=x,y=obj),size=3,color="white")
  p<-p+labs(title="Objective Model",x="x")
  
  
  ggsave(paste("../Figures/objm",num,".pdf",sep=""), width = 8, height = 6)
  p1<-p
  
  p<-ggplot(data=plotFrame)
  p<-p+geom_line(size=2,color="blue",aes(x=t,y=con))
  p<-p+geom_line(size=1,color="blue",linetype="dashed",aes(x=t,y=conl))
  p<-p+geom_line(size=1,color="blue",linetype="dashed",aes(x=t,y=conu))
  
  p<-p+geom_point(data=pointFrame,aes(x=x,y=con),size=5,color="blue")
  p<-p+geom_point(data=pointFrame,aes(x=x,y=con),size=3,color="white")
  p<-p+labs(title="Constraint Model",x="x")
  
  p<-p+geom_hline(yintercept=0)
  ggsave(paste("../Figures/conm",num,".pdf",sep=""), width = 8, height = 6)
  p2<-p
  
  
  p<-ggplot(data=plotFrame)
  p<-p+geom_line(aes(x=t,y=EI),size=2)
  #p<-p+geom_vline(xintercept = pointFrame$x,linetype="dashed",color="red",size=0.5)
  p<-p+labs(title="EI(x) (Objetcive function)",x="x")
  ggsave(paste("../Figures/ei",SORC$name,num,".pdf",sep=""), width = 8, height = 6)
  p3<-p
  
  p<-ggplot(data=plotFrame)
  p<-p+geom_line(aes(x=t,y=FI),size=2)
  #p<-p+geom_vline(xintercept = pointFrame$x,linetype="dashed",color="red",size=0.5)
  p<-p+labs(title="F(x) (Feasibility Measure)",x="x")
  ggsave(paste("../Figures/fi",SORC$name,num,".pdf",sep=""), width = 8, height = 6)
  p4<-p
  
  p<-ggplot(data=plotFrame)
  p<-p+geom_line(aes(x=t,y=-EImod),size=2)
 # p<-p+geom_vline(xintercept = pointFrame$x,linetype="dashed",color="red",size=0.5)
  p<-p+labs(title="EImod(x) (Constrained EI)",x="x")
  p5<-p
  ggsave(paste("../Figures/eimod",SORC$name,num,".pdf",sep=""), width = 8, height = 6)
  
  newp<-grid.arrange(p1, p3,p2,p5, ncol = 2)
 # browser()
  ggsave(paste(SORC$name,num,".pdf",sep=""),newp)
  
 # browser()
  
  

}