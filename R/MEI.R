
MEI<-function (itemBank, item, x, theta, it, method="BM", 
     priorDist = "norm", priorPar = c(0, 1),D=1, range=c(-4,4),parInt=c(-4,4,33), infoType="observed") 
{
if (infoType!="Fisher" & infoType!="observed") stop("'infoType' must be either 'Fisher' or 'observed'",call.=FALSE)
   th<-theta
   itj<-rbind(it,itemBank$itemPar[item,])
   th0<-thetaEst(itj,c(x,0),D=D,method=method,priorDist=priorDist,priorPar=priorPar,range=range,parInt=parInt)
   th1<-thetaEst(itj,c(x,1),D=D,method=method,priorDist=priorDist,priorPar=priorPar,range=range,parInt=parInt)
   p1<-Pi(th,rbind(itemBank$itemPar[item,]),D=D)$Pi
   p0<-1-p1
if (infoType=="Fisher"){
Ij0<-Ii(th0,rbind(itemBank$itemPar[item,]),D=D)$Ii
Ij1<-Ii(th1,rbind(itemBank$itemPar[item,]),D=D)$Ii
}
else{
   Ij0<-OIi(th0,rbind(itemBank$itemPar[item,]),0,D=D)
   Ij1<-OIi(th0,rbind(itemBank$itemPar[item,]),1,D=D)
}
res<-p0*Ij0+p1*Ij1
return(as.numeric(res))}