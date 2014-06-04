semTheta<-function(thEst,it,x=NULL,model=NULL,D=1,method="BM",priorDist="norm",priorPar=c(0,1),parInt=c(-4,4,33)){
if (method=="EAP"){
RES<-eapSem(thEst,it,x=x,model=model,D=D,priorDist=priorDist,priorPar=priorPar,lower=parInt[1],upper=parInt[2],nqp=parInt[3])
}
else{
if (is.null(model)){
## Derivative of r0(th) for standard error
dr0<-function(th,it,D=1,method="BM",priorDist="norm",priorPar=c(0,1)){
if (method=="BM") res<-switch(priorDist,
norm=-1/priorPar[2]^2,
unif=0,
Jeffreys=(sum(Ii(th,it,D=D)$d2Ii)*sum(Ii(th,it,D=D)$Ii)-sum(Ii(th,it,D=D)$dIi)^2)/(2*sum(Ii(th,it,D=D)$Ii)^2))
else res<-switch(method,
ML=0,
#WL=(sum(Ji(th,it,D=D)$dJi)*sum(Ii(th,it,D=D)$Ii)-sum(Ji(th,it,D=D)$Ji)*sum(Ii(th,it,D=D)$dIi))/(2*sum(Ii(th,it,D=D)$Ii)^2))
WL=0)
return(res)}
info<-sum(Ii(thEst,it,D=D)$Ii)
res<--dr0(thEst,it,D=D,method=method,priorDist=priorDist,priorPar=priorPar)+info
RES<-1/sqrt(res)
}
else{
met<-switch(method,"ML"=1,"BM"=2,"WL"=3,"EAP"=4)
pd<-switch(priorDist,"norm"=1,"unif"=2,"Jeffreys"=3)
if (met==1 | (met==2 & pd==2)) optI<-sum(Ii(thEst,it,model=model,D=D)$Ii)
if (met==2 & pd==1) optI<-sum(Ii(thEst,it,model=model,D=D)$Ii)+1/priorPar[2]^2
if ((met==2 & pd==3) | met==3){
prI<-Ii(thEst,it,model=model,D=D)
prJ<-Ji(thEst,it,model=model,D=D)
if (met==2) optI<-sum(prI$Ii)+(sum(prI$dIi)^2-sum(prI$d2Ii)*sum(prI$Ii))/(2*sum(prI$Ii)^2)
# else optI<-sum(prI$Ii)+(sum(prJ$Ji)*sum(prI$dIi)-sum(prJ$dJi)*sum(prI$Ii))/(2*sum(prI$Ii)^2)
else optI<-sum(prI$Ii)
}
RES<-1/sqrt(optI)
}
}
return(RES)}
