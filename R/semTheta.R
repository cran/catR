semTheta<-function(thEst,it,x=NULL,D=1,method="BM",priorDist="norm",priorPar=c(0,1),eapPar=c(-4,4,20)){
if (method=="EAP"){
RES<-eapSem(thEst,it,x=x,D=D,priorDist=priorDist,priorPar=priorPar,lower=eapPar[1],upper=eapPar[2],nqp=eapPar[3])
}
else{
## function J(th) and first derivative
Ji<-function(th,it,D=1){
pr<-Pi(th,it,D=D)
P<-pr$Pi
Q<-1-P
dP<-pr$dPi
d2P<-pr$d2Pi
d3P<-pr$d3Pi
Ji<-dP*d2P/(P*Q)
dJi<-(P*Q*(d2P^2+dP*d3P)-dP^2*d2P*(Q-P))/(P^2*Q^2)
res<-list(Ji=Ji,dJi=dJi)
return(res)}
## Derivative of r0(th) for standard error
dr0<-function(th,it,D=1,method="BM",priorDist="norm",priorPar=c(0,1)){
if (method=="BM") res<-switch(priorDist,
norm=-1/priorPar[2]^2,
unif=0,
Jeffreys=(sum(Ii(th,it,D=D)$d2Ii)*sum(Ii(th,it,D=D)$Ii)-sum(Ii(th,it,D=D)$dIi)^2)/(2*sum(Ii(th,it,D=D)$Ii)^2))
else res<-switch(method,
ML=0,
WL=(sum(Ji(th,it,D=D)$dJi)*sum(Ii(th,it,D=D)$Ii)-sum(Ji(th,it,D=D)$Ji)*sum(Ii(th,it,D=D)$dIi))/(2*sum(Ii(th,it,D=D)$Ii)^2))
return(res)}
info<-sum(Ii(thEst,it,D=D)$Ii)
res<--dr0(thEst,it,D=D,method=method,priorDist=priorDist,priorPar=priorPar)+info
RES<-1/sqrt(res)
}
return(RES)}
