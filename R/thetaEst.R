thetaEst<-function(it,x,D=1,method="BM",priorDist="norm",priorPar=c(0,1),range=c(-4,4),parInt=c(-4,4,33)){
if (method=="EAP"){
res<-eapEst(it,x,D=D,priorDist=priorDist,priorPar=priorPar,lower=parInt[1],upper=parInt[2],nqp=parInt[3])
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
## Function r0(th) setting the derivative of the log-"prior" distribution
r0<-function(th,it,D=1,method="BM",priorDist="norm",priorPar=c(0,1)){
if (method=="BM") res<-switch(priorDist,
norm=(priorPar[1]-th)/priorPar[2]^2,
unif=0,
Jeffreys=sum(Ii(th,it,D=D)$dIi)/(2*sum(Ii(th,it,D=D)$Ii)))
else res<-switch(method,
ML=0,
WL=sum(Ji(th,it,D=D)$Ji)/(2*sum(Ii(th,it,D=D)$Ii)))
return(res)}
## Function r(th,x) setting the drivative of the log-likelihood
r<-function(th,it,x,D=1){
pr<-Pi(th,it,D=D)
P<-pr$Pi
Q<-1-P
dP<-pr$dPi
res<-sum(dP*(x-P)/(P*Q))
return(res)}
## Function T(th) yielding the ability estimate 
T<-function(th,it,x,D=1,method="BM",priorDist="norm",priorPar=c(0,1)){
r0(th,it,D=D,method=method,priorDist=priorDist,priorPar=priorPar)+
r(th,it,x,D=D)
}
if (method=="BM" & priorDist=="unif") f<-function(th) T(th,it,x,D=D,method="ML")
else f<-function(th) T(th,it,x,D=D,method=method,priorDist=priorDist,priorPar=priorPar)
if (method=="BM" & priorDist=="unif") RANGE<-priorPar
else RANGE<-range
if (f(RANGE[1])<0 & f(RANGE[2])<0) res<-RANGE[1]
else{
if (f(RANGE[1])>0 & f(RANGE[2])>0) res<-RANGE[2]
else res<-uniroot(f,RANGE)$root
}
}
return(res)}

