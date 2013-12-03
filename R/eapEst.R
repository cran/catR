eapEst<-function(it,x,D=1,priorDist="norm",priorPar=c(0,1),lower=-4,upper=4,nqp=33){
L<-function(th,it,x) prod(Pi(th,it,D=D)$Pi^x*(1-Pi(th,it,D=D)$Pi)^(1-x))
g<-function(s){
res<-NULL
for (i in 1:length(s)) res[i]<-switch(priorDist,
norm=s[i]*dnorm(s[i],priorPar[1],priorPar[2])*L(s[i],it,x),
unif=s[i]*dunif(s[i],priorPar[1],priorPar[2])*L(s[i],it,x),
Jeffreys=s[i]*sqrt(sum(Ii(s[i],it,D=D)$Ii))*L(s[i],it,x))
return(res)}
h<-function(s){
res<-NULL
for (i in 1:length(s)) res[i]<-switch(priorDist,
norm=dnorm(s[i],priorPar[1],priorPar[2])*L(s[i],it,x),
unif=dunif(s[i],priorPar[1],priorPar[2])*L(s[i],it,x),
Jeffreys=sqrt(sum(Ii(s[i],it,D=D)$Ii))*L(s[i],it,x))
return(res)}
X<-seq(from=lower,to=upper,length=nqp)
Y1<-g(X)
Y2<-h(X)
RES<-integrate.xy(X,Y1)/integrate.xy(X,Y2)
return(RES)}


