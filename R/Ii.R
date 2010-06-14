Ii<-function(th,it,D=1){
pr<-Pi(th,it,D=D)
P<-pr$Pi
Q<-1-P
dP<-pr$dPi
d2P<-pr$d2Pi
d3P<-pr$d3Pi
Ii<-dP^2/(P*Q)
dIi<-dP*(2*P*Q*d2P-dP^2*(Q-P))/(P^2*Q^2)
d2Ii<-(2*P*Q*(d2P^2+dP*d3P)-2*dP^2*d2P*(Q-P))/(P^2*Q^2)-(3*P^2*Q*dP^2*d2P-P*dP^4*(2*Q-P))/(P^4*Q^2)+(3*P*Q^2*dP^2*d2P-Q*dP^4*(Q-2*P))/(P^2*Q^4)
res<-list(Ii=Ii,dIi=dIi,d2Ii=d2Ii)
return(res)}

