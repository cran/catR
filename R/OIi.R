OIi<-function(th,it,x,D=1){
pr<-Pi(th,it,D=D)
P<-pr$Pi
Q<-1-P
dP<-pr$dPi
d2P<-pr$d2Pi
res<-(P*Q*dP^2-(x-P)*(P*Q*d2P+dP^2*(P-Q)))/(P^2*Q^2)
return(res)}

