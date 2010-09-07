EPV<-function (itemBank, item, x, theta, it, priorDist = "norm", 
    priorPar = c(0, 1), D = 1, parInt = c(-4,4, 33)) {

th <- theta
itj <- rbind(it, itemBank$itemPar[item, ])
p1 <- Pi(th, rbind(itemBank$itemPar[item, ]), D = D)$Pi
p0 <- 1 - p1

th0<-eapEst(itj,c(x,0),D=D,priorDist=priorDist,priorPar=priorPar,lower=parInt[1],upper=parInt[2],nqp=parInt[3])
th1<-eapEst(itj,c(x,1),D=D,priorDist=priorDist,priorPar=priorPar,lower=parInt[1],upper=parInt[2],nqp=parInt[3])

var0<-(eapSem(th0,itj,c(x,0),D=D,priorDist=priorDist,priorPar=priorPar,lower=parInt[1],upper=parInt[2],nqp=parInt[3]))^2
var1<-(eapSem(th1,itj,c(x,1),D=D,priorDist=priorDist,priorPar=priorPar,lower=parInt[1],upper=parInt[2],nqp=parInt[3]))^2
res<-p0*var0+p1*var1
return(as.numeric(res))}