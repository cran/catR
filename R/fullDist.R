fullDist <- function(th, it, method = "BM", priorDist="norm",
                     priorPar=c(0,1), weight = "Huber", 
                     tuCo   = 1, range  = c(-4 ,4), parInt = c(-4, 4, 33)){
### generation of all required binary response patterns 
dataGen <- function(n, model="1PL"){
if (model=="1PL"){
  res <- matrix(0, n + 1, n)
  for(i in 1:n) res[i + 1, 1:i] <- 1
}
else{
  res   <- matrix(NA, 2^n, n)
  for (i in 1:n) res[, i] <- gl(2, 2^(n - i), 2^n)
  res <- res - 1
}
  return(res)
}
### Lord-Wingersky algorithm
LW<-function(th,it,D=1){
P<-Pi(th,it,D=D)$Pi
Q<-1-P
res<-matrix(NA,nrow(it)+1,nrow(it))
res[1,1]<-Q[1]
res[2,1]<-P[1]
for (i in 2:nrow(it)){
for (j in 0:i){
if (j==0) res[j+1,i]<-res[j+1,i-1]*Q[i]
else{
if (j==i) res[j+1,i]<-res[j,i-1]*P[i]
else res[j+1,i]<-res[j,i-1]*P[i]+res[j+1,i-1]*Q[i]
}
}}
RES<-cbind(0:ncol(res),res[,ncol(res)])
return(RES)}
# main function
it<-rbind(it)
if (abs(mean(it[,1])-1)<1e-5 & var(it[,1])<1e-5) mod<-"1PL"
else mod<-"other"
data <- dataGen(nrow(it),model=mod)
if (mod=="1PL"){
res <- matrix(NA, nrow(data), 1 + length(th))
  for(i in 1:nrow(data)){
    res[i, 1] <- thetaEst(it, data[i, ], method = method, 
                          priorDist = priorDist,
                          priorPar = priorPar, weight = weight,
                          tuCo = tuCo, range = range)
  }
  for(j in 1:length(th)){
    res[, j + 1] <- LW(th[j], it)[, 2]
  }
}
else{
  res <- matrix(NA, nrow(data), 1 + length(th))
  for (i in 1:nrow(data)){
    res[i, 1] <- thetaEst(it, data[i, ], method = method,
                          priorDist = priorDist,
                          priorPar = priorPar, weight = weight,
                          tuCo = tuCo, range = range)
    for (j in 1:length(th)){
      pi <- Pi(th[j], it)$Pi
      qi <- 1 - pi
      res[i, 1 + j] <- prod(pi^data[i, ] * qi^(1 - data[i, ]))
    }
  }
}
  return(res)
}

