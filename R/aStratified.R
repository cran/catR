

aStratified<-function(itemBank,K,model=NULL){
if (!is.null(model)){
if (sum(model==c("GRM","MGRM","GPCM"))==0) stop("IRT model is not allowed for a-starified sampling",call.=FALSE)
}
if (length(K)>1) KK<-K
else {
nr<-floor(nrow(itemBank)/K)
KK<- rep(nr,K)
if (sum(KK)<nrow(itemBank)) KK[length(KK)]<-nrow(itemBank)-sum(KK[1:(length(KK)-1)])
}
ind<-cumsum(KK)
ind<-c(-1,ind)
aSort<-rank(itemBank[,1])
res<-rep(NA,nrow(itemBank))
for (i in 1:(length(ind)-1)) {
IND<-which(aSort>(ind[i]) & aSort<=ind[i+1])
res[IND]<-paste("Group",i,sep="")
}
return(res)}
