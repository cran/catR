startItems<-function(it,fixItems=NULL,seed=NULL,nrItems=1,theta=0,bw=4,range=c(-4,4)){
b<-it[,2]
if (is.null(fixItems)==FALSE){
items<-fixItems
par<-it[fixItems,]
bOpt<-NA
}
else{
if (is.null(seed)==FALSE){
set.seed(seed)
repeat{
items<-round(runif(nrItems,0,1)*nrow(it))+1
if (length(unique(items))==nrItems) break
}
par<-it[items,]
bOpt<-NA
}
else{
if (nrItems==1) bOpt<-theta
else {
bOpt<-seq(from=theta-bw,to=theta+bw,length=nrItems)
bOpt[1]<-max(c(bOpt[1],range[1]))
bOpt[nrItems]<-min(c(bOpt[nrItems],range[2]))
}
items=NULL
ind<-rep(1,nrow(it))
for (i in 1:length(bOpt)){
ind[items]<-0
keep<-(1:nrow(it))[ind==1]
items<-c(items,min(keep[abs(b[keep]-bOpt[i])==min(abs(b[keep]-bOpt[i]))]))
}
items=sort(items)
par=it[items,] 
}
}
res<-list(items=items,par=par,bOpt=bOpt)
return(res)
}
