nextItem<-function(itemBank,theta,out=NULL,method="info"){
if (method=="info"){
ind<-(1:length(itemBank$theta))[abs(itemBank$theta-theta)==min(abs(itemBank$theta-theta))]
info<-itemBank$infoTab[ind,]
items<-rep(1,nrow(itemBank$itemPar))
items[out]<-0
if (length(ind)==1) {
keep<-(1:length(info))[items==1]
select<-min((1:length(info[keep]))[info[keep]==max(info[keep])])
res<-list(item=keep[select],par=itemBank$itemPar[keep[select],],info=max(info[keep]),method=method)
}
else {
keep<-(1:ncol(info))[items==1]
pr<-NULL
for (j in 1:length(ind)){
select<-(1:length(info[j,keep]))[info[j,keep]==max(info[j,keep])]
pr<-rbind(pr,c(keep[select],info[j,keep[select]]))
}
ind2<-min((1:nrow(pr))[pr[,2]==max(pr[,2])])
res<-list(item=pr[ind2,1],par=itemBank$itemPar[pr[ind2,1],],info=max(pr[,2]),method=method)
}
}
else{
if (method=="Owen"){
indic<-rep(1,nrow(itemBank$itemPar))
indic[out]<-0
items<-(1:nrow(itemBank$itemPar))[indic==1]
ind<-items[abs(itemBank$itemPar[items,2]-theta)==min(abs(itemBank$itemPar[items,2]-theta))]
if (length(ind)>1) ind<-min(ind)
th<-(1:length(itemBank$theta))[abs(itemBank$theta-theta)==min(abs(itemBank$theta-theta))]
if (length(th)>1) th<-min(th)
res<-list(item=ind,par=itemBank$itemPar[ind,],info=itemBank$infoTab[th,ind],method=method)
}
else stop("'method' must be either 'info' or 'Owen'",call.=FALSE)
}
return(res)}

