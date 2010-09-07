nextItem<-function (itemBank, theta=0, out = NULL, x=NULL, criterion="MFI", method = "BM", priorDist="norm", priorPar=c(0,1), D=1,range=c(-4,4),parInt=c(-4,4,33),infoType = "observed") 
{
crit<-switch(criterion,MFI="MFI",Owen="Owen",MLWI="MLWI",MPWI="MPWI",MEI="MEI",MEPV="MEPV",random="random")
if (is.null(crit)==TRUE) stop("invalid 'criterion' name",call.=FALSE)
if (crit == "MFI") {
        ind <- (1:length(itemBank$theta))[abs(itemBank$theta - 
            theta) == min(abs(itemBank$theta - theta))]
        info <- itemBank$infoTab[ind, ]
        items <- rep(1, nrow(itemBank$itemPar))
        items[out] <- 0
        if (length(ind) == 1) {
            keep <- (1:length(info))[items == 1]
            select <- min((1:length(info[keep]))[info[keep] == 
                max(info[keep])])
            res <- list(item = keep[select], par = itemBank$itemPar[keep[select], 
                ], info = max(info[keep]), criterion = criterion)
        }
        else {
            keep <- (1:ncol(info))[items == 1]
            pr <- NULL
            for (j in 1:length(ind)) {
                select <- (1:length(info[j, keep]))[info[j, keep] == 
                  max(info[j, keep])]
                pr <- rbind(pr, c(keep[select], info[j, keep[select]]))
            }
            ind2 <- min((1:nrow(pr))[pr[, 2] == max(pr[, 2])])
            res <- list(item = pr[ind2, 1], par = itemBank$itemPar[pr[ind2, 
                1], ], info = max(pr[, 2]), criterion = criterion)
        }
    }
if (crit == "Owen") {
            indic <- rep(1, nrow(itemBank$itemPar))
            indic[out] <- 0
            items <- (1:nrow(itemBank$itemPar))[indic == 1]
            ind <- items[abs(itemBank$itemPar[items, 2] - theta) == 
                min(abs(itemBank$itemPar[items, 2] - theta))]
            if (length(ind) > 1) 
                ind <- min(ind)
            th <- (1:length(itemBank$theta))[abs(itemBank$theta - 
                theta) == min(abs(itemBank$theta - theta))]
            if (length(th) > 1) 
                th <- min(th)
            res <- list(item = ind, par = itemBank$itemPar[ind, 
                ], info = itemBank$infoTab[th, ind], criterion = criterion)
          }
if (crit=="MLWI" | crit=="MPWI"){
if (length(out)==1) par<-matrix(itemBank$itemPar[out,],1,4)
else par<-itemBank$itemPar[out,]
ITEMS<-rep(1,nrow(itemBank$itemPar))
ITEMS[out]<-0
likInfo<-rep(0,nrow(itemBank$itemPar))
for (i in 1:nrow(itemBank$itemPar)){
if (ITEMS[i]==1) likInfo[i]<-MWI(itemBank,i,x,par,type=criterion,lower=parInt[1],upper=parInt[2],nqp=parInt[3],priorDist=priorDist,priorPar=priorPar)
}
ind<-min((1:length(ITEMS))[likInfo==max(likInfo)])
res<-list(item = ind, par = itemBank$itemPar[ind,], 
info = max(likInfo), criterion = criterion)
}
if (crit=="MEI"){
        items <- rep(1, nrow(itemBank$itemPar))
        items[out] <- 0
infos<-rep(0,length(items))
for (i in 1:length(items)){
if (items[i]>0) infos[i]<-MEI(itemBank,item=i,x=x,theta=theta,it=itemBank$itemPar[out,],method=method,priorDist=priorDist,priorPar=priorPar,D=D,range=range,parInt=parInt,infoType=infoType)
}
best<-min((1:length(items))[infos==max(infos)])
res<-list(item=best,par=itemBank$itemPar[best,],info=max(infos),criterion=criterion)
}

if (crit=="MEPV"){
        items <- rep(1, nrow(itemBank$itemPar))
        items[out] <- 0
epvs<-rep(1000,length(items))
for (i in 1:length(items)){
if (items[i]>0) epvs[i]<-EPV(itemBank,item=i,x=x,theta=theta,it=itemBank$itemPar[out,],priorDist=priorDist,priorPar=priorPar,D=D,parInt=parInt)
}
best<-min((1:length(items))[epvs==min(epvs)])
res<-list(item=best,par=itemBank$itemPar[best,],info=min(epvs),criterion=criterion)
}


if (crit=="random"){
        items <- rep(1, nrow(itemBank$itemPar))
        items[out] <- 0
gen<-as.integer(runif(1,0,1)*(sum(items)))+1
ind<-(1:nrow(itemBank$itemPar))[items>0][gen]
res<-list(item=ind,par=itemBank$itemPar[ind,],info=NA,criterion=criterion)
}
return(res)
}
