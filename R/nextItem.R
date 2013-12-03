nextItem<-function (itemBank, theta = 0, out = NULL, x = NULL, criterion = "MFI", 
    method = "BM", priorDist = "norm", priorPar = c(0, 1), D = 1, 
    range = c(-4, 4), parInt = c(-4, 4, 33), infoType = "observed",randomesque=1, cbControl=NULL) 
{
    crit <- switch(criterion, MFI = "MFI", Urry = "Urry", MLWI = "MLWI", 
        MPWI = "MPWI", MEI = "MEI", MEPV = "MEPV", random = "random")
    if (is.null(crit)) 
        stop("invalid 'criterion' name", call. = FALSE)
if (is.null(cbControl)) OUT<-out
else{
if (is.null(itemBank$cbGroup)) stop(paste("Item bank '", deparse(substitute(itemBank)),"' cannot be used for content balancing",sep = ""),call.=FALSE) 
if (!test.cbList(cbControl,itemBank)$test) stop(test.cbList(cbControl,itemBank)$message,call.=FALSE)
if (sum(cbControl$props)!=1) cbControl$props<-cbControl$props/sum(cbControl$props)
nrGroup<-length(cbControl$names)
if (is.null(out)) empProp<-rep(0,nrGroup)
else{
empProp<-NULL
for (i in 1:nrGroup) empProp[i]<-length(out[itemBank$cbGroup[out]==cbControl$names[i]])
empProp<-empProp/sum(empProp)
}
thProp<-cbControl$props
if (min(empProp)==0){
indGroup<-(1:nrGroup)[empProp==0]
selGroup<-ifelse(length(indGroup)==1,indGroup,sample(indGroup,1))
}
else{
indGroup<-(1:nrGroup)[(thProp-empProp)==max(thProp-empProp)]
selGroup<-ifelse(length(indGroup)==1,indGroup,sample(indGroup,1))
}
OUT<-unique(c(out,(1:length(itemBank$cbGroup))[itemBank$cbGroup!=cbControl$names[selGroup]]))
}
    if (crit == "MFI") {
ind <- (1:length(itemBank$theta))[abs(itemBank$theta - 
            theta) == min(abs(itemBank$theta - theta))]
if (length(ind)>1) ind<-min(ind)
        info <- itemBank$infoTab[ind, ]
        items <- rep(1, nrow(itemBank$itemPar))
        items[OUT] <- 0
infoVal<-sort(info[items==1],decreasing=TRUE)[min(c(randomesque,sum(items)))]
keep<-(1:length(info))[items==1 & info>=infoVal]
select<-ifelse(length(keep)==1,keep,sample(c(keep),1))
res <- list(item = select, par = itemBank$itemPar[select,], info = itemBank$infoTab[ind,select], criterion = criterion,randomesque=randomesque)
    }
    if (crit == "Urry") {
        indic <- rep(1, nrow(itemBank$itemPar))
        indic[OUT] <- 0
        th <- (1:length(itemBank$theta))[abs(itemBank$theta - 
            theta) == min(abs(itemBank$theta - theta))]
        if (length(th) > 1) 
            th <- min(th)
diff<-abs(itemBank$itemPar[, 2] - theta)
bVal<-sort(diff[indic==1])[min(c(randomesque,sum(indic)))]
keep<-(1:nrow(itemBank$itemPar))[indic == 1 & diff<=bVal]
select<-ifelse(length(keep)==1,keep,sample(keep,1))
res <- list(item = select, par = itemBank$itemPar[select, ], 
            info = itemBank$infoTab[th, select], criterion = criterion,randomesque=randomesque)
    }
    if (crit == "MLWI" | crit == "MPWI") {
        if (length(out) == 1) 
            par <- matrix(itemBank$itemPar[out, ], 1, 4)
        else par <- itemBank$itemPar[out, ]
        ITEMS <- rep(1, nrow(itemBank$itemPar))
        ITEMS[OUT] <- 0
        likInfo <- rep(0, nrow(itemBank$itemPar))
        for (i in 1:nrow(itemBank$itemPar)) {
            if (ITEMS[i] == 1) 
                likInfo[i] <- MWI(itemBank, i, x, par, type = criterion, 
                  lower = parInt[1], upper = parInt[2], nqp = parInt[3], 
                  priorDist = priorDist, priorPar = priorPar)
        }
likVal<-sort(likInfo,decreasing=TRUE)[min(c(randomesque,sum(ITEMS)))]
keep<-(1:length(ITEMS))[likInfo>=likVal]
select<-ifelse(length(keep)==1,keep,sample(keep,1))
res <- list(item = select, par = itemBank$itemPar[select, ], 
info = likInfo[select], criterion = criterion,randomesque=randomesque)
    }
    if (crit == "MEI") {
        items <- rep(1, nrow(itemBank$itemPar))
        items[OUT] <- 0
        infos <- rep(0, length(items))
        for (i in 1:length(items)) {
            if (items[i] > 0) 
                infos[i] <- MEI(itemBank, item = i, x = x, theta = theta, 
                  it = itemBank$itemPar[out, ], method = method, 
                  priorDist = priorDist, priorPar = priorPar, 
                  D = D, range = range, parInt = parInt, infoType = infoType)
        }
infoVal<-sort(infos,decreasing=TRUE)[min(c(randomesque,sum(items)))]
keep<-(1:nrow(itemBank$itemPar))[infos>=infoVal]
select<-ifelse(length(keep)==1,keep,sample(keep,1))
res<-list(item = select, par = itemBank$itemPar[select,], info = infos[select], criterion = criterion,randomesque=randomesque)
    }
    if (crit == "MEPV") {
        items <- rep(1, nrow(itemBank$itemPar))
        items[OUT] <- 0
        epvs <- rep(1000, length(items))
        for (i in 1:length(items)) {
            if (items[i] > 0) 
                epvs[i] <- EPV(itemBank, item = i, x = x, theta = theta, 
                  it = itemBank$itemPar[out, ], priorDist = priorDist, 
                  priorPar = priorPar, D = D, parInt = parInt)
        }
epVal<-sort(epvs)[min(c(randomesque,sum(items)))]
keep<-(1:nrow(itemBank$itemPar))[epvs<=epVal]
select<-ifelse(length(keep)==1,keep,sample(keep,1))
res <- list(item = select, par = itemBank$itemPar[select,], info =epvs[select], criterion = criterion,randomesque=randomesque)
    }
    if (crit == "random") {
        items <- rep(1, nrow(itemBank$itemPar))
        items[OUT] <- 0
        gen <- as.integer(runif(1, 0, 1) * (sum(items))) + 1
        ind <- (1:nrow(itemBank$itemPar))[items > 0][gen]
        res <- list(item = ind, par = itemBank$itemPar[ind, ], 
            info = NA, criterion = criterion, randomesque=randomesque)
    }
if (is.null(cbControl)) res[[6]]<-res[[7]]<-res[[8]]<-NA
else{
res[[6]]<-empProp
postProp<-NULL
for (i in 1:nrGroup) postProp[i]<-length(c(res$item,out)[itemBank$cbGroup[c(res$item,out)]==cbControl$names[i]])
res[[7]]<-postProp/sum(postProp)
res[[8]]<-thProp
}
names(res)[6:8]<-c("prior.prop","post.prop","cb.prop")
return(res)
}
