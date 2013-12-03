randomCAT<-function(trueTheta,itemBank,maxItems=50,cbControl=NULL,start=list(fixItems=NULL,seed=NULL,nrItems=1,theta=0,halfRange=2,startSelect="bOpt"),test=list(method="BM",priorDist="norm",priorPar=c(0,1),range=c(-4,4),D=1,parInt=c(-4,4,33),itemSelect="MFI",infoType="observed",randomesque=1),stop=list(rule="length",thr=20,alpha=0.05),final=list(method="BM",priorDist="norm",priorPar=c(0,1),range=c(-4,4),D=1,parInt=c(-4,4,33),alpha=0.05),save.output=FALSE,output=c("out","default")){
if (!testList(start,type="start")$test) stop(testList(start,type="start")$message,call.=FALSE)
if (!testList(test,type="test")$test) stop(testList(test,type="test")$message,call.=FALSE)
if (!testList(stop,type="stop")$test) stop(testList(stop,type="stop")$message,call.=FALSE)
if (!testList(final,type="final")$test) stop(testList(final,type="final")$message,call.=FALSE)
internalCAT<-function(){
startList<-list(fixItems=start$fixItems,seed=start$seed,nrItems=NULL,theta=NULL,halfRange=2,startSelect="bOpt")
startList$nrItems<-ifelse(is.null(start$nrItems),1,start$nrItems)
startList$theta<-ifelse(is.null(start$theta),0,start$theta)
startList$halfRange<-ifelse(is.null(start$halfRange),2,start$halfRange)
startList$startSelect<-ifelse(is.null(start$startSelect),"bOpt",start$startSelect)
start<-startList
testList<-list(method=NULL,priorDist=NULL,priorPar=c(0,1),range=c(-4,4),D=1,parInt=c(-4,4,33),itemSelect="MFI",infoType="observed",randomesque=1)
testList$method<-ifelse(is.null(test$method),"BM",test$method)
testList$priorDist<-ifelse(is.null(test$priorDist),"norm",test$priorDist)
if(!is.null(test$priorPar)){
testList$priorPar[1]<-test$priorPar[1]
testList$priorPar[2]<-test$priorPar[2]
}
if(!is.null(test$range)){
testList$range[1]<-test$range[1]
testList$range[2]<-test$range[2]
}
testList$D<-ifelse(is.null(test$D),1,test$D)
if (!is.null(test$parInt)){
testList$parInt[1]<-test$parInt[1]
testList$parInt[2]<-test$parInt[2]
testList$parInt[3]<-test$parInt[3]
}
testList$itemSelect<-ifelse(is.null(test$itemSelect),"MFI",test$itemSelect)
testList$infoType<-ifelse(is.null(test$infoType),"observed",test$infoType)
testList$randomesque<-ifelse(is.null(test$randomesque),1,test$randomesque)
test<-testList
stopList<-list(rule=NULL,thr=20,alpha=0.05)
stopList$rule<-ifelse(is.null(stop$rule),"length",stop$rule)
stopList$thr<-ifelse(is.null(stop$thr),20,stop$thr)
stopList$alpha<-ifelse(is.null(stop$alpha),0.05,stop$alpha)
stop<-stopList
if (stop$rule=="length" & stop$thr>maxItems) stop("'maxItems' is smaller than test length criterion 'stop$thr'",call.=FALSE)
finalList<-list(method=NULL,priorDist=NULL,priorPar=c(0,1),range=c(-4,4),D=1,parInt=c(-4,4,33),alpha=0.05)
finalList$method<-ifelse(is.null(final$method),"BM",final$method)
finalList$priorDist<-ifelse(is.null(final$priorDist),"norm",final$priorDist)
if(is.null(final$priorPar)==FALSE){
finalList$priorPar[1]<-final$priorPar[1]
finalList$priorPar[2]<-final$priorPar[2]
}
if(is.null(final$range)==FALSE){
finalList$range[1]<-final$range[1]
finalList$range[2]<-final$range[2]
}
finalList$D<-ifelse(is.null(final$D),1,final$D)
if (is.null(final$parInt)==FALSE){
finalList$parInt[1]<-final$parInt[1]
finalList$parInt[2]<-final$parInt[2]
finalList$parInt[3]<-final$parInt[3]
}
finalList$alpha<-ifelse(is.null(final$alpha),0.05,final$alpha)
final<-finalList
pr0<-startItems(itemBank=itemBank,fixItems=start$fixItems,seed=start$seed,nrItems=start$nrItems,theta=start$theta,halfRange=start$halfRange,startSelect=start$startSelect)
ITEMS<-pr0$items
PAR<-rbind(pr0$par)
PATTERN<-rbinom(length(ITEMS),1,Pi(trueTheta,PAR,D=test$D)$Pi)
TH<-thetaEst(PAR,PATTERN,D=test$D,method=test$method,priorDist=test$priorDist,priorPar=test$priorPar,range=test$range,parInt=test$parInt)
SETH<-semTheta(TH,PAR,x=PATTERN,D=test$D,method=test$method,priorDist=test$priorDist,priorPar=test$priorPar,parInt=test$parInt)
thProv<-TH
if (stop$rule=="length") maxLength<-min(c(maxItems,stop$thr))
else maxLength<-maxItems
if (stop$rule=="classification" & (TH-qnorm(1-stop$alpha/2)*SETH>=stop$thr | TH+qnorm(1-stop$alpha/2)*SETH<=stop$thr)){
finalEst<-thetaEst(PAR,PATTERN,D=final$D,method=final$method,priorDist=final$priorDist,priorPar=final$priorPar,range=final$range,parInt=final$parInt)
seFinal<-semTheta(finalEst,PAR,x=PATTERN,D=final$D,method=final$method,priorDist=final$priorDist,priorPar=final$priorPar,parInt=final$parInt)
confIntFinal<-c(finalEst-qnorm(1-final$alpha/2)*seFinal,finalEst+qnorm(1-final$alpha/2)*seFinal)
endWarning<-FALSE
RES<-list(trueTheta=trueTheta,maxItems=maxItems,testItems=ITEMS,itemPar=PAR,pattern=PATTERN,thetaProv=TH,seProv=SETH,thFinal=finalEst,seFinal=seFinal,ciFinal=confIntFinal,startFixItems=start$fixItems,startSeed=start$seed,startNrItems=start$nrItems,startTheta=start$theta,startHalfRange=start$halfRange,startThStart=pr0$thStart,startSelect=start$startSelect,provMethod=test$method,provDist=test$priorDist,provPar=test$priorPar,provRange=test$range,provD=test$D,itemSelect=test$itemSelect,infoType=test$infoType,randomesque=test$randomesque,cbControl=cbControl,cbGroup=itemBank$cbGroup,stopRule=stop$rule,stopThr=stop$thr,stopAlpha=stop$alpha,endWarning=endWarning,finalMethod=final$method,finalDist=final$priorDist,finalPar=final$priorPar,finalRange=final$range,finalD=final$D,finalAlpha=final$alpha,save.output=save.output,output=output)
class(RES)<-"cat"
}
else{
repeat{
pr<-nextItem(itemBank,thProv,out=ITEMS,x=PATTERN,criterion=test$itemSelect,method=test$method,parInt=test$parInt,priorDist=test$priorDist,priorPar=test$priorPar,infoType=test$infoType,D=test$D,range=test$range,randomesque=test$randomesque,cbControl=cbControl)
ITEMS<-c(ITEMS,pr$item)
PAR<-rbind(PAR,pr$par)
PATTERN<-c(PATTERN,rbinom(1,1,Pi(trueTheta,rbind(pr$par),D=test$D)$Pi))
thProv<-thetaEst(PAR,PATTERN,D=test$D,method=test$method,priorDist=test$priorDist,priorPar=test$priorPar,range=test$range,parInt=test$parInt)
TH<-c(TH,thProv)
seProv<-semTheta(thProv,PAR,x=PATTERN,D=test$D,method=test$method,priorDist=test$priorDist,priorPar=test$priorPar,parInt=test$parInt)
SETH<-c(SETH,seProv)
if ((length(ITEMS)>=maxLength) | (stop$rule=="precision" & seProv<=stop$thr) | (stop$rule=="classification" & (thProv-qnorm(1-stop$alpha/2)*seProv>=stop$thr | thProv+qnorm(1-stop$alpha/2)*seProv<=stop$thr))) break
}
finalEst<-thetaEst(PAR,PATTERN,D=final$D,method=final$method,priorDist=final$priorDist,priorPar=final$priorPar,range=final$range,parInt=final$parInt)
seFinal<-semTheta(finalEst,PAR,x=PATTERN,D=final$D,method=final$method,priorDist=final$priorDist,priorPar=final$priorPar,parInt=final$parInt)
confIntFinal<-c(finalEst-qnorm(1-final$alpha/2)*seFinal,finalEst+qnorm(1-final$alpha/2)*seFinal)
if ((stop$rule=="length" & length(ITEMS)<stop$thr) | (stop$rule=="precision" & seProv>stop$thr) | (stop$rule=="classification" & thProv-qnorm(1-stop$alpha/2)*seProv<stop$thr & thProv+qnorm(1-stop$alpha/2)*seProv>stop$thr)) endWarning<-TRUE
else endWarning<-FALSE
RES<-list(trueTheta=trueTheta,maxItems=maxItems,testItems=ITEMS,itemPar=PAR,pattern=PATTERN,thetaProv=TH,seProv=SETH,thFinal=finalEst,seFinal=seFinal,ciFinal=confIntFinal,startFixItems=start$fixItems,startSeed=start$seed,startNrItems=start$nrItems,startTheta=start$theta,startHalfRange=start$halfRange,startThStart=pr0$thStart,startSelect=start$startSelect,provMethod=test$method,provDist=test$priorDist,provPar=test$priorPar,provRange=test$range,provD=test$D,itemSelect=test$itemSelect,infoType=test$infoType,randomesque=test$randomesque,cbControl=cbControl,cbGroup=itemBank$cbGroup,stopRule=stop$rule,stopThr=stop$thr,stopAlpha=stop$alpha,endWarning=endWarning,finalMethod=final$method,finalDist=final$priorDist,finalPar=final$priorPar,finalRange=final$range,finalD=final$D,finalAlpha=final$alpha,save.output=save.output,output=output)
class(RES)<-"cat"
}
return(RES)
}
   resToReturn <- internalCAT()
    if (save.output) {
        if (output[2] == "default") 
            wd <- paste(getwd(), "/", sep = "")
        else wd <- output[2]
        fileName <- paste(wd, output[1], ".txt", sep = "")
        capture.output(resToReturn, file = fileName)
    }
    return(resToReturn)
}

print.cat<-function (x, ...) 
{
    cat("Random generation of a CAT response pattern", "\n", 
        "\n")
    cat(" True ability level:", round(x$trueTheta, 2), "\n", 
        "\n")
    cat(" Starting parameters:", "\n")
    if (is.null(x$startFixItems)) 
        nr1 <- x$startNrItems
    else nr1 <- length(x$startFixItems)
    cat("   Number of early items:", nr1, "\n")
    if (!is.null(x$startFixItems)) 
        met1 <- "Chosen by administrator"
    else {
        if (!is.null(x$startSeed)) 
            met1 <- "Random selection in item bank"
        else {
            if (x$startSelect == "bOpt") {
                if (x$startNrItems == 1) 
                  met1 <- "matching item difficulty to starting ability"
                else met1 <- "matching item difficulties to starting abilities"
            }
            else {
                if (x$startNrItems == 1) 
                  met1 <- "maximum informative item for starting ability"
                else met1 <- "maximum informative items for starting abilities"
            }
        }
    }
    if (nr1 == 1) 
        cat("   Early item selection:", met1, "\n")
    else cat("   Early items selection:", met1, "\n")
    if (!is.null(x$startFixItems)) {
        if (length(x$startFixItems) == 1) 
            met1bis <- paste("   Item administered: ", x$startFixItems, 
                sep = "")
        else {
            met1bis <- paste("   Items administered: ", x$startFixItems[1], 
                sep = "")
            if (length(x$startFixItems) == 2) 
                met1bis <- paste(met1bis, " and ", x$startFixItems[2], 
                  sep = "")
            else {
                for (i in 2:(length(x$startFixItems) - 1)) met1bis <- paste(met1bis, 
                  ", ", x$startFixItems[i], sep = "")
                met1bis <- paste(met1bis, " and ", x$startFixItems[length(x$startFixItems)], 
                  sep = "")
            }
        }
        cat(met1bis, "\n")
    }
    if (is.null(x$startFixItems) & !is.null(x$startSeed)) {
        if (x$startNrItems == 1) 
            met1bis <- paste("   Item administered: ", x$testItems[1], 
                sep = "")
        else {
            met1bis <- paste("   Items administered: ", x$testItems[1], 
                sep = "")
            if (x$startNrItems == 2) 
                met1bis <- paste(met1bis, " and ", x$testItems[2], 
                  sep = "")
            else {
                for (i in 2:(x$startNrItems - 1)) met1bis <- paste(met1bis, 
                  ", ", x$testItems[i], sep = "")
                met1bis <- paste(met1bis, " and ", x$testItems[x$startNrItems], 
                  sep = "")
            }
        }
        cat(met1bis, "\n")
    }
    if (is.null(x$startFixItems) & is.null(x$startSeed)) {
        if (x$startNrItems == 1) 
            met1bis <- paste("   Starting ability: ", x$startThStart, 
                sep = "")
        else {
            met1bis <- paste("   Starting abilities: ", sort(x$startThStart)[1], 
                sep = "")
            if (length(x$startThStart) == 2) 
                met1bis <- paste(met1bis, " and ", sort(x$startThStart)[2], 
                  sep = "")
            else {
                for (i in 2:(length(x$startThStart) - 1)) met1bis <- paste(met1bis, 
                  ", ", sort(x$startThStart)[i], sep = "")
                met1bis <- paste(met1bis, " and ", sort(x$startThStart)[length(x$startThStart)], 
                  sep = "")
            }
        }
        cat(met1bis, "\n")
        if (is.null(x$startFixItems) & is.null(x$startSeed)){
if (x$startNrItems > 2) {
            met1ter <- paste("   Order of starting abilities administration: ", 
                x$startThStart[1], sep = "")
            for (i in 2:(length(x$startThStart) - 1)) met1ter <- paste(met1ter, 
                ", ", x$startThStart[i], sep = "")
            met1ter <- paste(met1ter, " and ", x$startThStart[length(x$startThStart)], 
                sep = "")
            cat(met1ter, "\n")
}
else{
if (x$startNrItems == 2) {
met1ter <- paste("   Order of starting abilities administration: ", 
                x$startThStart[1]," and ",x$startThStart[2], sep = "")
            cat(met1ter, "\n")
}
}
        }
    }
    cat("\n", "Adaptive test parameters:", "\n")
    itemSel <- switch(x$itemSelect, MFI = "maximum Fisher information", 
        Urry = "Urry's procedure", MLWI = "Maximum likelihood weighted information (MLWI)", 
        MPWI = "Maximum posterior weighted information (MPWI)", 
        MEI = "Maximum expected information (MEI)", MEPV = "Minimum Expected Posterior Variance (MEPV)", 
        random = "Random selection")
    cat("   Next item selection method:", itemSel, "\n")
    if (x$itemSelect == "MEI") {
        infTyp <- switch(x$infoType, observed = "observed information function", 
            Fisher = "Fisher information function")
        cat("   Type of information:", infTyp, "\n")
    }
    met2 <- switch(x$provMethod, BM = "Bayes modal (MAP) estimator", 
        WL = "Weighted likelihood estimator", ML = "Maximum likelihood estimator", 
        EAP = "Expected a posteriori (EAP) estimator")
    if (x$provMethod == "BM" | x$provMethod == "EAP" | x$itemSelect == 
        "MLWI" | x$itemSelect == "MPWI") {
        met3 <- switch(x$provDist, norm = paste("N(", round(x$provPar[1], 
            2), ",", round(x$provPar[2]^2, 2), ") prior", sep = ""), 
            unif = paste("U(", round(x$provPar[1], 2), ",", round(x$provPar[2], 
                2), ") prior", sep = ""), Jeffreys = "Jeffreys' prior")
    }
    if (x$provMethod == "ML") 
        ra1 <- paste("[", round(x$provRange[1], 2), ",", round(x$provRange[2], 
            2), "]", sep = "")
    cat("   Provisional ability estimator:", met2, "\n")
    if (x$provMethod == "BM" | x$provMethod == "EAP" | x$itemSelect == 
        "MPWI") 
        cat("   Provisional prior ability distribution:", met3, 
            "\n")
    if (x$provMethod == "ML") 
        cat("   Provisional range of ability values:", ra1, "\n")
    cat("\n")
    cat(" Stopping rule:", "\n")
    met4 <- switch(x$stopRule, length = "length of test", precision = "precision of ability estimate", 
        classification = paste("classification based on ", 100 * 
            (1 - x$stopAlpha), "% confidence interval", sep = ""))
    cat("   Stopping criterion:", met4, "\n")
    switch(x$stopRule, precision = cat("   Maximum SE value:", 
        round(x$stopThr, 2), "\n"), classification = cat("   Classification threshold:", 
        round(x$stopThr, 2), "\n"))
    cat("   Maximum test length:", ifelse(x$stopRule == "length", 
        min(c(x$maxItems, x$stopThr)), x$maxItems), "items", 
        "\n")
    cat("\n", "Item exposure control:", "\n")
    cat("   Method: 'randomesque'","\n")
    cat("   Number of 'randomesque' items: ",x$randomesque,"\n",sep="")

    cat("\n", "Content balancing control:", "\n")
    if (is.null(x$cbControl)) cat("   No control for content balancing","\n")
    else{
    cat("   Expected proportions of items per subgroup:","\n","\n")
    mat<-rbind(round(x$cbControl$props/sum(x$cbControl$props),3))
    rownames(mat)<-""
    colnames(mat)<-x$cbControl$names
    print(format(mat, justify = "right"), quote = FALSE)
    cat("\n")
}
    cat("\n", "Adaptive test details:", "\n")
    mat <- rbind(as.character(1:length(x$testItems)), as.character(x$testItems), 
        round(x$pattern, 0))
    nra <- length(x$pattern) - length(x$thetaProv)
    mat <- rbind(mat, c(rep(NA, nra), round(x$thetaProv, 3)), 
        c(rep(NA, nra), round(x$seProv, 3)))
    rownames(mat) <- c("Nr", "Item", "Resp.", "Est.", "SE")
    colnames(mat) <- rep("", ncol(mat))
    print(format(mat, justify = "right"), quote = FALSE)
    cat("\n")
    if (x$endWarning == TRUE) 
        cat("WARNING: stopping rule was not satisfied after", 
            x$maxItems, "items!", "\n", "\n")
    cat(" Final results:", "\n")
    met <- switch(x$finalMethod, BM = "Bayes modal (MAP) estimator", 
        WL = "Weighted likelihood estimator", ML = "Maximum likelihood estimator", 
        EAP = "Expected a posteriori (EAP) estimator")
    if (x$finalMethod == "BM" | x$finalMethod == "EAP") {
        met2 <- switch(x$finalDist, norm = paste("N(", round(x$finalPar[1], 
            2), ",", round(x$finalPar[2]^2, 2), ") prior", sep = ""), 
            unif = paste("U(", round(x$finalPar[1], 2), ",", 
                round(x$finalPar[2], 2), ") prior", sep = ""), 
            Jeffreys = "Jeffreys' prior")
    }
    if (x$finalMethod == "ML") 
        ra1 <- paste("[", round(x$finalRange[1], 2), ",", round(x$finalRange[2], 
            2), "]", sep = "")
    cat("   Length of adaptive test:", length(x$testItems), "items", 
        "\n")
    cat("   Final ability estimator:", met, "\n")
    if (x$finalMethod == "BM" | x$finalMethod == "EAP") 
        cat("   Final prior distribution:", met2, "\n")
    if (x$finalMethod == "ML") 
        cat("   Final range of ability values:", ra1, "\n")
    cat("   Final ability estimate (SE):", round(x$thFinal, 3), 
        paste("(", round(x$seFinal, 3), ")", sep = ""), "\n")
    cat(paste("   ", (1 - x$finalAlpha) * 100, "% confidence interval: [", 
        round(x$ciFinal[1], 3), ",", round(x$ciFinal[2], 3), 
        "]", sep = ""), "\n")
    if (x$stopRule == "classification") {
        if (x$ciFinal[1] > x$stopThr) 
            mess <- paste("ability is larger than ", round(x$stopThr, 
                2), sep = "")
        else {
            if (x$ciFinal[2] < x$stopThr) 
                mess <- paste("ability is smaller than ", round(x$stopThr, 
                  2), sep = "")
            else mess <- paste("ability is not different from ", 
                round(x$stopThr, 2), sep = "")
        }
        cat("   Final subject classification:", mess, "\n")
    }
    if (!is.null(x$cbControl)){
    cat("   Proportions of items per subgroup (expected and observed):","\n","\n")
    mat<-rbind(round(x$cbControl$props/sum(x$cbControl$props),3))
    nr<-NULL
    for (i in 1:length(x$cbControl$names)) nr[i]<-length(x$testItems[x$cbGroup[x$testItems]==x$cbControl$names[i]])
    nr<-nr/sum(nr)
    mat<-rbind(mat,round(nr,3))
    rownames(mat)<-c("Exp.","Obs.")
    colnames(mat)<-x$cbControl$names
    print(format(mat, justify = "right"), quote = FALSE)
    cat("\n")
    cat("   Items administered per subgroup:","\n","\n")
    for (i in 1:length(x$cbControl$names)) {
    if (length(x$testItems[x$cbGroup[x$testItems]==x$cbControl$names[i]])==0) mess<-"none"
    else{
    its<-sort(x$testItems[x$cbGroup[x$testItems]==x$cbControl$names[i]])
    mess<-its[1]
    if (length(its)>1){
    for (j in 2:length(its)) mess<-paste(mess,", ",its[j],sep="")
}
}
    cat("   ",x$cbControl$names[i],": ",mess,"\n",sep="")
}
}
  if (!x$save.output) 
        cat("\n","Output was not captured!", "\n")
    else {
        if (x$output[2] == "default") 
            wd <- paste(getwd(), "/", sep = "")
        else wd <- x$output[2]
        fileName <- paste(wd, x$output[1], ".txt", sep = "")
        cat("\n","Output was captured and saved into file", "\n", 
            " '", fileName, "'", "\n", "on ",as.character(Sys.Date()),"\n","\n", sep = "")
    }

}



plot.cat<-function (x, ci = FALSE, alpha = 0.05, trueTh = TRUE, classThr = NULL, save.plot = FALSE, 
    save.options = c("plot", "default", "pdf"), ...) 
{
 internalCAT<-function(){
   res <- x
    if (is.logical(trueTh) == FALSE) 
        stop("'trueTh' must be either TRUE or FALSE", call. = FALSE)
    if (is.logical(ci) == FALSE) 
        stop("'ci' must be either TRUE or FALSE", call. = FALSE)
    if (is.numeric(classThr) == FALSE & is.null(classThr) == 
        FALSE) 
        stop("'classThr' must be either a  numeric threshold or NULL", 
            call. = FALSE)
    X <- 1:length(res$testItems)
    nra <- length(res$pattern) - length(res$thetaProv)
    Y <- c(rep(NA, nra), res$thetaProv)
    r1 <- res$thetaProv - qnorm(1 - alpha/2) * res$seProv
    r2 <- res$thetaProv + qnorm(1 - alpha/2) * res$seProv
    vectRange <- c(res$thetaProv, res$trueTheta)
    if (ci == TRUE) 
        vectRange <- c(vectRange, r1, r2)
    if (is.null(classThr) == FALSE) 
        vectRange <- c(vectRange, classThr)
    ra <- range(vectRange)
    ra[1] <- ra[1] - 0.2
    ra[2] <- ra[2] + 0.2
    r1 <- c(rep(NA, nra), r1)
    r2 <- c(rep(NA, nra), r2)
    plot(X, Y, type = "o", xlab = "Item", ylab = "Ability estimate", 
        ylim = ra, cex = 0.7)
    if (ci == TRUE) {
        for (i in 1:length(X)) {
            lines(rep(i, 2), c(r1[i], r2[i]), lty = 3)
            lines(c(i - 0.2, i + 0.2), rep(r1[i], 2))
            lines(c(i - 0.2, i + 0.2), rep(r2[i], 2))
        }
    }
    if (trueTh == TRUE) 
        abline(h = res$trueTheta)
    if (is.null(classThr) == FALSE) 
        abline(h = classThr, lty = 2)
}
    internalCAT()
    if (save.plot) {
        plotype <- NULL
        if (save.options[3] == "pdf") 
            plotype <- 1
        if (save.options[3] == "jpeg") 
            plotype <- 2
        if (is.null(plotype)) 
            cat("Invalid plot type (should be either 'pdf' or 'jpeg').", 
                "\n", "The plot was not captured!", "\n")
        else {
            if (save.options[2] == "default") 
                wd <- paste(getwd(), "/", sep = "")
            else wd <- save.options[2]
            nameFile <- paste(wd, save.options[1], switch(plotype, 
                `1` = ".pdf", `2` = ".jpg"), sep = "")
            if (plotype == 1) {
                {
                  pdf(file = nameFile)
                  internalCAT()
                }
                dev.off()
            }
            if (plotype == 2) {
                {
                  jpeg(filename = nameFile)
                  internalCAT()
                }
                dev.off()
            }
            cat("The plot was captured and saved into", "\n", 
                " '", nameFile, "'", "\n", "\n", sep = "")
        }
    }
    else cat("The plot was not captured!", "\n", sep = "")
}