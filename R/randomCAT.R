randomCAT<-function(trueTheta,itemBank,maxItems=50,start=list(fixItems=NULL,seed=NULL,nrItems=1,theta=0,bw=4,range=c(-4,4)),test=list(method="BM",priorDist="norm",priorPar=c(0,1),range=c(-4,4),D=1,eapPar=c(-4,4,20),itemSelect="info"),stop=list(rule="length",thr=20,alpha=0.05),final=list(method="BM",priorDist="norm",priorPar=c(0,1),range=c(-4,4),D=1,eapPar=c(-4,4,20),alpha=0.05)){
if (testList(start,type="start")$test==FALSE) stop(testList(start,type="start")$message,call.=FALSE)
if (testList(test,type="test")$test==FALSE) stop(testList(test,type="test")$message,call.=FALSE)
if (testList(stop,type="stop")$test==FALSE) stop(testList(stop,type="stop")$message,call.=FALSE)
if (testList(final,type="final")$test==FALSE) stop(testList(final,type="final")$message,call.=FALSE)
startList<-list(fixItems=start$fixItems,seed=start$seed,nrItems=NULL,theta=NULL,bw=NULL,range=c(-4,4))
startList$nrItems<-ifelse(is.null(start$nrItems),1,start$nrItems)
startList$theta<-ifelse(is.null(start$theta),0,start$theta)
startList$bw<-ifelse(is.null(start$bw),4,start$bw)
if (is.null(start$range)==FALSE){
startList$range[1]<-start$range[1]
startList$range[2]<-start$range[2]
}
start<-startList
testList<-list(method=NULL,priorDist=NULL,priorPar=c(0,1),range=c(-4,4),D=1,eapPar=c(-4,4,20),itemSelect="info")
testList$method<-ifelse(is.null(test$method),"BM",test$method)
testList$priorDist<-ifelse(is.null(test$priorDist),"norm",test$priorDist)
if(is.null(test$priorPar)==FALSE){
testList$priorPar[1]<-test$priorPar[1]
testList$priorPar[2]<-test$priorPar[2]
}
if(is.null(test$range)==FALSE){
testList$range[1]<-test$range[1]
testList$range[2]<-test$range[2]
}
testList$D<-ifelse(is.null(test$D),1,test$D)
if (is.null(test$eapPar)==FALSE){
testList$eapPar[1]<-test$eapPar[1]
testList$eapPar[2]<-test$eapPar[2]
testList$eapPar[3]<-test$eapPar[3]
}
testList$itemSelect<-ifelse(is.null(test$itemSelect),"info",test$itemSelect)
test<-testList
stopList<-list(rule=NULL,thr=20,alpha=0.05)
stopList$rule<-ifelse(is.null(stop$rule),"length",stop$rule)
stopList$thr<-ifelse(is.null(stop$thr),20,stop$thr)
stopList$alpha<-ifelse(is.null(stop$alpha),0.05,stop$alpha)
stop<-stopList
finalList<-list(method=NULL,priorDist=NULL,priorPar=c(0,1),range=c(-4,4),D=1,eapPar=c(-4,4,20),alpha=0.05)
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
if (is.null(final$eapPar)==FALSE){
finalList$eapPar[1]<-final$eapPar[1]
finalList$eapPar[2]<-final$eapPar[2]
finalList$eapPar[3]<-final$eapPar[3]
}
finalList$alpha<-ifelse(is.null(final$alpha),0.05,final$alpha)
final<-finalList
pr0<-startItems(itemBank$itemPar,fixItems=start$fixItems,seed=start$seed,nrItems=start$nrItems,theta=start$theta,bw=start$bw,range=start$range)
ITEMS<-pr0$items
PAR<-rbind(pr0$par)
PATTERN<-rbinom(length(ITEMS),1,Pi(trueTheta,PAR,D=test$D)$Pi)
TH<-thetaEst(PAR,PATTERN,D=test$D,method=test$method,priorDist=test$priorDist,priorPar=test$priorPar,range=test$range,eapPar=test$eapPar)
SETH<-semTheta(TH,PAR,x=PATTERN,D=test$D,method=test$method,priorDist=test$priorDist,priorPar=test$priorPar,eapPar=test$eapPar)
thProv<-TH
if (stop$rule=="length") maxLength<-min(c(maxItems,stop$thr))
else maxLength<-maxItems
if (stop$rule=="classification" & (TH-qnorm(1-stop$alpha/2)*SETH>=stop$thr | TH+qnorm(1-stop$alpha/2)*SETH<=stop$thr)){
finalEst<-thetaEst(PAR,PATTERN,D=final$D,method=final$method,priorDist=final$priorDist,priorPar=final$priorPar,range=final$range,eapPar=final$eapPar)
seFinal<-semTheta(finalEst,PAR,x=PATTERN,D=final$D,method=final$method,priorDist=final$priorDist,priorPar=final$priorPar,eapPar=final$eapPar)
confIntFinal<-c(finalEst-qnorm(1-final$alpha/2)*seFinal,finalEst+qnorm(1-final$alpha/2)*seFinal)
endWarning<-FALSE
RES<-list(trueTheta=trueTheta,maxItems=maxItems,testItems=ITEMS,itemPar=PAR,pattern=PATTERN,thetaProv=TH,seProv=SETH,thFinal=finalEst,seFinal=seFinal,ciFinal=confIntFinal,startFixItems=start$fixItems,startSeed=start$seed,startNrItems=start$nrItems,startTheta=start$theta,startBw=start$bw,startbOpt=pr0$bOpt,startRange=start$range,provMethod=test$method,provDist=test$priorDist,provPar=test$priorPar,provRange=test$range,provD=test$D,stopRule=stop$rule,stopThr=stop$thr,stopAlpha=stop$alpha,endWarning=endWarning,finalMethod=final$method,finalDist=final$priorDist,finalPar=final$priorPar,finalRange=final$range,finalD=final$D,finalAlpha=final$alpha)
class(RES)<-"cat"
}
else{
repeat{
pr<-nextItem(itemBank,thProv,out=ITEMS,method=test$itemSelect)
ITEMS<-c(ITEMS,pr$item)
PAR<-rbind(PAR,pr$par)
PATTERN<-c(PATTERN,rbinom(1,1,Pi(trueTheta,rbind(pr$par),D=test$D)$Pi))
thProv<-thetaEst(PAR,PATTERN,D=test$D,method=test$method,priorDist=test$priorDist,priorPar=test$priorPar,range=test$range,eapPar=test$eapPar)
TH<-c(TH,thProv)
seProv<-semTheta(thProv,PAR,x=PATTERN,D=test$D,method=test$method,priorDist=test$priorDist,priorPar=test$priorPar,eapPar=test$eapPar)
SETH<-c(SETH,seProv)
if ((length(ITEMS)>=maxLength) | (stop$rule=="precision" & seProv<=stop$thr) | (stop$rule=="classification" & (thProv-qnorm(1-stop$alpha/2)*seProv>=stop$thr | thProv+qnorm(1-stop$alpha/2)*seProv<=stop$thr))) break
}
finalEst<-thetaEst(PAR,PATTERN,D=final$D,method=final$method,priorDist=final$priorDist,priorPar=final$priorPar,range=final$range,eapPar=final$eapPar)
seFinal<-semTheta(finalEst,PAR,x=PATTERN,D=final$D,method=final$method,priorDist=final$priorDist,priorPar=final$priorPar,eapPar=final$eapPar)
confIntFinal<-c(finalEst-qnorm(1-final$alpha/2)*seFinal,finalEst+qnorm(1-final$alpha/2)*seFinal)
if ((stop$rule=="length" & length(ITEMS)<stop$thr) | (stop$rule=="precision" & seProv>stop$thr) | (stop$rule=="classification" & thProv-qnorm(1-stop$alpha/2)*seProv<stop$thr & thProv+qnorm(1-stop$alpha/2)*seProv>stop$thr)) endWarning<-TRUE
else endWarning<-FALSE
RES<-list(trueTheta=trueTheta,maxItems=maxItems,testItems=ITEMS,itemPar=PAR,pattern=PATTERN,thetaProv=TH,seProv=SETH,thFinal=finalEst,seFinal=seFinal,ciFinal=confIntFinal,startFixItems=start$fixItems,startSeed=start$seed,startNrItems=start$nrItems,startTheta=start$theta,startBw=start$bw,startbOpt=pr0$bOpt,startRange=start$range,provMethod=test$method,provDist=test$priorDist,provPar=test$priorPar,provRange=test$range,provD=test$D,itemSelect=test$itemSelect,stopRule=stop$rule,stopThr=stop$thr,stopAlpha=stop$alpha,endWarning=endWarning,finalMethod=final$method,finalDist=final$priorDist,finalPar=final$priorPar,finalRange=final$range,finalD=final$D,finalAlpha=final$alpha)
class(RES)<-"cat"
}
return(RES)
}


print.cat<-function(x, ...){
cat("Random generation of a CAT response pattern","\n","\n")
cat(" True ability level:",round(x$trueTheta,2),"\n","\n")
cat(" Starting parameters:","\n")
if (is.null(x$startFixItems)==TRUE) nr1<-x$startNrItems
else nr1<-length(x$startFixItems)
cat("   Number of early items:",nr1,"\n")
if (is.null(x$startFixItems)==FALSE) met1<-"Chosen by administrator"
else {
if (is.null(x$startSeed)==FALSE) met1<-"Random selection in item bank"
else {
if (x$startNrItems==1) met1<-"with difficulty level as close as possible" 
else{
met1<-"with difficulty levels as close as possible"
}
}
}
if (nr1==1) cat("   Early item selection:",met1,"\n")
else cat("   Early items selection:",met1,"\n")
if (is.null(x$startFixItems)==FALSE){
if (length(x$startFixItems)==1) met1bis<-paste("    (item administered: ",x$startFixItems,")",sep="")
else{
met1bis<-paste("    (item administered: ",x$startFixItems[1],sep="")
for (i in 2:length(x$startFixItems)) met1bis<-paste(met1bis," - ",x$startFixItems[i],sep="")
met1bis<-paste(met1bis,")",sep="")
}
cat(met1bis,"\n")
}
if (is.null(x$startFixItems)==TRUE & is.null(x$startSeed)==TRUE){
if (x$startNrItems==1) met1bis<-paste("    to value ",x$startTheta,sep="")
else{
met1bis<-paste("    to values ",x$startbOpt[1],sep="")
if (length(x$startbOpt)==2) met1bis<-paste(met1bis," and ",x$startbOpt[2],sep="")
else{
for (i in 2:(length(x$startbOpt)-1)) met1bis<-paste(met1bis,", ",x$startbOpt[i],sep="")
met1bis<-paste(met1bis," and ",x$startbOpt[length(x$startbOpt)],sep="")
}
}
cat(met1bis,"\n")
}
cat("\n","Adaptive test parameters:","\n")
itemSel<-switch(x$itemSelect,info="maximum information criterion",Owen="Owen's approximate Bayes procedure")
cat("   Next item selection method:",itemSel,"\n")
met2<-switch(x$provMethod,
BM="Bayes modal (MAP) estimator",
WL="Weighted likelihood estimator",
ML="Maximum likelihood estimator",
EAP="Expected a posteriori (EAP) estimator")
if (x$provMethod=="BM" | x$provMethod=="EAP"){
met3<-switch(x$provDist,
norm=paste("N(",round(x$provPar[1],2),",",round(x$provPar[2]^2,2),") prior",sep=""),
unif=paste("U(",round(x$provPar[1],2),",",round(x$provPar[2],2),") prior",sep=""),
Jeffreys="Jeffreys' prior")
}
if (x$provMethod=="ML") ra1<-paste("[",round(x$provRange[1],2),",",round(x$provRange[2],2),"]",sep="")
cat("   Provisional ability estimator:",met2,"\n")
if (x$provMethod=="BM" | x$provMethod=="EAP") cat("   Provisional prior distribution:",met3,"\n")
if (x$provMethod=="ML") cat("   Provisional range of ability values:",ra1,"\n")
cat("\n")
cat(" Stopping rule:","\n")
met4<-switch(x$stopRule,
length="length of test",
precision="precision of ability estimate",
classification=paste("classification based on ",100*(1-x$stopAlpha),"% confidence interval",sep=""))
cat("   Stopping criterion:",met4,"\n")
switch(x$stopRule,
length=cat("   Minimum test length:",x$stopThr,"items","\n"),
precision=cat("   Maximum SE value:",round(x$stopThr,2),"\n"),
classification=cat("   Classification threshold:",round(x$stopThr,2),"\n"))
cat("   Maximum test length:",x$maxItems,"items","\n")
mat<-rbind(as.character(1:length(x$testItems)),as.character(x$testItems),round(x$pattern,0))
nra<-length(x$pattern)-length(x$thetaProv)
mat<-rbind(mat,c(rep(NA,nra),round(x$thetaProv,3)),c(rep(NA,nra),round(x$seProv,3)))
rownames(mat)<-c("Nr","Item","Resp.","Est.","SE")
colnames(mat)<-rep("",ncol(mat))
cat("\n","Adaptive test details:","\n")
print(format(mat,justify="right"),quote=FALSE)
cat("\n")
if (x$endWarning==TRUE) cat("WARNING: stopping rule was not satisfied after",x$maxItems,"items!","\n","\n")
cat(" Final results:","\n")
met<-switch(x$finalMethod,
BM="Bayes modal (MAP) estimator",
WL="Weighted likelihood estimator",
ML="Maximum likelihood estimator",
EAP="Expected a posteriori (EAP) estimator")
if (x$finalMethod=="BM" | x$finalMethod=="EAP"){
met2<-switch(x$finalDist,
norm=paste("N(",round(x$finalPar[1],2),",",round(x$finalPar[2]^2,2),") prior",sep=""),
unif=paste("U(",round(x$finalPar[1],2),",",round(x$finalPar[2],2),") prior",sep=""),
Jeffreys="Jeffreys' prior")
}
if (x$finalMethod=="ML") ra1<-paste("[",round(x$finalRange[1],2),",",round(x$finalRange[2],2),"]",sep="")
cat("   Length of adaptive test:",length(x$testItems),"items","\n")
cat("   Final ability estimator:",met,"\n")
if (x$finalMethod=="BM" | x$finalMethod=="EAP") cat("   Final prior distribution:",met2,"\n")
if (x$finalMethod=="ML") cat("   Final range of ability values:",ra1,"\n")
cat("   Final ability estimate (SE):",round(x$thFinal,3),paste("(",round(x$seFinal,3),")",sep=""),"\n")
cat(paste("   ",(1-x$finalAlpha)*100,"% confidence interval: [",round(x$ciFinal[1],3),",",round(x$ciFinal[2],3),"]",sep=""),"\n")
if (x$stopRule=="classification"){
if (x$ciFinal[1]>x$stopThr) mess<-paste("ability is larger than ",round(x$stopThr,2),sep="")
else{
if (x$ciFinal[2]<x$stopThr) mess<-paste("ability is smaller than ",round(x$stopThr,2),sep="")
else mess<-paste("ability is not different from ",round(x$stopThr,2),sep="")
}
cat("   Final subject classification:",mess,"\n","\n")
}
}


plot.cat<-function(x,ci=TRUE,alpha=0.05,thr=NULL, ...){
if (is.null(thr)==FALSE & is.numeric(thr)==FALSE) stop("'thr' must be either numeric or NULL",call.=FALSE)
X<-1:length(x$testItems)
nra<-length(x$pattern)-length(x$thetaProv)
Y<-c(rep(NA,nra),x$thetaProv)
r1<-x$thetaProv-qnorm(1-alpha/2)*x$seProv
r2<-x$thetaProv+qnorm(1-alpha/2)*x$seProv	
if (ci==TRUE) ra<-range(c(r1,r2,thr))
else ra<-range(c(x$thetaProv,thr))
r1<-c(rep(NA,nra),r1)
r2<-c(rep(NA,nra),r2)
plot(X,Y,type="o",xlab="Item",ylab="Ability estimate",ylim=ra,cex=0.7)
if (ci==TRUE){
for (i in 1:length(X)){
lines(rep(i,2),c(r1[i],r2[i]))
lines(c(i-0.2,i+0.2),rep(r1[i],2))
lines(c(i-0.2,i+0.2),rep(r2[i],2))
}
if (is.null(thr)==FALSE) abline(h=thr,lty=2)
}
}

