createItemBank<-function(items=100,model="4PL",thMin=-4,thMax=4,step=0.01, seed=1, D=1){
if (is.matrix(items)==TRUE) itemPar<-items
else{
set.seed(seed)
b<-rnorm(items)
if (model!="1PL") a<-rnorm(items,1,0.2)
else a<-rep(1,items)
if (model=="3PL" | model=="4PL") c<-runif(items,0,0.25)
else c<-rep(0,items)
if (model=="4PL") d<-runif(items,0.75,1)
else d<-rep(1,items)
itemPar<-cbind(a,b,c,d)
}
colnames(itemPar)<-c("a","b","c","d")
theta<-seq(from=thMin,to=thMax,by=step)
infoTab<-matrix(NA,length(theta),nrow(itemPar))
for (i in 1:length(theta)) infoTab[i,]<-Ii(theta[i],itemPar,D=D)$Ii
res<-list(itemPar=itemPar,theta=theta,infoTab=infoTab)
class(res)<-"itBank"
return(res)}
