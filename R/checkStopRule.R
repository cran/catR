checkStopRule<-function(th,se,N,it=NULL,model=NULL,D=1,stop){
res<-FALSE
res.rule<-NULL
for (i in 1:length(stop$rule)){
ind<-switch(stop$rule[i],length=1,precision=2,classification=3,minInfo=4)
if (ind==1){
if (N>=stop$thr[i]){
res<-TRUE
res.rule<-c(res.rule,stop$rule[i])
}}
if (ind==2){
if (se<=stop$thr[i]){
res<-TRUE
res.rule<-c(res.rule,stop$rule[i])
}}
if (ind==3){
if (th-qnorm(1-stop$alpha/2)*se>=stop$thr[i] | th+qnorm(1-stop$alpha/2)*se<=stop$thr[i]){
res<-TRUE
res.rule<-c(res.rule,stop$rule[i])
}}
if (ind==4){
info<-Ii(th,it,model=model,D=D)$Ii
if (max(info)<=stop$thr[i]){
res<-TRUE
res.rule<-c(res.rule,stop$rule[i])
}}
}
return(list(decision=res,rule=res.rule))
}
