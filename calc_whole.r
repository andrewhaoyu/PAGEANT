source('solverf4.r')
source('PowerCalc_Rare.r')
source('PowerCalc_Rare_whole.r')
source('trgamma.r')

calcGenomeLevel = function(K,m,grid,epr,alphaT,Total,CASE=NULL,CONTROL=NULL,PC=NA,TEST = test,QT=NULL,nameEsseble=NULL,JJ=NA){
Esp = epr*K


alphav=seq(10^(-10),1,length=1000000)
PP = 1-pgamma(0.00001/epr,alphav,alphav)
topval = max(PP)
alM = alphav[which(PP==topval)]
twoSolutions = PP[1000000]<topval
Prp =seq(1/2*epr,topval,length.out=grid)
Obj = transform(Total)

EVTV = NULL

PM1 = NULL
PL1 = NULL

PM2 = NULL
PL2 = NULL

PM3 = NULL
PL3 = NULL

JV = NULL
MAFV = NULL
EffectV = NULL


if (grid==20){perm=2500}
if (grid==50){perm=5000}
if (grid==100){perm=10000}

for (e in 1:grid){
cat('Running ',e,'th EV distribution \n')
pr = Prp[e]

if (!(twoSolutions & pr>PP[1000000])){
#cat('One solution','\n')
c1=proc.time()
aa = which(PP<pr)
ss = PP[aa[length(aa)]]
ss = max(ss,0.00000001) 
aa = which(PP>=pr)
ee=  PP[aa[1]]
sol = get_solution(epr,pr,ss,ee)
c2=proc.time()
cat('Sol:',c2-c1,'\n')

alpha = sol$est[1]
beta = sol$est[2]
#cat(alpha,beta,alpha/beta*K,'\n')
if (!is.na(alpha) & !is.na(beta) & ((alpha/beta*K)<=0.5)){
i=1
kEV=NULL
b1 = proc.time()
while (i<=2500){
EVV = rtrunc(K,'gamma',0,0.01,alpha,beta) 
if(sum(EVV)<=0.5){
kEV = c(kEV,sum(EVV))
i=i+1
}
}
#cat(mean(kEV),'\n')
if (abs(mean(kEV)-Esp)/(sd(kEV)/sqrt(5000))<3){

b2=proc.time()
cat('Est',b2-b1,'\n')
EVTV= c(EVTV,mean(kEV))
a1=proc.time()
result <- get_AproxQ2(alpha,beta,K,alphaT,Total,Obj,perm,CASE=CASE,CONTROL=CONTROL,PC=PC,TEST = TEST,QT=QT,nameEsseble=NULL,JJ=JJ) 

JV = c(result$JV)
MAFV = c(result$MAFV)
EffectV = c(EffectV,mean(result$Effect>0.001)*K)


result$a[result$a==0]=10^(-16)
result$a[result$a==1]=1-10^(-16)
PM1 = c(PM1,mean(result$a))
PL1 = c(PL1,mean(log(1-result$a)))

result$b[result$b==0]=10^(-16)
result$b[result$b==1]=1-10^(-16)

PM2 = c(PM2,mean(result$b))
PL2 = c(PL2,mean(log(1-result$b)))
result$c[result$c==0]=10^(-16)
result$c[result$c==1]=1-10^(-16)
PM3 = c(PM3,mean(result$c))
PL3 =c(PL3,mean(log(1-result$c)))

}
}
}else{
aa = which(PP<pr & (1:length(PP))<which(PP==topval))
s1 = alphav[length(aa)]
e1 = alphav[length(aa)+1]
aa = which(PP<pr & (1:length(PP))>=which(PP==topval))
e2 = alphav[aa[1]]
s2 = alphav[aa[1]-1]
#cat('Two solutions','\n')
sol1 = get_solution(epr,pr,(s1),(e1))
sol2= get_solution(epr,pr,(s2),(e2))

alpha = sol1$est[1]
beta = sol1$est[2]
#cat(alpha,beta,alpha/beta*K,'\n')
if (!is.na(alpha) & !is.na(beta) & ((alpha/beta*K)<=0.5)){
i=1
b1 = proc.time()
kEV=NULL
while (i<=5000){
EVV = rtrunc(K,'gamma',0,0.01,alpha,beta) 
if(sum(EVV)<=0.5){
kEV = c(kEV,sum(EVV))
i=i+1
}
}

if (abs(mean(kEV)-Esp)/(sd(kEV)/sqrt(5000))<3){

b2=proc.time()
cat('Est',b2-b1,'\n')
EVTV= c(EVTV,mean(kEV))
a1=proc.time()
result <- get_AproxQ2(alpha,beta,K,alphaT,Total,Obj,perm,CASE=CASE,CONTROL=CONTROL,PC=PC,TEST = TEST,QT=QT,nameEsseble=NULL,JJ=JJ) 
 
JV = c(result$JV)
MAFV = c(result$MAFV)
EffectV = c(EffectV,mean(result$Effect>0.001)*K)
 
result$a[result$a==0]=10^(-16)
result$a[result$a==1]=1-10^(-16)
PM1 = c(PM1,mean(result$a))
PL1 = c(PL1,mean(log(1-result$a)))

result$b[result$b==0]=10^(-16)
result$b[result$b==1]=1-10^(-16)

PM2 = c(PM2,mean(result$b))
PL2 = c(PL2,mean(log(1-result$b)))
result$c[result$c==0]=10^(-16)
result$c[result$c==1]=1-10^(-16)
PM3 = c(PM3,mean(result$c))
PL3 =c(PL3,mean(log(1-result$c)))

}
}

alpha = sol2$est[1]
beta = sol2$est[2]
if (!is.na(alpha) & !is.na(beta) & ((alpha/beta*K)<=0.5)){
i=1
b1 = proc.time()
kEV=NULL
while (i<=5000){
EVV = rtrunc(K,'gamma',0,0.01,alpha,beta) 
if(sum(EVV)<=0.5){
kEV = c(kEV,sum(EVV))
i=i+1
}
}

if (abs(mean(kEV)-Esp)/(sd(kEV)/sqrt(5000))<3){

b2=proc.time()
cat('Est',b2-b1,'\n')
EVTV= c(EVTV,mean(kEV))
a1=proc.time()
result <- get_AproxQ2(alpha,beta,K,alphaT,Total,Obj,perm,CASE=CASE,CONTROL=CONTROL,PC=PC,TEST = TEST,QT=QT,nameEsseble=NULL,JJ=JJ) 

JV = c(result$JV)
MAFV = c(result$MAFV)
EffectV = c(EffectV,mean(result$Effect>0.001)*K)

result$a[result$a==0]=10^(-16)
result$a[result$a==1]=1-10^(-16)
PM1 = c(PM1,mean(result$a))
PL1 = c(PL1,mean(log(1-result$a)))

result$b[result$b==0]=10^(-16)
result$b[result$b==1]=1-10^(-16)

PM2 = c(PM2,mean(result$b))
PL2 = c(PL2,mean(log(1-result$b)))
result$c[result$c==0]=10^(-16)
result$c[result$c==1]=1-10^(-16)
PM3 = c(PM3,mean(result$c))
PL3 =c(PL3,mean(log(1-result$c)))

}

}

}
}
mm=m
Prb1=NULL
Prb2=NULL
Prb3=NULL
L = length(PM1)
for (i in 1:L){
prb1=NULL
prb2=NULL
prb3=NULL
for (m in 0:K){ 
prb1 = c(prb1,choose(K,m)*PM1[i]^m*exp((K-m)*PL1[i]))
prb2 = c(prb2,choose(K,m)*PM2[i]^m*exp((K-m)*PL2[i]))
prb3 = c(prb3,choose(K,m)*PM3[i]^m*exp((K-m)*PL3[i]))
}
Prb1 = rbind(Prb1,prb1)
Prb2 = rbind(Prb2,prb2)
Prb3 = rbind(Prb3,prb3)
}
E1 = PM1*K
E2 = PM2*K
E3 = PM3*K

modelS1_E = list(min=min(E1),max=max(E1))
modelS2_E = list(min=min(E2),max=max(E2))
modelS3_E = list(min=min(E3),max=max(E3))
m=mm
if(m==0){
  Prb1M = Prb1[,1]
  Prb2M = Prb2[,1]
  Prb3M = Prb3[,1]
}else{
  Prb1M = rowSums(Prb1[,1:(m+1)])
  Prb2M = rowSums(Prb2[,1:(m+1)])
  Prb3M = rowSums(Prb3[,1:(m+1)])
}

#
# USE FOLOWING VARIABLES FOR HISTOGRAM !!!
# USE FOLOWING VARIABLES FOR HISTOGRAM !!!
# USE FOLOWING VARIABLES FOR HISTOGRAM !!!
# MAFV-  Empirical distibution of MAF
# JV = Empirical distibution of Size of a Locus
# EffectV - Empirical distibution for # Loci with EV>0.1%


Obj = list(modelS1_E=modelS1_E,modelS2_E=modelS2_E,modelS3_E=modelS3_E,prb1M =list(min=min(Prb1M),max=max(Prb1M)),prb2M =list(min=min(Prb2M),max=max(Prb2M)),prb3M =list(min=min(Prb3M),max=max(Prb3M))) 

places <- 2
S1 <- round(as.numeric(c(Obj[[1]][1],Obj[[1]][2],Obj[[4]][1],Obj[[4]][2])),places)
S2 <- round(as.numeric(c(Obj[[2]][1],Obj[[2]][2],Obj[[5]][1],Obj[[5]][2])),places)
S3 <- round(as.numeric(c(Obj[[3]][1],Obj[[3]][2],Obj[[6]][1],Obj[[6]][2])),places)
Item <- c("Min expected number of discoveries ","Max expected number of discoveries",paste0("min probability of ",m," or less discoveries"),paste0("max probability of ",m," or less discoveries"))
result <- data.frame(Item=Item,S1=S1,S2=S2,S3=S3)
 colnames(result) <- c("","Scenario S1","Scenario S2","Scenario S3")
 
 data1 <- data.frame(MAFV=MAFV)
 p1 <- ggplot(data1,aes(MAFV))+
   geom_histogram(aes(x=data1$MAFV,y=(..count..)/sum(..count..)),
                  fill="dodgerblue4",
                  alpha=0.75
   )+
   theme_bw()+
   theme_new()+
   labs(title="The Histogram of Minor Allele Frequency (MAF)",x="Empirical distibution of Minor Allele Frequency (MAF)")
 data2 <- data.frame(JV =JV)
 p2 <- ggplot(data2,aes(JV))+
   geom_histogram(aes(x=data2$JV,y=(..count..)/sum(..count..)),
                  fill="#c0392b",
                  alpha=0.75
   )+
   theme_bw()+
   theme_new()+
   labs(title="The Histogram of Size of a Locus",x="Size of a Locus",y="Proportion")
 data3 <- data.frame(EffectV = EffectV)
 p3 <- ggplot(data3,aes(EffectV))+
   geom_histogram(aes(x=data3$EffectV,y=(..count..)/sum(..count..)),
                  fill="gray15",
                  alpha=0.75
   )+
   theme_bw()+
   theme_new()+
   labs(title="The Histogram of Size of a Locus",x="Size of a Locus",y="Proportion")
 
 return (list(result,p1,p2,p3))
}

