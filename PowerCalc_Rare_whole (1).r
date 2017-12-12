get_AproxQ2  = function(al,be,K,alpha,Total,Obj,perm,CASE=NULL,CONTROL=NULL,PC=NULL,TEST = 'SKAT',QT='CC',nameEsseble=NULL,JJ=NULL){
PowerEindP = NULL
PowerBindP = NULL
PowerRelBandP = NULL
EE = NULL
n = CASE*CONTROL/Total

if (TEST=='SKAT'){
aa=1
bb=25
}
if (TEST=='Calpha'){
aa=1
bb=1
}
if (TEST=='Hotelling'){
aa=1/2
bb=1/2
}
#cat(aa,bb,'\n')


if (QT == 'QT'){
n = Total # sample size
}
ptm <- proc.time()
#Obj = transform(Total) #Get appropriate MAF distribution of J and MAF


MAFL = 1/Total
MAFU = 0.01
g = (2*MAFU*(1-MAFU)*log10(MAFU)^2 - 2*MAFL*(1-MAFL)*log10(MAFL)^2)/(MAFU - MAFL)
a = 2*MAFU*(1-MAFU)*log10(MAFU)^2 - MAFU*g



VPNames = Obj$Names
TJ = Obj$SizeGene
TJJ = Obj$SizeGeneV
DJ = dist(TJ)
meanP = Obj$mean
MAF = Obj$maf

EindP =NULL
BindP = NULL
ErelP = NULL

for (i in 1:perm){
#cat(i,'\n')

EVV = 0.51
while (sum(EVV)>0.5){
#EVV = rtrunc(K*2,'gamma',0,0.01,al,be) 
EVV = rgamma(K+100,al,be)
EVV = EVV[EVV<0.01]
EVV = EVV[1:K]
#cat(length(EVV),'\n')
}
#cat(sum(EVV),al,be,'\n')
E = EVV[1]
EE = c(EE,E)
#cat(E,'\n')
if (is.na(JJ)){
if (is.null(nameEsseble)){
#
# Draw sample
#
#cat('Gene is not  Present','\n')
J = as.numeric(sampleD(DJ,names(TJ),1))

}else{
kkk = which(VPNames==nameEsseble)
if (length(kkk)==0){
#
# Draw sample
#
#cat('Gene is not  Present','\n')
J = as.numeric(sampleD(DJ,names(TJ),1))
}else{
#cat('Gene is Present','\n')
J = TJJ[kkk]
#cat(J,'\n')
}
}
}else{
#cat('JJ',JJ,'\n')
J = JJ
}
#cat(J,'\n')
#
# Draw MAF
#

pj = as.numeric(sample(MAF,J,replace=TRUE))
#pj =runif(J)*0.1
m = E*n/J/(a+g*mean(pj))
gm = g*m

if (is.na(PC)){
pEinP = calcPower(pj,n*E,n,aa,bb,alpha)
pBinP = calcPowerBPind(pj,n*E,n,aa,bb,alpha)
ErelP = c(ErelP,calcPowerD(pj,n*E,n,gm,aa,bb,alpha))
}else{
SC = ceiling(J*PC)

pEinP = calcPower_SC(pj,SC,n*E,n,aa,bb,alpha)
pBinP = calcPowerBPind_SC(pj,SC,n*E,n,aa,bb,alpha)
ErelP = c(ErelP,calcPowerD_SC(pj,SC,n*E,n,gm,aa,bb,alpha))
}


BindP = c(BindP,pBinP)
EindP = c(EindP,pEinP)
}
PowerRelBandP = rbind(PowerRelBandP,ErelP)
PowerBindP = rbind(PowerBindP,BindP)
PowerEindP = rbind(PowerEindP,EindP)
#cat(sdd,'\n')

ptm1 <- proc.time()
cat('CPU time Total:',ptm1-ptm,'\n')
return (list(a=PowerEindP,b=PowerBindP,c =PowerRelBandP,ee=EE)) 
}


