##
##
## Average over MAF
## Function of EV and J
##
##
source('skat_calculator2.r')
source('skat_calculator2sc.r')
source('sampler.r')

theme_new <- function(base_size = 12, base_family = "Helvetica"){
  theme_bw(base_size = base_size, base_family = base_family) %+replace%
    theme(plot.title=element_text(size=12,face="bold.italic",hjust = 0.5),
          axis.title.x = element_text(size = 10,face = "bold",hjust = 0.5),
          axis.title.y =  element_text(size = 10,face = "bold",hjust = 0.5),
          #panel.background = element_blank(),
          #axis.line = element_line(colour = "black")
    )
}




# QT indicator for quantative trait

#
# E  can be a vector
#


#
# Delete nameEsseble
#
get_Aprox2 = function(EV,alpha,Total,CASE=NULL,CONTROL=NULL,PC=NULL,QT='CC',nameEsseble=NULL,JJ=NULL){
PowerEindP = NULL
PowerBindP = NULL
PowerRelBandP = NULL
n = CASE*CONTROL/Total
if (QT == 'QT'){
n = Total # sample size
}
ptm <- proc.time()
Obj = transform(Total) #Get appropriate MAF distribution of J and MAF


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

LL = length(EV)

for (j in 1:LL){
E  = EV[j]

EindP =NULL
BindP = NULL
ErelP = NULL

for (i in 1:1000){
#cat(i,'\n')

if (is.null(JJ)){
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

if (is.null(PC)){
pEinP = calcPower(pj,n*E,n,alpha)
pBinP = calcPowerBPind(pj,n*E,n,alpha)
ErelP = c(ErelP,calcPowerD(pj,n*E,n,gm,alpha))
}else{
SC = ceiling(J*PC)

pEinP = calcPower_SC(pj,SC,n*E,n,alpha)
pBinP = calcPowerBPind_SC(pj,SC,n*E,n,alpha)
ErelP = c(ErelP,calcPowerD_SC(pj,SC,n*E,n,gm,alpha))
}


BindP = c(BindP,pBinP)
EindP = c(EindP,pEinP)
}
PowerRelBandP = rbind(PowerRelBandP,ErelP)
PowerBindP = rbind(PowerBindP,BindP)
PowerEindP = rbind(PowerEindP,EindP)
}
#cat(sdd,'\n')
ptm1 <- proc.time()
cat('CPU time Total:',ptm1-ptm,'\n')
BrelMean <- apply(PowerRelBandP,1,mean)
BinMean <- apply(PowerBindP,1,mean)
EinMean <- apply(PowerEindP,1,mean)

EV.Length <- length(EV)

EV_x <- rep(EV*100,3)
group <- c(rep("Genetic Arc I",EV.Length),rep("Genetic Arc II",EV.Length),
           rep("Genetic Arc III",EV.Length))
Power <- c(BinMean,EinMean,BrelMean)
data <- data.frame(EV_x,group,Power)
p <- ggplot(data)+geom_line(aes(x=EV_x,y=Power,colour=group))+
  #scale_x_continuous(limits=limits.x,expand = c(0,0),breaks = breaks.x)+
  scale_y_continuous(expand = c(0,0),limits=c(0,1))+
  #scale_y_continuous(expand = c(0,0),limits=c(0,(density.y.max+0.5)))+
  theme_new()+
  labs(title="The Mean Power Distribution",y="Mean Power",x="Variance Explained (Percent)")

return(list(a=BrelMean,b=BinMean,c =EinMean,p))
}


transform=function(Total,MAFU=0.01){
ptm <- proc.time()

if (Total<5000){
cat(1,'\n')
A = readRDS('whole_N4375.rds')
}else{
if (Total<8000){
cat(2,'\n')
A = readRDS('whole_N8750.rds')
}else{
if (Total<16000){
cat(3,'\n')
A = readRDS('whole_N17500.rds')
}else{
if (Total<32000){
cat(4,'\n')
A = readRDS('whole_N35000.rds')
}else{
A = readRDS('whole.rds')
}
}
}
}
mafV = A[,2]
namesV = A[,1]

mafV = as.numeric(mafV)
ptm1 <- proc.time()
cat('CPU time load:',ptm1-ptm,'\n')
MAFL = 1/Total
KKK = (mafV >= MAFL) & (mafV <= MAFU) 
ptm <- proc.time()
nn = table(namesV[KKK])
names = names(nn)
SizeGene = table(as.vector(nn))
names = names[SizeGene>1]
SizeGene = SizeGene[SizeGene>1]
ptm1 <- proc.time()
cat('CPU time create tables:',ptm1-ptm,'\n')
return(list(SizeGene=SizeGene,mean=mean(mafV[KKK]),Names=names,SizeGeneV = nn,maf = mafV[KKK]))
}


