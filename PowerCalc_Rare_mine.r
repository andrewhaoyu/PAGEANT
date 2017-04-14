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
          axis.title.x = element_text(size = 10,face = "bold",hjust=0.5),
          axis.title.y =  element_text(size = 10,face = "bold",hjust=0.5),
          #panel.background = element_blank(),
          #axis.line = element_line(colour = "black")
    )
}

# QT indicator for quantative trait


#
# Delete nameEsseble
#
get_Aprox = function(E,alpha,Total,CASE=NULL,CONTROL=NULL,PC=NULL,QT='CC',nameEsseble=NULL,JJ=NULL){
PowerEindP = NULL
PowerBindP = NULL
PowerRelBandP = NULL
n = CASE*CONTROL/Total
if (QT == 'QT'){
n = Total # sample size
}
#ptm <- proc.time()
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
PowerRelBandP = c(PowerRelBandP,calcPowerD(pj,n*E,n,gm,alpha))
}else{
SC = ceiling(J*PC)

pEinP = calcPower_SC(pj,SC,n*E,n,alpha)
pBinP = calcPowerBPind_SC(pj,SC,n*E,n,alpha)
PowerRelBandP = c(PowerRelBandP,calcPowerD_SC(pj,SC,n*E,n,gm,alpha))
}


PowerBindP = c(PowerBindP,pBinP)
PowerEindP = c(PowerEindP,pEinP)
}
mmEinP = round(mean(PowerEindP),4)
qunatEinP = round(quantile(PowerEindP),4)[2:4] #25%, median, 75%

mmBinP = round(mean(PowerBindP),4)
qunatBinP = round(quantile(PowerBindP),4)[2:4] #25%, median, 75%

mmBrelP = round(mean(PowerRelBandP),4)

qunatBrelP =  round(quantile(PowerRelBandP),4)[2:4]

density.Ein <- density(PowerEindP)
density.Bin <- density(PowerBindP)
density.Brel <- density(PowerRelBandP)

density.y.max <- max(c(quantile(density.Ein$y),quantile(density.Bin$y),quantile(density.Brel$y)))
density.x.min <- min(c(qunatEinP[1],qunatBinP[1],qunatBrelP[1]))
density.x.max <- max(c(qunatEinP[3],qunatBinP[3],qunatBrelP[3]))

limits.x <- c(0,1)
breaks.x <- seq(0,1,0.1)
if((density.x.max-density.x.min)<0.1)
{
  limits.x.low <- round(density.x.min,1)
  limits.x.high <- round(density.x.max,1)
  if(limits.x.low==limits.x.high){
    limits.x.low <- limits.x.low-0.01
  }
  limits.x <- c(limits.x.low,limits.x.high)
  breaks.x <- seq(limits.x.low,limits.x.high,0.001)
}

Power <- c(PowerBindP,PowerEindP,PowerRelBandP)
Method <- c(rep("Geneic Arc I",length(PowerBindP)),rep("Genetic Arc II",length(PowerEindP)),rep("Genetic Arc III",length(PowerRelBandP)))
data <- data.frame(Power,Method)

p <- ggplot(data,aes(x=Power))+
  geom_density(aes(group=Method,colour=Method,fill=Method),alpha=0.3,adjust=1)+
  scale_x_continuous(limits=limits.x,expand = c(0,0),breaks = breaks.x)+
  scale_y_continuous(expand = c(0,0),limits=c(0,(density.y.max+0.5)))+
  #scale_y_continuous(expand = c(0,0),limits=c(0,(density.y.max+0.5)))+
  theme_new()+
  labs(title="The Density of Power",y="Density")
# p <- ggplot(data, aes(x=Power)) +
#   geom_histogram(aes(y=..density..,fill=Method,colour=Method,group=Method),binwidth=.02, alpha=.3, position="identity")+
#   scale_x_continuous(limits=c(0,1.02),expand = c(0,0),breaks = seq(0,1.02,0.1))+
#   #scale_y_continuous(expand = c(0,0),limits=c(0,15))+
#   theme_new()+
#   labs(title="The Density of Power",y="Density")
  
#cat(sdd,'\n')
#ptm1 <- proc.time()
#cat('CPU time Total:',ptm1-ptm,'\n')
return (list(quantile=c(meanPowerBinP=mmBinP,quantBinP = qunatBinP,meanPowerEinP = mmEinP,quantEinP=qunatEinP,meanPowerBrelP = mmBrelP,qunatBrelP = qunatBrelP),p=p,PowerEindP,PowerBindP,PowerRelBandP)) 
}




transform=function(Total,MAFU=0.01){
#ptm <- proc.time()

if (Total<5000){
#cat(1,'\n')
A = readRDS('whole_N4375.rds')
}else{
if (Total<8000){
#cat(2,'\n')
A = readRDS('whole_N8750.rds')
}else{
if (Total<16000){
#cat(3,'\n')
A = readRDS('whole_N17500.rds')
}else{
if (Total<32000){
#cat(4,'\n')
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
#ptm1 <- proc.time()
#cat('CPU time load:',ptm1-ptm,'\n')
MAFL = 1/Total
KKK = (mafV >= MAFL) & (mafV <= MAFU) 
#ptm <- proc.time()
nn = table(namesV[KKK])
names = names(nn)
SizeGene = table(as.vector(nn))
names = names[SizeGene>1]
SizeGene = SizeGene[SizeGene>1]
#ptm1 <- proc.time()
#cat('CPU time create tables:',ptm1-ptm,'\n')
return(list(SizeGene=SizeGene,mean=mean(mafV[KKK]),Names=names,SizeGeneV = nn,maf = mafV[KKK]))
}


