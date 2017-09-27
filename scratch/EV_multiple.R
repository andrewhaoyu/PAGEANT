PowerEindP = NULL
PowerBindP = NULL
PowerRelBandP = NULL
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
  n = 1 # sample size
  Ratio=1
}
# ptm <- proc.time()
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
    JJ.vec[i] <- J
    pj = as.numeric(sample(MAF,J,replace=TRUE))
    #pj =runif(J)*0.1
    m = E*n/J/(a+g*mean(pj))
    gm = g*m
    
    if (is.na(PC)){
      pEinP = calcPowerSample(pj,1,n*E,PowerThr,aa,bb,alpha,Ratio)
      pBinP = calcPowerBPindSample(pj,1,n*E,PowerThr,aa,bb,alpha,Ratio)
      ErelP = c(ErelP,calcPowerDSample(pj,1,n*E,PowerThr,gm,aa,bb,alpha,Ratio))
    }else{
      
      pEinP = calcPowerSample(pj,PC,n*E,PowerThr,aa,bb,alpha,Ratio)
      pBinP = calcPowerBPindSample(pj,PC,n*E,PowerThr,aa,bb,alpha,Ratio)
      ErelP = c(ErelP,calcPowerDSample(pj,PC,n*E,PowerThr,gm,aa,bb,alpha,Ratio))
    }
    
    BindP = c(BindP,pBinP)
    EindP = c(EindP,pEinP)
  }
  ErelPX = ErelP[ErelP!=Ratio*10^4]
  ErelP[ErelP==Ratio*10^4]=mean(ErelPX)
  ErelPX = EindP[EindP!=Ratio*10^4]
  EindP[EindP==Ratio*10^4]=mean(ErelPX)
  ErelPX = BindP[BindP!=Ratio*10^4]
  BindP[BindP==Ratio*10^4]=mean(ErelPX)
  PowerRelBandP = rbind(PowerRelBandP,ErelP)
  PowerBindP = rbind(PowerBindP,BindP)
  PowerEindP = rbind(PowerEindP,EindP)
  save(PowerRelBandP,PowerBindP,PowerEindP,file='5')
}
#cat(sdd,'\n')
#ptm1 <- proc.time()
#cat('CPU time Total:',ptm1-ptm,'\n')
BrelMean <- apply(PowerRelBandP,1,mean)
BinMean <- apply(PowerBindP,1,mean)
EinMean <- apply(PowerEindP,1,mean)

EV.Length <- length(EV)

# EV_x <- rep(EV*100,3)
# group <- c(rep("Senarario 1",EV.Length),rep("Senarario 2",EV.Length),
#            rep("Senarario 3",EV.Length))
# Power <- c(EinMean,BinMean,BrelMean)
# data <- data.frame(EV_x,group,Power)
data <- data.frame(EV*100,EinMean,BinMean,BrelMean)
colnames(data) <- c("EV","EindP","BindP","ErelP")
#y.up <- min(1,(max(Power)+0.1))
p1 <- ggplot(data)+geom_line(aes(x=EV,y=EindP), colour="#c0392b")+
  theme_bw()+
  theme_new()+
  labs(title="The Mean Sample Size Distribution of Senarario 1",y="Mean Sample Size",x="Variance Explained (Percent)")
p2 <- ggplot(data)+geom_line(aes(x=EV,y=BindP), colour="dodgerblue4")+
  theme_bw()+
  theme_new()+
  labs(title="The Mean Sample Size Distribution of Senarario 2",y="Mean Sample Size",x="Variance Explained (Percent)")
p3 <- ggplot(data)+geom_line(aes(x=EV,y=ErelP), colour="chartreuse4")+
  theme_bw()+
  theme_new()+
  labs(title="The Mean Sample Size Distribution of Senarario 3",y="Mean Sample Size",x="Variance Explained (Percent)")
data <- data.frame(JJ=JJ.vec)
p4 <- ggplot(data,aes(data$JJ))+
  geom_histogram(aes(x=data$JJ,y=..density..),
                 fill="gray15",
                 alpha=0.75
  )+
  theme_bw()+
  theme_new()+
  labs(title="The Histogram of Parameters",x="Parameters")




MeanPower.combine <- data.frame(EV=EV*100,GeneI =EinMean,GeneII = BinMean,GeneIII=BrelMean)
colnames(MeanPower.combine) <- c("EV(Percent)",
                                 "Senarario 1 Mean Sample Size",
                                 "Senarario 2 Mean Sample Size",
                                 "Senarario 3 Mean Sample Size"
)
MeanPower.combine[,2:4] <- round(MeanPower.combine[,2:4],3)
MeanPower.combine <- MeanPower.combine[c(1,3,5,7,9),]
return(list(MeanPower.combine,p1=p1,p2=p2,p3=p3,p4=p4))