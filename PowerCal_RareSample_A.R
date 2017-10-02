##
##
## Average over MAF
## Function of EV and J
##
##

source('skat_sample_calculatorsc.r')
source('sampler.r')


###ceil a funciton up to a 10^places digit
ceilupto <- function(number,places=2){
  (ceiling(number/10^places))*10^places
}

theme_new <- function(base_size = 12, base_family = "Helvetica"){
  theme(
    plot.title=element_text(size=12,face="bold.italic",hjust = 0.5),
    axis.title.x = element_text(size = 10,face = "bold",hjust = 0.5),
    axis.title.y =  element_text(size = 10,face = "bold",hjust = 0.5)
    
    #panel.background = element_blank(),
    #axis.line = element_line(colour = "black")
  )
}


# QT indicator for quantative trait


get_Aprox_Sample <- function(EV,PowerThr,alpha,PC=NA,TEST = 'SKAT',QT='CC',nameEsseble=NULL,JJ=NA,PRC=NA,ONESNP=F){
  Ratio = 2
  Total=10^7
  CASE=10^7
  CONTROL=1
  JJ.vec <- rep(0,1000)
  if (!is.na(JJ)){
    if (JJ == 1){ONESNP=T}
  }
  if(ONESNP==T){
    JJ=1
    Total=10^3
    CASE=10^3
    if(length(EV)==1&TEST!='Burden Test'){
      PowerEindP = NULL
      PowerBindP = NULL
      PowerRelBandP = NULL
      n = CASE*CONTROL/Total
      n=1
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
        #save(PowerRelBandP,PowerBindP,PowerEindP,file='1')
      }
      
      mmEinP = ceilupto(mean(PowerEindP))
      qunatEinP = ceilupto(quantile(PowerEindP))[2:4] #25%, median, 75%
      
      
      mmBinP = ceilupto(mean(PowerBindP))
      qunatBinP = ceilupto(quantile(PowerBindP))[2:4] #25%, median, 75%
      
      mmBrelP = ceilupto(mean(PowerRelBandP))
      
      qunatBrelP =  ceilupto(quantile(PowerRelBandP))[2:4]
      
      data <- data.frame(PowerEindP=t(PowerEindP),
                         PowerBindP=t(PowerBindP),
                         PowerRelBandP=t(PowerRelBandP),
                         JJ.vec = JJ.vec)
      colnames(data) <- c("EindP","BindP","ErelP","JJ")
      
      
      p1 <- ggplot(data,aes(data$EindP))+
        geom_histogram(aes(x=data$EindP,y=(..count..)/sum(..count..)),
                       fill="#c0392b",
                       alpha=0.75
        )+
        theme_bw()+
        theme_new()+
        labs(title="The Histogram of Sample Size of Senarario S1",x="Sample Size",y="Proportion")
      p2 <- ggplot(data,aes(data$BindP))+
        geom_histogram(aes(x=data$BindP,y=(..count..)/sum(..count..)),
                       fill="dodgerblue4",
                       alpha=0.75
        )+
        theme_bw()+
        theme_new()+
        labs(title="The Histogram of Sample Size of Senarario S2",x="Sample Size",y="Proportion")
      p3 <- ggplot(data,aes(data$ErelP))+
        geom_histogram(aes(x=data$ErelP,y=(..count..)/sum(..count..)),
                       fill="chartreuse4",
                       alpha=0.75
        )+
        theme_bw()+
        theme_new()+
        labs(title="The Histogram of Sample Size of Senarario S3",x="Sample Size",y="Proportion")
      p4 <- ggplot(data,aes(data$JJ))+
        geom_histogram(aes(x=data$JJ,y=(..count..)/sum(..count..)),
                       fill="gray15",
                       alpha=0.75
        )+
        theme_bw()+
        theme_new()+
        labs(title="The Histogram of Total number of variants (J)",x="Total number of variants (J)",y="Proportion")
      
      Gene_Arc_1 <- c(mmEinP,qunatEinP)
      Gene_Arc_2 <- c(mmBinP,qunatBinP)
      Gene_Arc_3 <- c(mmBrelP,qunatBrelP)
      
      Power_Dist <- c("Mean Sample Size","25% Quantile of Sample Size","Medium of Sample Size",
                      "75% Quantile of Sample Size")
      combind.result <- data.frame(Power_Dist,Gene_Arc_1,Gene_Arc_2,Gene_Arc_3)
      colnames(combind.result) <- c("Sample Size Distribution","Senarario S1","Senarario S2",
                                    "Senarario S3")
      
      
      #return (list(combind.result,p1=NULL,p2=NULL,p3=NULL,p4=p4))
      return (list(combind.result,p1=p1,p2=p2,p3=p3,p4=p4))
      
    }
    else if(length(EV)>1&TEST!='Burden Test'){
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
            SC = ceiling(J*PC)
            
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
        #save(PowerRelBandP,PowerBindP,PowerEindP,file='2')
      }
      BrelMean <- ceilupto(apply(PowerRelBandP,1,mean))
      BinMean <- ceilupto(apply(PowerBindP,1,mean))
      EinMean <- ceilupto(apply(PowerEindP,1,mean))
      
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
        labs(title="The Mean Sample Size Distribution of Senarario 1",y="Mean Sample Size",x="Variance Explained (Percent)",y="Proportion")
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
        geom_histogram(aes(x=data$JJ,y=(..count..)/sum(..count..)),
                       fill="gray15",
                       alpha=0.75
        )+
        theme_bw()+
        theme_new()+
        labs(title="The Histogram of total number of variants (J)",x="total number of variants (J)", y = "Proportion")
      
      
      
      
      MeanPower.combine <- data.frame(EV=EV*100,GeneI =EinMean,GeneII = BinMean,GeneIII=BrelMean)
      colnames(MeanPower.combine) <- c("EV(Percent)",
                                       "Senarario 1 Mean Sample Size",
                                       "Senarario 2 Mean Sample Size",
                                       "Senarario 3 Mean Sample Size"
      )
      MeanPower.combine[,2:4] <- round(MeanPower.combine[,2:4],3)
      MeanPower.combine <- MeanPower.combine[c(1,3,5,7,9),]
      return(list(MeanPower.combine,p1=p1,p2=p2,p3=p3,p4=p4))
    }
    
  }else{
    if(length(EV)==1&TEST!='Burden Test'){
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
      n=1
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
      LL = length(EV)
      
      
      
      for (j in 1:LL){
        E  = EV[j]
        #cat('E',E,'\n')
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
        #save(PowerRelBandP,PowerBindP,PowerEindP,file='4')
      }
      #cat(sdd,'\n')
      #ptm1 <- proc.time()
      #cat('CPU time Total:',ptm1-ptm,'\n')
      
      # mmEinP = round(mean(PowerEindP),4)
      # qunatEinP = round(quantile(PowerEindP),4)[2:4] #25%, median, 75%
      # 
      # mmBinP = round(mean(PowerBindP),4)
      # qunatBinP = round(quantile(PowerBindP),4)[2:4] #25%, median, 75%
      # 
      # mmBrelP = round(mean(PowerRelBandP),4)
      
      # qunatBrelP =  round(quantile(PowerRelBandP),4)[2:4]
      mmEinP = ceilupto(mean(PowerEindP))
      qunatEinP = ceilupto(quantile(PowerEindP))[2:4] #25%, median, 75%
      
      #mmBinP = round(mean(PowerBindP),4)
      mmBinP = ceilupto(mean(PowerBindP))
      qunatBinP = ceilupto(quantile(PowerBindP))[2:4] #25%, median, 75%
      
      mmBrelP = ceilupto(mean(PowerRelBandP))
      
      qunatBrelP =  ceilupto(quantile(PowerRelBandP))[2:4]
      
      # Power <- c(PowerEindP,PowerBindP,PowerRelBandP)
      # Method <- c(rep("Senarario S1",length(PowerEindP)),rep("Senarario S2",length(PowerBindP)),rep("Senarario S3",length(PowerRelBandP)))
      # data <- data.frame(Power,Method)
      # 
      data <- data.frame(PowerEindP=t(PowerEindP),
                         PowerBindP=t(PowerBindP),
                         PowerRelBandP=t(PowerRelBandP),
                         JJ.vec = JJ.vec)
      colnames(data) <- c("EindP","BindP","ErelP","JJ")
      
      
      p1 <- ggplot(data,aes(data$EindP))+
        geom_histogram(aes(x=data$EindP,y=..count../(sum(..count..))),
                       fill="#c0392b",
                       alpha=0.75
        )+
        theme_bw()+
        theme_new()+
        labs(title="The Histogram of Sample Size of Senarario S1",x="Sample Size",y="Proportion")
      p2 <- ggplot(data,aes(data$BindP))+
        geom_histogram(aes(x=data$BindP,y=..count../(sum(..count..))),
                       fill="dodgerblue4",
                       alpha=0.75
        )+
        theme_bw()+
        theme_new()+
        labs(title="The Histogram of Sample Size of Senarario S2",x="Sample Size",y="Proportion")
      p3 <- ggplot(data,aes(data$ErelP))+
        geom_histogram(aes(x=data$ErelP,y=..count../(sum(..count..))),
                       fill="chartreuse4",
                       alpha=0.75
        )+
        theme_bw()+
        theme_new()+
        labs(title="The Histogram of Sample Size of Senarario S3",x="Sample Size",y="Proportion")
      p4 <- ggplot(data,aes(data$JJ))+
        geom_histogram(aes(x=data$JJ,y=..count../(sum(..count..))),
                       fill="gray15",
                       alpha=0.75
        )+
        theme_bw()+
        theme_new()+
        labs(title="The Histogram of Total number of variants (J)",x="Total number of variants (J)",y="Proportion")
      
      
      # p <- ggplot(data,aes(x=Power))+ geom_histogram(aes(group=Method,colour=Method,fill=Method),alpha=0.3,adjust=1)+
      #scale_x_continuous(limits=limits.x,expand = c(0,0),breaks = breaks.x)+
      #scale_y_continuous(expand = c(0,0),limits=c(0,(density.y.max+0.5)))+
      
      #theme_new()+
      #labs(title="The Histogram of Sample Size",x="Sample Size",y="Proportion")
      # p <- ggplot(data, aes(x=Power)) +
      #   geom_histogram(aes(y=..count../(sum(..count..)),fill=Method,colour=Method,group=Method),binwidth=.02, alpha=.3, position="identity")+
      #   scale_x_continuous(limits=c(0,1.02),expand = c(0,0),breaks = seq(0,1.02,0.1))+
      #   #scale_y_continuous(expand = c(0,0),limits=c(0,15))+
      #   theme_new()+
      #   labs(title="The Density of Power",y="Density")
      
      Gene_Arc_1 <- c(mmEinP,qunatEinP)
      Gene_Arc_2 <- c(mmBinP,qunatBinP)
      Gene_Arc_3 <- c(mmBrelP,qunatBrelP)
      
      Power_Dist <- c("Mean Sample Size","25% Quantile of Sample Size","Medium of Sample Size",
                      "75% Quantile of Sample Size")
      combind.result <- data.frame(Power_Dist,Gene_Arc_1,Gene_Arc_2,Gene_Arc_3)
      colnames(combind.result) <- c("Sample Size Distribution","Senarario S1","Senarario S2",
                                    "Senarario S3")
      
      
      return (list(combind.result,p1=p1,p2=p2,p3=p3,p4=p4))
      #return (list(combind.result,grid.arrange(p1,p2,p3,p4,ncol=2)))
    }
    else if(length(EV)>1&TEST!='Burden Test'){
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
      # group <- c(rep("Senarario S1",EV.Length),rep("Senarario S2",EV.Length),
      #            rep("Senarario S3",EV.Length))
      # Power <- c(EinMean,BinMean,BrelMean)
      # data <- data.frame(EV_x,group,Power)
      BrelMean <- ceilupto(apply(PowerRelBandP,1,mean))
      BinMean <- ceilupto(apply(PowerBindP,1,mean))
      EinMean <- ceilupto(apply(PowerEindP,1,mean))
      
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
        labs(title="The Mean Sample Size Distribution of Senarario 1",y="Mean Sample Size",x="Variance Explained (Percent)",y="Proportion")
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
        geom_histogram(aes(x=data$JJ,y=(..count..)/sum(..count..)),
                       fill="gray15",
                       alpha=0.75
        )+
        theme_bw()+
        theme_new()+
        labs(title="The Histogram of total number of variants (J)",x="total number of variants (J)",y="Proportion")
      
      
      
      
      MeanPower.combine <- data.frame(EV=EV*100,GeneI =EinMean,GeneII = BinMean,GeneIII=BrelMean)
      colnames(MeanPower.combine) <- c("EV(Percent)",
                                       "Senarario 1 Mean Sample Size",
                                       "Senarario 2 Mean Sample Size",
                                       "Senarario 3 Mean Sample Size"
      )
      MeanPower.combine[,2:4] <- round(MeanPower.combine[,2:4],3)
      MeanPower.combine <- MeanPower.combine[c(1,3,5,7,9),]
      return(list(MeanPower.combine,p1=p1,p2=p2,p3=p3,p4=p4))
    }
    else if(TEST=='Burden Test'&length(EV)==1) {
      
      Power= NULL
      
      n = CASE*CONTROL/Total
      if (QT == 'QT'){
        n = 1 # sample size
        Ratio=1
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
      LL = length(EV)
      
      
      
      for (j in 1:LL){
        E  = EV[j]
        
        Pow =NULL
        
        
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
          SC = ceiling(J*PC)
          SP = ceiling(PRC*SC)
          pEP = calcPowerLSample(J,PC,PRC,E,PowerThr,alpha,Ratio)
          
          Pow = c(Pow,pEP)
        }
        PowerX = Pow[Pow!=Ratio*10^4]
        Pow[Pow==Ratio*10^4]=mean(PowerX)
        Power = rbind(Power,Pow)
      }
      
      PowerBindP <- PowerEindP <- PowerRelBandP <- Power
      save(PowerRelBandP,PowerBindP,PowerEindP,file='6')
      mmEinP = ceilupto(mean(PowerEindP))
      qunatEinP = ceilupto(quantile(PowerEindP))[2:4] #25%, median, 75%
      
      
      mmBinP = ceilupto(mean(PowerBindP))
      qunatBinP = ceilupto(quantile(PowerBindP))[2:4] #25%, median, 75%
      
      mmBrelP = ceilupto(mean(PowerRelBandP))
      
      qunatBrelP =  ceilupto(quantile(PowerRelBandP))[2:4]
      
      data <- data.frame(PowerEindP=t(PowerEindP),
                         PowerBindP=t(PowerBindP),
                         PowerRelBandP=t(PowerRelBandP),
                         JJ.vec = JJ.vec)
      colnames(data) <- c("EindP","BindP","ErelP","JJ")
      
      
      p1 <- ggplot(data,aes(data$EindP))+
        geom_histogram(aes(x=data$EindP,y=..count../(sum(..count..))),
                       fill="#c0392b",
                       alpha=0.75
        )+
        theme_bw()+
        theme_new()+
        labs(title="The Histogram of Sample Size of Senarario S1",x="Sample Size",y="Proportion")
      p2 <- ggplot(data,aes(data$BindP))+
        geom_histogram(aes(x=data$BindP,y=..count../(sum(..count..))),
                       fill="dodgerblue4",
                       alpha=0.75
        )+
        theme_bw()+
        theme_new()+
        labs(title="The Histogram of Sample Size of Senarario S2",x="Sample Size",y="Proportion")
      p3 <- ggplot(data,aes(data$ErelP))+
        geom_histogram(aes(x=data$ErelP,y=..count../(sum(..count..))),
                       fill="chartreuse4",
                       alpha=0.75
        )+
        theme_bw()+
        theme_new()+
        labs(title="The Histogram of Sample Size of Senarario S3",x="Sample Size",y="Proportion")
      p4 <- ggplot(data,aes(data$JJ))+
        geom_histogram(aes(x=data$JJ,y=..count../(sum(..count..))),
                       fill="gray15",
                       alpha=0.75
        )+
        theme_bw()+
        theme_new()+
        labs(title="The Histogram of Total number of variants (J)",x="Total number of variants (J)",y="Proportion")
      
      Gene_Arc_1 <- c(mmEinP,qunatEinP)
      Gene_Arc_2 <- c(mmBinP,qunatBinP)
      Gene_Arc_3 <- c(mmBrelP,qunatBrelP)
      
      Power_Dist <- c("Mean Sample Size","25% Quantile of Sample Size","Medium of Sample Size",
                      "75% Quantile of Sample Size")
      combind.result <- data.frame(Power_Dist,Gene_Arc_1,Gene_Arc_2,Gene_Arc_3)
      colnames(combind.result) <- c("Sample Size Distribution","Senarario S1","Senarario S2",
                                    "Senarario S3")
      
      
      return (list(combind.result,p1=p1,p2=p2,p3=p3,p4=p4))
    }else if(TEST=='Burden Test'&(length(EV)>1)){
      Power= NULL
      
      n = CASE*CONTROL/Total
      if (QT == 'QT'){
        n = 1 # sample size
        Ratio=1
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
      LL = length(EV)
      
      
      
      for (j in 1:LL){
        E  = EV[j]
        
        Pow =NULL
        
        
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
          SC = ceiling(J*PC)
          SP = ceiling(PRC*SC)
          pEP = calcPowerLSample(J,SC,SP,E,PowerThr,alpha,Ratio)
          
          Pow = c(Pow,pEP)
          
        }
        PowerX = Pow[Pow!=Ratio*10^4]
        Pow[Pow==Ratio*10^4]=mean(PowerX)
        Power = rbind(Power,Pow)
        Power = rbind(Power,Pow)
        
      }
      
      PowerRelBandP <- PowerBindP <- PowerEindP <- Power
      save(PowerRelBandP,PowerBindP,PowerEindP,file='6')
      BrelMean <- apply(PowerRelBandP,1,mean)
      BinMean <- apply(PowerBindP,1,mean)
      EinMean <- apply(PowerEindP,1,mean)
      
      EV.Length <- length(EV)
      
      # EV_x <- rep(EV*100,3)
      # group <- c(rep("Senarario S1",EV.Length),rep("Senarario S2",EV.Length),
      #            rep("Senarario S3",EV.Length))
      # Power <- c(EinMean,BinMean,BrelMean)
      # data <- data.frame(EV_x,group,Power)
      data <- data.frame(EV*100,EinMean,BinMean,BrelMean)
      colnames(data) <- c("EV","EindP","BindP","ErelP")
      #y.up <- min(1,(max(Power)+0.1))
      p1 <- ggplot(data)+geom_line(aes(x=EV,y=EindP), colour="#c0392b")+
        theme_bw()+
        theme_new()+
        labs(title="The Mean Sample Size Distribution of Senarario S1",y="Mean Sample Size",x="Variance Explained (Percent)")
      p2 <- ggplot(data)+geom_line(aes(x=EV,y=BindP), colour="dodgerblue4")+
        theme_bw()+
        theme_new()+
        labs(title="The Mean Sample Size Distribution of Senarario S2",y="Mean Sample Size",x="Variance Explained (Percent)")
      p3 <- ggplot(data)+geom_line(aes(x=EV,y=ErelP), colour="chartreuse4")+
        theme_bw()+
        theme_new()+
        labs(title="The Mean Sample Size Distribution of Senarario S3",y="Mean Sample Size",x="Variance Explained (Percent)")
      data <- data.frame(JJ=JJ.vec)
      p4 <- ggplot(data,aes(data$JJ))+
        geom_histogram(aes(x=data$JJ,y=..count../(sum(..count..))),
                       fill="gray15",
                       alpha=0.75
        )+
        theme_bw()+
        theme_new()+
        labs(title="The Histogram of Total number of variants (J)",x="Total number of variants (J)",y="Proportion")
      
      
      
      
      MeanPower.combine <- data.frame(EV=EV*100,GeneI =EinMean,GeneII = BinMean,GeneIII=BrelMean)
      colnames(MeanPower.combine) <- c("EV(Percent)",
                                       "Senarario S1 Mean Sample Size",
                                       "Senarario S2 Mean Sample Size",
                                       "Senarario S3 Mean Sample Size"
      )
      MeanPower.combine[,2:4] <- round(MeanPower.combine[,2:4],3)
      MeanPower.combine <- MeanPower.combine[c(1,3,5,7,9),]
      return(list(MeanPower.combine,p1=p1,p2=p2,p3=p3,p4=p4))
    }
  }
}

calcPowerL =function(J,SC,SP,E,n,alpha){
  nc1 = n*E*abs(SC-2*SP)^2/J/SC
  ctt = qchisq(1-alpha,df=1)
  power = 1 - pchisq(ctt,df=1,nc=nc1)
  return (power) 
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

