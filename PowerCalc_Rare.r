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
  theme(
    plot.title=element_text(size=12,face="bold.italic",hjust = 0.5),
    axis.title.x = element_text(size = 10,face = "bold",hjust = 0.5),
    axis.title.y =  element_text(size = 10,face = "bold",hjust = 0.5)
    
    #panel.background = element_blank(),
    #axis.line = element_line(colour = "black")
  )
}




get_Aprox <- function(EV,alpha,Total,CASE=NULL,CONTROL=NULL,PC=NA,TEST = 'SKAT',QT='CC',nameEsseble=NULL,JJ=NA,PRC=NA,ONESNP=F){
  JJ.vec <- rep(0,1000)
  if(ONESNP==T){
    JJ=1
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
      }
      #cat(sdd,'\n')
      #ptm1 <- proc.time()
      #cat('CPU time Total:',ptm1-ptm,'\n')
      
      mmEinP = round(mean(PowerEindP),4)
      qunatEinP = round(quantile(PowerEindP),4)[2:4] #25%, median, 75%
      
      mmBinP = round(mean(PowerBindP),4)
      qunatBinP = round(quantile(PowerBindP),4)[2:4] #25%, median, 75%
      
      mmBrelP = round(mean(PowerRelBandP),4)
      
      qunatBrelP =  round(quantile(PowerRelBandP),4)[2:4]
      
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
        labs(title="The Histogram of Power of Senarario S1",x="Power",y="Proportion")
      p2 <- ggplot(data,aes(data$BindP))+
        geom_histogram(aes(x=data$BindP,y=(..count..)/sum(..count..)),
                       fill="dodgerblue4",
                       alpha=0.75
        )+
        theme_bw()+
        theme_new()+
        labs(title="The Histogram of Power of Senarario S2",x="Power",y="Proportion")
      p3 <- ggplot(data,aes(data$ErelP))+
        geom_histogram(aes(x=data$ErelP,y=(..count..)/sum(..count..)),
                       fill="chartreuse4",
                       alpha=0.75
        )+
        theme_bw()+
        theme_new()+
        labs(title="The Histogram of Power of Senarario S3",x="Power",y="Proportion")
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
      
      Power_Dist <- c("Mean Power","25% Quantile of Power","Medium of Power",
                      "75% Quantile of Power")
      combind.result <- data.frame(Power_Dist,Gene_Arc_1,Gene_Arc_2,Gene_Arc_3)
      colnames(combind.result) <- c("Power Distribution","Senarario S1","Senarario S2",
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
        n = Total # sample size
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
      }
      #cat(sdd,'\n')
      #ptm1 <- proc.time()
      #cat('CPU time Total:',ptm1-ptm,'\n')
      BrelMean <- apply(PowerRelBandP,1,mean)
      BinMean <- apply(PowerBindP,1,mean)
      EinMean <- apply(PowerEindP,1,mean)
      
      EV.Length <- length(EV)
      
      EV_x <- rep(EV*100,3)
      group <- c(rep("Genetic Arc I",EV.Length),rep("Genetic Arc II",EV.Length),
                 rep("Genetic Arc III",EV.Length))
      Power <- c(EinMean,BinMean,BrelMean)
      data <- data.frame(EV_x,group,Power)
      y.up <- min(1,(max(Power)+0.1))
      p <- ggplot(data)+geom_line(aes(x=EV_x,y=Power,colour=group))+
        scale_y_continuous(expand = c(0,0),limits=c(0,y.up))+
        #scale_y_continuous(expand = c(0,0),limits=c(0,(density.y.max+0.5)))+
        theme_new()+
        labs(title="The Mean Power Distribution",y="Mean Power",x="Variance Explained (Percent)")
      MeanPower.combine <- data.frame(EV=EV*100,GeneI =EinMean,GeneII = BinMean,GeneIII=BrelMean)
      colnames(MeanPower.combine) <- c("EV(Percent)",
                                       "Genetic Arc I Mean Power",
                                       "Genetic Arc II Mean Power",
                                       "Genetic Arc III Mean Power"
      )
      MeanPower.combine[,2:4] <- round(MeanPower.combine[,2:4],3)
      MeanPower.combine <- MeanPower.combine[c(1,3,5,7,9),]
      return(list(MeanPower.combine,p))
    }
    else if(TEST=='Burden Test'&length(EV)==1) {
      
      Power= NULL
      
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
          pEP = calcPowerL(J,SC,SP,E,n,alpha)
          
          Pow = c(Pow,pEP)
        }
        Power = rbind(Power,Pow)
      }
      
      PowerBindP <- PowerEindP <- PowerRelBandP <- Power
      
      
      
      BrelMean <- apply(PowerRelBandP,1,mean)
      BinMean <- apply(PowerBindP,1,mean)
      EinMean <- apply(PowerEindP,1,mean)
      
      EV.Length <- length(EV)
      
      
      EV.Length <- length(EV)
      
      data <- data.frame(EV*100,EinMean,BinMean,BrelMean)
      colnames(data) <- c("EV","EindP","BindP","ErelP")
      #y.up <- min(1,(max(Power)+0.1))
      p1 <- ggplot(data)+geom_line(aes(x=EV,y=EindP), colour="#c0392b")+
        theme_bw()+
        theme_new()+
        labs(title="The Mean Power Distribution of Senarario 1",y="Mean Power",x="Variance Explained (Percent)",y="Proportion")
      p2 <- ggplot(data)+geom_line(aes(x=EV,y=BindP), colour="dodgerblue4")+
        theme_bw()+
        theme_new()+
        labs(title="The Mean Power Distribution of Senarario 2",y="Mean Power",x="Variance Explained (Percent)")
      p3 <- ggplot(data)+geom_line(aes(x=EV,y=ErelP), colour="chartreuse4")+
        theme_bw()+
        theme_new()+
        labs(title="The Mean Power Distribution of Senarario 3",y="Mean Power",x="Variance Explained (Percent)")
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
                                       "Senarario 1 Mean Power",
                                       "Senarario 2 Mean Power",
                                       "Senarario 3 Mean Power"
      )
      MeanPower.combine[,2:4] <- round(MeanPower.combine[,2:4],3)
      MeanPower.combine <- MeanPower.combine[c(1,3,5,7,9),]
      return(list(MeanPower.combine,p1=p1,p2=p2,p3=p3,p4=p4))
    }else if(TEST=='Burden Test'&(length(EV)>1)){
      Power= NULL
      
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
          pEP = calcPowerL(J,SC,SP,E,n,alpha)
          
          Pow = c(Pow,pEP)
        }
        Power = rbind(Power,Pow)
        
      }
      
      PowerRelBandP <- PowerBindP <- PowerEindP <- Power
      
      
      BrelMean <- apply(PowerRelBandP,1,mean)
      BinMean <- apply(PowerBindP,1,mean)
      EinMean <- apply(PowerEindP,1,mean)
      
      EV.Length <- length(EV)
      
      
      EV.Length <- length(EV)
      
      data <- data.frame(EV*100,EinMean,BinMean,BrelMean)
      colnames(data) <- c("EV","EindP","BindP","ErelP")
      #y.up <- min(1,(max(Power)+0.1))
      p1 <- ggplot(data)+geom_line(aes(x=EV,y=EindP), colour="#c0392b")+
        theme_bw()+
        theme_new()+
        labs(title="The Mean Power Distribution of Senarario 1",y="Mean Power",x="Variance Explained (Percent)",y="Proportion")
      p2 <- ggplot(data)+geom_line(aes(x=EV,y=BindP), colour="dodgerblue4")+
        theme_bw()+
        theme_new()+
        labs(title="The Mean Power Distribution of Senarario 2",y="Mean Power",x="Variance Explained (Percent)")
      p3 <- ggplot(data)+geom_line(aes(x=EV,y=ErelP), colour="chartreuse4")+
        theme_bw()+
        theme_new()+
        labs(title="The Mean Power Distribution of Senarario 3",y="Mean Power",x="Variance Explained (Percent)")
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
                                       "Senarario 1 Mean Power",
                                       "Senarario 2 Mean Power",
                                       "Senarario 3 Mean Power"
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
      }
      #cat(sdd,'\n')
      #ptm1 <- proc.time()
      #cat('CPU time Total:',ptm1-ptm,'\n')
      mmEinP = round(mean(PowerEindP),4)
      qunatEinP = round(quantile(PowerEindP),4)[2:4] #25%, median, 75%
      
      mmBinP = round(mean(PowerBindP),4)
      qunatBinP = round(quantile(PowerBindP),4)[2:4] #25%, median, 75%
      
      mmBrelP = round(mean(PowerRelBandP),4)
      
      qunatBrelP =  round(quantile(PowerRelBandP),4)[2:4]
      
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
        labs(title="The Histogram of Power of Senarario S1",x="Power",y="Proportion")
      p2 <- ggplot(data,aes(data$BindP))+
        geom_histogram(aes(x=data$BindP,y=(..count..)/sum(..count..)),
                       fill="dodgerblue4",
                       alpha=0.75
        )+
        theme_bw()+
        theme_new()+
        labs(title="The Histogram of Power of Senarario S2",x="Power",y="Proportion")
      p3 <- ggplot(data,aes(data$ErelP))+
        geom_histogram(aes(x=data$ErelP,y=(..count..)/sum(..count..)),
                       fill="chartreuse4",
                       alpha=0.75
        )+
        theme_bw()+
        theme_new()+
        labs(title="The Histogram of Power of Senarario S3",x="Power",y="Proportion")
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
      
      Power_Dist <- c("Mean Power","25% Quantile of Power","Medium of Power",
                      "75% Quantile of Power")
      combind.result <- data.frame(Power_Dist,Gene_Arc_1,Gene_Arc_2,Gene_Arc_3)
      colnames(combind.result) <- c("Power Distribution","Senarario S1","Senarario S2",
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
        n = Total # sample size
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
      }
      #cat(sdd,'\n')
      #ptm1 <- proc.time()
      #cat('CPU time Total:',ptm1-ptm,'\n')
      
      
      BrelMean <- apply(PowerRelBandP,1,mean)
      BinMean <- apply(PowerBindP,1,mean)
      EinMean <- apply(PowerEindP,1,mean)
      
      EV.Length <- length(EV)
      
      
      EV.Length <- length(EV)
      
      data <- data.frame(EV*100,EinMean,BinMean,BrelMean)
      colnames(data) <- c("EV","EindP","BindP","ErelP")
      #y.up <- min(1,(max(Power)+0.1))
      p1 <- ggplot(data)+geom_line(aes(x=EV,y=EindP), colour="#c0392b")+
        theme_bw()+
        theme_new()+
        labs(title="The Mean Power Distribution of Senarario 1",y="Mean Power",x="Variance Explained (Percent)",y="Proportion")
      p2 <- ggplot(data)+geom_line(aes(x=EV,y=BindP), colour="dodgerblue4")+
        theme_bw()+
        theme_new()+
        labs(title="The Mean Power Distribution of Senarario 2",y="Mean Power",x="Variance Explained (Percent)")
      p3 <- ggplot(data)+geom_line(aes(x=EV,y=ErelP), colour="chartreuse4")+
        theme_bw()+
        theme_new()+
        labs(title="The Mean Power Distribution of Senarario 3",y="Mean Power",x="Variance Explained (Percent)")
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
                                       "Senarario 1 Mean Power",
                                       "Senarario 2 Mean Power",
                                       "Senarario 3 Mean Power"
      )
      MeanPower.combine[,2:4] <- round(MeanPower.combine[,2:4],3)
      MeanPower.combine <- MeanPower.combine[c(1,3,5,7,9),]
      return(list(MeanPower.combine,p1=p1,p2=p2,p3=p3,p4=p4))
    }
    else if(TEST=='Burden Test'&length(EV)==1) {
      
      Power= NULL
      
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
          pEP = calcPowerL(J,SC,SP,E,n,alpha)
          
          Pow = c(Pow,pEP)
        }
        Power = rbind(Power,Pow)
      }
      
      PowerBindP <- PowerEindP <- PowerRelBandP <- Power
      
      mmEinP = round(mean(PowerEindP),4)
      qunatEinP = round(quantile(PowerEindP),4)[2:4] #25%, median, 75%
      
      mmBinP = round(mean(PowerBindP),4)
      qunatBinP = round(quantile(PowerBindP),4)[2:4] #25%, median, 75%
      
      mmBrelP = round(mean(PowerRelBandP),4)
      
      qunatBrelP =  round(quantile(PowerRelBandP),4)[2:4]
      
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
        labs(title="The Histogram of Power of Senarario S1",x="Power",y="Proportion")
      p2 <- ggplot(data,aes(data$BindP))+
        geom_histogram(aes(x=data$BindP,y=(..count..)/sum(..count..)),
                       fill="dodgerblue4",
                       alpha=0.75
        )+
        theme_bw()+
        theme_new()+
        labs(title="The Histogram of Power of Senarario S2",x="Power",y="Proportion")
      p3 <- ggplot(data,aes(data$ErelP))+
        geom_histogram(aes(x=data$ErelP,y=(..count..)/sum(..count..)),
                       fill="chartreuse4",
                       alpha=0.75
        )+
        theme_bw()+
        theme_new()+
        labs(title="The Histogram of Power of Senarario S3",x="Power",y="Proportion")
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
      
      Power_Dist <- c("Mean Power","25% Quantile of Power","Medium of Power",
                      "75% Quantile of Power")
      combind.result <- data.frame(Power_Dist,Gene_Arc_1,Gene_Arc_2,Gene_Arc_3)
      colnames(combind.result) <- c("Power Distribution","Senarario S1","Senarario S2",
                                    "Senarario S3")
      
      
      #return (list(combind.result,p1=NULL,p2=NULL,p3=NULL,p4=p4))
      return (list(combind.result,p1=p1,p2=p2,p3=p3,p4=p4))
    }else if(TEST=='Burden Test'&(length(EV)>1)){
      Power= NULL
      
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
          pEP = calcPowerL(J,SC,SP,E,n,alpha)
          
          Pow = c(Pow,pEP)
        }
        Power = rbind(Power,Pow)
        
      }
      
      PowerRelBandP <- PowerBindP <- PowerEindP <- Power
      
      
      
      BrelMean <- apply(PowerRelBandP,1,mean)
      BinMean <- apply(PowerBindP,1,mean)
      EinMean <- apply(PowerEindP,1,mean)
      
      EV.Length <- length(EV)
      
      
      EV.Length <- length(EV)
      
      data <- data.frame(EV*100,EinMean,BinMean,BrelMean)
      colnames(data) <- c("EV","EindP","BindP","ErelP")
      #y.up <- min(1,(max(Power)+0.1))
      p1 <- ggplot(data)+geom_line(aes(x=EV,y=EindP), colour="#c0392b")+
        theme_bw()+
        theme_new()+
        labs(title="The Mean Power Distribution of Senarario 1",y="Mean Power",x="Variance Explained (Percent)",y="Proportion")
      p2 <- ggplot(data)+geom_line(aes(x=EV,y=BindP), colour="dodgerblue4")+
        theme_bw()+
        theme_new()+
        labs(title="The Mean Power Distribution of Senarario 2",y="Mean Power",x="Variance Explained (Percent)")
      p3 <- ggplot(data)+geom_line(aes(x=EV,y=ErelP), colour="chartreuse4")+
        theme_bw()+
        theme_new()+
        labs(title="The Mean Power Distribution of Senarario 3",y="Mean Power",x="Variance Explained (Percent)")
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
                                       "Senarario 1 Mean Power",
                                       "Senarario 2 Mean Power",
                                       "Senarario 3 Mean Power"
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

