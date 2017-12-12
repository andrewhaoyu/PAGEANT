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
  pj.vec <- rep(0,1000)
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
        n = Total # Power
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
          pj.vec[i] <- pj[1]
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
      
      
      data <- data.frame(JJ.vec = JJ.vec,
                         pj.vec = pj.vec)
      colnames(data) <- c("JJ","pj")
      p1 <- ggplot(data,aes(data$JJ))+
        geom_histogram(aes(x=data$JJ,y=(..count..)/sum(..count..)),
                       fill="dodgerblue4",
                       alpha=0.75
        )+
        theme_bw()+
        theme_new()+
        labs(title="The Histogram of Total number of variants (J)",x="Total number of variants (J)",y="Proportion")
      p2 <- ggplot(data,aes(pj))+
        geom_histogram(aes(x=data$pj,y=(..count..)/sum(..count..)),
                       fill="#c0392b",
                       alpha=0.75
        )+
        theme_bw()+
        theme_new()+
        labs(title="The Histogram of Minor Allele Frequency (MAF)",x="Minor Allele Frequency (MAF)",y="Proportion")
      
      p3 <- ggplot()+geom_blank()+
        theme_bw()+
        theme(
          panel.border = element_blank()
        )
      
      
      
      Gene_Arc_1 <- c(mmEinP,qunatEinP)
      Gene_Arc_2 <- c(mmBinP,qunatBinP)
      Gene_Arc_3 <- c(mmBrelP,qunatBrelP)
      
      Power_Dist <- c("Mean Power","25% Quantile of Power","Medium of Power",
                      "75% Quantile of Power")
      combind.result <- data.frame(Power_Dist,Gene_Arc_1,Gene_Arc_2,Gene_Arc_3)
      colnames(combind.result) <- c("Power Distribution","Scenario S1","Scenario S2",
                                    "Scenario S3")
      
      
      #return (list(combind.result,p1=NULL,p2=NULL,p3=NULL,p4=p4))
      return (list(combind.result,p1=p1,p2=p2,p3=p3))
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
        n = Total # Power
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
          pj.vec[i] <- pj[1]
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
      
      data <- data.frame(JJ=JJ.vec,pj.vec =pj.vec)
      colnames(data) <- c("JJ","pj")
      p1 <- ggplot(data,aes(data$JJ))+
        geom_histogram(aes(x=data$JJ,y=(..count..)/sum(..count..)),
                       fill="dodgerblue4",
                       alpha=0.75
        )+
        theme_bw()+
        theme_new()+
        labs(title="The Histogram of total number of variants (J)",x="total number of variants (J)",y="Proportion")
      p2 <- ggplot(data,aes(pj))+
        geom_histogram(aes(x=data$pj,y=(..count..)/sum(..count..)),
                       fill="#c0392b",
                       alpha=0.75
        )+
        theme_bw()+
        theme_new()+
        labs(title="The Histogram of Minor Allele Frequency (MAF)",x="Minor Allele Frequency (MAF)",y="Proportion")
      p3 <- ggplot()+geom_blank()+
        theme_bw()+
        theme(
          panel.border = element_blank()
        )
      
      
      
      
      MeanPower.combine <- data.frame(EV=EV*100,GeneI =EinMean,GeneII = BinMean,GeneIII=BrelMean)
      colnames(MeanPower.combine) <- c("EV(Percent)",
                                       "Scenario 1 Mean Power",
                                       "Scenario 2 Mean Power",
                                       "Scenario 3 Mean Power"
      )
      MeanPower.combine[,2:4] <- round(MeanPower.combine[,2:4],3)
      #MeanPower.combine <- MeanPower.combine[c(1,3,5,7,9,11),]
      return(list(MeanPower.combine,p1=p1,p2=p2,p3=p3))
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
        n = Total # Power
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
          pj.vec[i] <- pj[1]
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
      
      data <- data.frame(JJ.vec = JJ.vec,
                         pj.vec = pj.vec)
      colnames(data) <- c("JJ","pj")
      p1 <- ggplot(data,aes(data$JJ))+
        geom_histogram(aes(x=data$JJ,y=(..count..)/sum(..count..)),
                       fill="dodgerblue4",
                       alpha=0.75
        )+
        theme_bw()+
        theme_new()+
        labs(title="The Histogram of Total number of variants (J)",x="Total number of variants (J)",y="Proportion")
      p2 <- ggplot(data,aes(pj))+
        geom_histogram(aes(x=data$pj,y=(..count..)/sum(..count..)),
                       fill="#c0392b",
                       alpha=0.75
        )+
        theme_bw()+
        theme_new()+
        labs(title="The Histogram of Minor Allele Frequency (MAF)",x="Minor Allele Frequency (MAF)",y="Proportion")
      
      p3 <- ggplot()+geom_blank()+
        theme_bw()+
        theme(
          panel.border = element_blank()
        )
      
      
      
      Gene_Arc_1 <- c(mmEinP,qunatEinP)
      Gene_Arc_2 <- c(mmBinP,qunatBinP)
      Gene_Arc_3 <- c(mmBrelP,qunatBrelP)
      
      Power_Dist <- c("Mean Power","25% Quantile of Power","Medium of Power",
                      "75% Quantile of Power")
      combind.result <- data.frame(Power_Dist,Gene_Arc_1,Gene_Arc_2,Gene_Arc_3)
      colnames(combind.result) <- c("Power Distribution","Scenario S1","Scenario S2",
                                    "Scenario S3")
      
      
      #return (list(combind.result,p1=NULL,p2=NULL,p3=NULL,p4=p4))
      return (list(combind.result,p1=p1,p2=p2,p3=p3))
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
        n = Total # Power
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
          pj.vec[i] <- pj[1]
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
      
      data <- data.frame(JJ=JJ.vec,pj.vec =pj.vec)
      colnames(data) <- c("JJ","pj")
      p1 <- ggplot(data,aes(data$JJ))+
        geom_histogram(aes(x=data$JJ,y=(..count..)/sum(..count..)),
                       fill="dodgerblue4",
                       alpha=0.75
        )+
        theme_bw()+
        theme_new()+
        labs(title="The Histogram of total number of variants (J)",x="total number of variants (J)",y="Proportion")
      p2 <- ggplot(data,aes(pj))+
        geom_histogram(aes(x=data$pj,y=(..count..)/sum(..count..)),
                       fill="#c0392b",
                       alpha=0.75
        )+
        theme_bw()+
        theme_new()+
        labs(title="The Histogram of Minor Allele Frequency (MAF)",x="Minor Allele Frequency (MAF)",y="Proportion")
      p3 <- ggplot()+geom_blank()+
        theme_bw()+
        theme(
          panel.border = element_blank()
        )
      
      
      
      
      MeanPower.combine <- data.frame(EV=EV*100,GeneI =EinMean,GeneII = BinMean,GeneIII=BrelMean)
      colnames(MeanPower.combine) <- c("EV(Percent)",
                                       "Scenario 1 Mean Power",
                                       "Scenario 2 Mean Power",
                                       "Scenario 3 Mean Power"
      )
      MeanPower.combine[,2:4] <- round(MeanPower.combine[,2:4],3)
      #MeanPower.combine <- MeanPower.combine[c(1,3,5,7,9,11),]
      return(list(MeanPower.combine,p1=p1,p2=p2,p3=p3))
    }
    else if(TEST=='Burden Test'&length(EV)==1) {
      
      Power= NULL
      
      n = CASE*CONTROL/Total
      if (QT == 'QT'){
        n = Total # Power
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
      
      data <- data.frame(JJ.vec = JJ.vec,
                         pj.vec = pj.vec)
      colnames(data) <- c("JJ","pj")
      p1 <- ggplot(data,aes(data$JJ))+
        geom_histogram(aes(x=data$JJ,y=(..count..)/sum(..count..)),
                       fill="dodgerblue4",
                       alpha=0.75
        )+
        theme_bw()+
        theme_new()+
        labs(title="The Histogram of Total number of variants (J)",x="Total number of variants (J)",y="Proportion")
      p2 <- ggplot()+geom_blank()+
        theme_bw()+
        theme(
          panel.border = element_blank()
        )
      
      p3 <- ggplot()+geom_blank()+
        theme_bw()+
        theme(
          panel.border = element_blank()
        )
      
      
      
      Gene_Arc_1 <- c(mmEinP,qunatEinP)
      Gene_Arc_2 <- c(mmBinP,qunatBinP)
      Gene_Arc_3 <- c(mmBrelP,qunatBrelP)
      
      Power_Dist <- c("Mean Power","25% Quantile of Power","Medium of Power",
                      "75% Quantile of Power")
      combind.result <- data.frame(Power_Dist,Gene_Arc_1,Gene_Arc_2,Gene_Arc_3)
      colnames(combind.result) <- c("Power Distribution","Scenario S1","Scenario S2",
                                    "Scenario S3")
      
      
      #return (list(combind.result,p1=NULL,p2=NULL,p3=NULL,p4=p4))
      return (list(combind.result,p1=p1,p2=p2,p3=p3))
    }else if(TEST=='Burden Test'&(length(EV)>1)){
      Power= NULL
      
      n = CASE*CONTROL/Total
      if (QT == 'QT'){
        n = Total # Power
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
      
      data <- data.frame(JJ=JJ.vec,pj.vec =pj.vec)
      colnames(data) <- c("JJ","pj")
      p1 <- ggplot(data,aes(data$JJ))+
        geom_histogram(aes(x=data$JJ,y=(..count..)/sum(..count..)),
                       fill="dodgerblue4",
                       alpha=0.75
        )+
        theme_bw()+
        theme_new()+
        labs(title="The Histogram of total number of variants (J)",x="total number of variants (J)",y="Proportion")
      p2 <- ggplot()+geom_blank()+
        theme_bw()+
        theme(
          panel.border = element_blank()
        )
      p3 <- ggplot()+geom_blank()+
        theme_bw()+
        theme(
          panel.border = element_blank()
        )
      
      
      
      
      MeanPower.combine <- data.frame(EV=EV*100,GeneI =EinMean,GeneII = BinMean,GeneIII=BrelMean)
      colnames(MeanPower.combine) <- c("EV(Percent)",
                                       "Scenario 1 Mean Power",
                                       "Scenario 2 Mean Power",
                                       "Scenario 3 Mean Power"
      )
      MeanPower.combine[,2:4] <- round(MeanPower.combine[,2:4],3)
      #MeanPower.combine <- MeanPower.combine[c(1,3,5,7,9,11),]
      return(list(MeanPower.combine,p1=p1,p2=p2,p3=p3))
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

