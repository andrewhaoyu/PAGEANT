##
##
## Average over MAF
## Function of EV and J
##
##
source('skat_sample_calculatorsc.r')
source('sampler.r')

theme_new <- function(base_size = 12, base_family = "Helvetica"){
  theme_bw(base_size = base_size, base_family = base_family) %+replace%
    theme(plot.title=element_text(size=12,face="bold.italic",hjust = 0.5,vjust= 10),
          axis.title.x = element_text(size = 10,face = "bold",hjust = 0.5),
          axis.title.y =  element_text(size = 10,face = "bold",hjust = 0.5),
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
		save(PowerRelBandP,PowerBindP,PowerEindP,file='1')
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
      density.Ein <- density(PowerEindP)
      density.Bin <- density(PowerBindP)
      density.Brel <- density(PowerRelBandP)

      density.y.max <- max(c(quantile(density.Ein$y),quantile(density.Bin$y),quantile(density.Brel$y)))
      density.x.min <- min(c(density.Ein$x,density.Bin$x,density.Brel$x))
      density.x.max <- max(c(density.Ein$x,density.Bin$x,density.Brel$x))

      limits.x <- c(0,1)
      breaks.x <- seq(0,1,0.1)
      if((density.x.max-density.x.min)<0.1)
      {
        limits.x.low <- round(density.x.min,1)
        limits.x.high <- round(density.x.max,1)
        if(limits.x.low==limits.x.high){
          limits.x.low <- max((limits.x.low-0.01),0)
        }
        limits.x <- c(limits.x.low,limits.x.high)
        breaks.x <- seq(limits.x.low,limits.x.high,0.001)
      }

      Power <- c(PowerEindP,PowerBindP,PowerRelBandP)
      Method <- c(rep("Geneic Arc I",length(PowerEindP)),rep("Genetic Arc II",length(PowerBindP)),rep("Genetic Arc III",length(PowerRelBandP)))
      data <- data.frame(Power,Method)

      p <- ggplot(data,aes(x=Power))+
        geom_density(aes(group=Method,colour=Method,fill=Method),alpha=0.3,adjust=1)+
        scale_x_continuous(limits=limits.x,expand = c(0,0),breaks = breaks.x)+
        scale_y_continuous(expand = c(0,0),limits=c(0,(density.y.max+0.5)))+
        theme_new()+
        labs(title="The Density of Power",y="Density",x="Power")+
       +theme(title ="The Density of Power" )
      # p <- ggplot(data, aes(x=Power)) +
      #   geom_histogram(aes(y=..density..,fill=Method,colour=Method,group=Method),binwidth=.02, alpha=.3, position="identity")+
      #   scale_x_continuous(limits=c(0,1.02),expand = c(0,0),breaks = seq(0,1.02,0.1))+
      #   #scale_y_continuous(expand = c(0,0),limits=c(0,15))+
      #   theme_new()+
      #   labs(title="The Density of Power",y="Density")

      Gene_Arc_1 <- c(mmEinP)
      Gene_Arc_2 <- c(mmBinP)
      Gene_Arc_3 <- c(mmBrelP)
      Power_Dist <- c("Power")
      combind.result <- data.frame(Power_Dist,Gene_Arc_1,Gene_Arc_2,Gene_Arc_3)
      colnames(combind.result) <- c("Power Distribution","Genetic Arc I","Genetic Arc II",
                                    "Genetic Arc III")


      return (list(combind.result,p=NULL))
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
		save(PowerRelBandP,PowerBindP,PowerEindP,file='2')
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

          SC = ceiling(J*PC)
          SP = ceiling(PRC*SC)
          pEP = calcPowerSampleL(J,PC,PRC,E,PowerThr,alpha,Ratio)

          Pow = c(Pow,pEP)
        }
		PowerX = Pow[Pow!=Ratio*10^4]
		Pow[Pow==Ratio*10^4]=mean(PowerX)
        Power = rbind(Power,Pow)
      }

      PowerBindP <- PowerEindP <- PowerRelBandP <- Power
		save(PowerRelBandP,PowerBindP,PowerEindP,file='3')
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
      density.x.min <- min(c(density.Ein$x,density.Bin$x,density.Brel$x))
      density.x.max <- max(c(density.Ein$x,density.Bin$x,density.Brel$x))

      limits.x <- c(0,1)
      breaks.x <- seq(0,1,0.1)
      if((density.x.max-density.x.min)<0.1)
      {
        limits.x.low <- round(density.x.min,1)
        limits.x.high <- round(density.x.max,1)
        if(limits.x.low==limits.x.high){
          limits.x.low <- max((limits.x.low-0.01),0)
        }
        limits.x <- c(limits.x.low,limits.x.high)
        breaks.x <- seq(limits.x.low,limits.x.high,0.001)
      }

      Power <- c(PowerEindP,PowerBindP,PowerRelBandP)
      Method <- c(rep("Geneic Arc I",length(PowerEindP)),rep("Genetic Arc II",length(PowerBindP)),rep("Genetic Arc III",length(PowerRelBandP)))
      data <- data.frame(Power,Method)

      p <- ggplot(data,aes(x=Power))+
        geom_density(aes(group=Method,colour=Method,fill=Method),alpha=0.3,adjust=1)+
        scale_x_continuous(limits=limits.x,expand = c(0,0),breaks = breaks.x)+
        scale_y_continuous(expand = c(0,0),limits=c(0,(density.y.max+0.5)))+
        #scale_y_continuous(expand = c(0,0),limits=c(0,(density.y.max+0.5)))+
        theme_new()+
        labs(title="The Density of Power",y="Density",x="Power")
      # p <- ggplot(data, aes(x=Power)) +
      #   geom_histogram(aes(y=..density..,fill=Method,colour=Method,group=Method),binwidth=.02, alpha=.3, position="identity")+
      #   scale_x_continuous(limits=c(0,1.02),expand = c(0,0),breaks = seq(0,1.02,0.1))+
      #   #scale_y_continuous(expand = c(0,0),limits=c(0,15))+
      #   theme_new()+
      #   labs(title="The Density of Power",y="Density")

      Gene_Arc_1 <- c(mmEinP)
      Gene_Arc_2 <- c(mmBinP)
      Gene_Arc_3 <- c(mmBrelP)
      Power_Dist <- c("Power")
      combind.result <- data.frame(Power_Dist,Gene_Arc_1,Gene_Arc_2,Gene_Arc_3)
      colnames(combind.result) <- c("Power Distribution","Genetic Arc I","Genetic Arc II",
                                    "Genetic Arc III")
      return (list(combind.result,p=NULL))
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

          SC = ceiling(J*PC)
          SP = ceiling(PRC*SC)
          pEP = calcPowerLSample(J,SC,SP,E,PowerThr,alpha,Ratio)

          Pow = c(Pow,pEP)
        }
		PowerX = Pow[Pow!=Ratio*10^4]
		Pow[Pow==Ratio*10^4]=mean(PowerX)
        Power = rbind(Power,Pow)

      }

      PowerRelBandP <- PowerBindP <- PowerEindP <- Power
	  save(PowerRelBandP,PowerBindP,PowerEindP,file='2')
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
        #scale_x_continuous(limits=limits.x,expand = c(0,0),breaks = breaks.x)+
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
		save(PowerRelBandP,PowerBindP,PowerEindP,file='4')
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
      density.Ein <- density(PowerEindP)
      density.Bin <- density(PowerBindP)
      density.Brel <- density(PowerRelBandP)

      density.y.max <- max(c(quantile(density.Ein$y),quantile(density.Bin$y),quantile(density.Brel$y)))
      density.x.min <- min(c(density.Ein$x,density.Bin$x,density.Brel$x))
      density.x.max <- max(c(density.Ein$x,density.Bin$x,density.Brel$x))

      limits.x <- c(0,1)
      breaks.x <- seq(0,1,0.1)
      if((density.x.max-density.x.min)<0.1)
      {
        limits.x.low <- round(density.x.min,1)
        limits.x.high <- round(density.x.max,1)
        if(limits.x.low==limits.x.high){
          limits.x.low <- max((limits.x.low-0.01),0)
        }
        limits.x <- c(limits.x.low,limits.x.high)
        breaks.x <- seq(limits.x.low,limits.x.high,0.001)
      }

      Power <- c(PowerEindP,PowerBindP,PowerRelBandP)
      Method <- c(rep("Geneic Arc I",length(PowerEindP)),rep("Genetic Arc II",length(PowerBindP)),rep("Genetic Arc III",length(PowerRelBandP)))
      data <- data.frame(Power,Method)

      p <- ggplot(data,aes(x=Power))+ geom_density(aes(group=Method,colour=Method,fill=Method),alpha=0.3,adjust=1)+
        scale_x_continuous(limits=limits.x,expand = c(0,0),breaks = breaks.x)+
        scale_y_continuous(expand = c(0,0),limits=c(0,(density.y.max+0.5)))+
       
       theme_new()+
      labs(title="The Density of Power",y="Density",x="Power")
      # p <- ggplot(data, aes(x=Power)) +
      #   geom_histogram(aes(y=..density..,fill=Method,colour=Method,group=Method),binwidth=.02, alpha=.3, position="identity")+
      #   scale_x_continuous(limits=c(0,1.02),expand = c(0,0),breaks = seq(0,1.02,0.1))+
      #   #scale_y_continuous(expand = c(0,0),limits=c(0,15))+
      #   theme_new()+
      #   labs(title="The Density of Power",y="Density")

      Gene_Arc_1 <- c(mmEinP,qunatEinP)
      Gene_Arc_2 <- c(mmBinP,qunatBinP)
      Gene_Arc_3 <- c(mmBrelP,qunatBrelP)
      Power_Dist <- c("Mean Power","25% Quantile of Power","Medium of Power",
                      "75% Quantile of Power")
      combind.result <- data.frame(Power_Dist,Gene_Arc_1,Gene_Arc_2,Gene_Arc_3)
      colnames(combind.result) <- c("Power Distribution","Genetic Arc I","Genetic Arc II",
                                    "Genetic Arc III")


      return (list(combind.result,p=p))
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

      EV_x <- rep(EV*100,3)
      group <- c(rep("Genetic Arc I",EV.Length),rep("Genetic Arc II",EV.Length),
                 rep("Genetic Arc III",EV.Length))
      Power <- c(EinMean,BinMean,BrelMean)
      data <- data.frame(EV_x,group,Power)
      y.up <- min(1,(max(Power)+0.1))
      p <- ggplot(data)+geom_line(aes(x=EV_x,y=Power,colour=group))+
        #scale_x_continuous(limits=limits.x,expand = c(0,0),breaks = breaks.x)+
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
      density.x.min <- min(c(density.Ein$x,density.Bin$x,density.Brel$x))
      density.x.max <- max(c(density.Ein$x,density.Bin$x,density.Brel$x))

      limits.x <- c(0,1)
      breaks.x <- seq(0,1,0.1)
      if((density.x.max-density.x.min)<0.1)
      {
        limits.x.low <- round(density.x.min,1)
        limits.x.high <- round(density.x.max,1)
        if(limits.x.low==limits.x.high){
          limits.x.low <- max((limits.x.low-0.01),0)
        }
        limits.x <- c(limits.x.low,limits.x.high)
        breaks.x <- seq(limits.x.low,limits.x.high,0.001)
      }

      Power <- c(PowerEindP,PowerBindP,PowerRelBandP)
      Method <- c(rep("Geneic Arc I",length(PowerEindP)),rep("Genetic Arc II",length(PowerBindP)),rep("Genetic Arc III",length(PowerRelBandP)))
      data <- data.frame(Power,Method)

      p <- ggplot(data,aes(x=Power))+
        geom_density(aes(group=Method,colour=Method,fill=Method),alpha=0.3,adjust=1)+
        scale_x_continuous(limits=limits.x,expand = c(0,0),breaks = breaks.x)+
        scale_y_continuous(expand = c(0,0),limits=c(0,(density.y.max+0.5)))+
        #scale_y_continuous(expand = c(0,0),limits=c(0,(density.y.max+0.5)))+
        theme_new()+
        labs(title="The Density of Power",y="Density",x="Power")
      # p <- ggplot(data, aes(x=Power)) +
      #   geom_histogram(aes(y=..density..,fill=Method,colour=Method,group=Method),binwidth=.02, alpha=.3, position="identity")+
      #   scale_x_continuous(limits=c(0,1.02),expand = c(0,0),breaks = seq(0,1.02,0.1))+
      #   #scale_y_continuous(expand = c(0,0),limits=c(0,15))+
      #   theme_new()+
      #   labs(title="The Density of Power",y="Density")

      Gene_Arc_1 <- c(mmEinP,qunatEinP)
      Gene_Arc_2 <- c(mmBinP,qunatBinP)
      Gene_Arc_3 <- c(mmBrelP,qunatBrelP)
      Power_Dist <- c("Mean Power","25% Quantile of Power","Medium of Power",
                      "75% Quantile of Power")
      combind.result <- data.frame(Power_Dist,Gene_Arc_1,Gene_Arc_2,Gene_Arc_3)
      colnames(combind.result) <- c("Power Distribution","Genetic Arc I","Genetic Arc II",
                                    "Genetic Arc III")
      return (list(combind.result,p=p))
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

      EV_x <- rep(EV*100,3)
      group <- c(rep("Genetic Arc I",EV.Length),rep("Genetic Arc II",EV.Length),
                 rep("Genetic Arc III",EV.Length))
      Power <- c(EinMean,BinMean,BrelMean)
      data <- data.frame(EV_x,group,Power)
      y.up <- min(1,(max(Power)+0.1))
      p <- ggplot(data)+geom_line(aes(x=EV_x,y=Power,colour=group))+
        #scale_x_continuous(limits=limits.x,expand = c(0,0),breaks = breaks.x)+
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

