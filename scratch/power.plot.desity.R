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
  geom_density(aes(x=data$EindP,y=..density..),
                 fill="#c0392b",
                 alpha=0.75
  )+
  theme_bw()+
  theme_new()+
  labs(title="The Histogram of Power of Senarario S1",x="Power",y="Density")
p2 <- ggplot(data,aes(data$BindP))+
  geom_density(aes(x=data$BindP,y=..density..),
                 fill="dodgerblue4",
                 alpha=0.75
  )+
  theme_bw()+
  theme_new()+
  labs(title="The Histogram of Power of Senarario S2",x="Power",y="Density")
p3 <- ggplot(data,aes(data$ErelP))+
  geom_density(aes(x=data$ErelP,y=..density..),
                 fill="chartreuse4",
                 alpha=0.75
  )+
  theme_bw()+
  theme_new()+
  labs(title="The Histogram of Power of Senarario S3",x="Power",y="Density")
p4 <- ggplot(data,aes(data$JJ))+
  geom_density(aes(x=data$JJ,y=..density..),
                 fill="gray15",
                 alpha=0.75
  )+
  theme_bw()+
  theme_new()+
  labs(title="The Histogram of Total number of variants (J)",x="Total number of variants (J)",y="Density")

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