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
  geom_histogram(aes(x=data$EindP,y=..density..),
                 fill="#c0392b",
                 alpha=0.75
  )+
  theme_bw()+
  theme_new()+
  labs(title="The Histogram of Sample Size of Senarario 1",x="Sample Size")
p2 <- ggplot(data,aes(data$BindP))+
  geom_histogram(aes(x=data$BindP,y=..density..),
                 fill="dodgerblue4",
                 alpha=0.75
  )+
  theme_bw()+
  theme_new()+
  labs(title="The Histogram of Sample Size of Senarario 2",x="Sample Size")
p3 <- ggplot(data,aes(data$ErelP))+
  geom_histogram(aes(x=data$ErelP,y=..density..),
                 fill="chartreuse4",
                 alpha=0.75
  )+
  theme_bw()+
  theme_new()+
  labs(title="The Histogram of Sample Size of Senarario 3",x="Sample Size")
p4 <- ggplot(data,aes(data$JJ))+
  geom_histogram(aes(x=data$JJ,y=..density..),
                 fill="gray15",
                 alpha=0.75
  )+
  theme_bw()+
  theme_new()+
  labs(title="The Histogram of Parameters",x="Parameters")

Gene_Arc_1 <- c(mmEinP,qunatEinP)
Gene_Arc_2 <- c(mmBinP,qunatBinP)
Gene_Arc_3 <- c(mmBrelP,qunatBrelP)

Power_Dist <- c("Mean Sample Size","25% Quantile of Sample Size","Medium of Sample Size",
                "75% Quantile of Sample Size")
combind.result <- data.frame(Power_Dist,Gene_Arc_1,Gene_Arc_2,Gene_Arc_3)
colnames(combind.result) <- c("Sample Size Distribution","Senarario 1","Senarario 2",
                              "Senarario 3")


return (list(combind.result,p1=p1,p2=p2,p3=p3,p4=p4))


