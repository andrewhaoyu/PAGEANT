Total = 10000
ncase = 5000
ncont = 5000
alpha= 0.1

EV  =0.01 # 1%
source('PowerCalc_Rare.r')
result <- get_Aprox(EV,alpha,Total,QT='QT') #case control
result <- get_Aprox2(EV,alpha,Total,CASE=ncase,CONTROL=ncont)

PowerEindP <- result[[3]]
PowerBindP <- result[[4]]
PowerRelBandP <- result[[5]]
qunatEinP = round(quantile(PowerEindP),4)[2:4] #25%, median, 75%
qunatBinP = round(quantile(PowerBindP),4)[2:4] #25%, median, 75%
qunatBrelP =  round(quantile(PowerRelBandP),4)[2:4]
density.y.max <- max(c(quantile(density.Ein$y,0.95),quantile(density.Bin$y,0.95),quantile(density.Brel$y,0.95)))
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
  theme_new()+
  labs(title="The Density of Power",y="Density")

p2 <- ggplot(data,aes(x=Power,fill=method))+
  geom_histogram(binwidth = .01,alpha=.5,position="identity")
p3 <- ggplot(data,aes(x=Power,fill=method))+geom_density()
  # scale_x_continuous(limits = c(0,1),expand = c(0,0))+
  # scale_y_continuous(limits = c(0,3.5),expand = c(0,0))+
  # theme_new()+
  # labs(title="The Density of Power",y="Density")

ggplot(dat, aes(x=rating, fill=cond)) +
  geom_histogram(binwidth=.5, alpha=.5, position="identity") +
  geom_vline(data=cdat, aes(xintercept=rating.mean,  colour=cond),
             linetype="dashed", size=1)








PowerBindP <- PowerEindP <- PowerRelBandP <- Power

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
