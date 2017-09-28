source('PowerCalc_RareSample_A.R')

# Total = 10000
# ncase = 5000
# ncont = 5000
alpha= 0.0001
#EV_low <- 0.01
#EV_high <- 0.05
#EV <- seq(EV_low,EV_high,(EV_high-EV_low)/10)
EV <- 0.5/100
PowerThr <- 0.8

result <- 
  get_Aprox_Sample(EV,PowerThr,alpha,PRC=0.1,TEST = 'Burden Test',QT='CC',nameEsseble=NULL,JJ=NA,ONESNP=T,PC=0.1)
grid.arrange(result[[2]],result[[4]],result[[3]],result[[5]],ncol=2)
