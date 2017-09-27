source('PowerCalc_RareSample_A.R')

# Total = 10000
# ncase = 5000
# ncont = 5000
alpha= 0.0001
EV_low <- 0.01
EV_high <- 0.05
EV <- seq(EV_low,EV_high,(EV_high-EV_low)/10)
EV <- 0.5/100
PowerThr <- 0.8

get_Aprox_Sample(EV,PowerThr,alpha,PC=NA,TEST = 'SKAT',QT='CC',nameEsseble=NULL,JJ=NA,PRC=NA,ONESNP=F)
