#
# Power calculations are for "most used" gene based test called SKAT (http://www.hsph.harvard.edu/skat/)
# It provides average power. We average out all unknown parameters related to gene-phenotype architecture.
#
# Power calculations can be done also for linear + various weights, C-lapha, Hotelling, etc.
# 
#


#
# Essential parameters: EV=% of variation in a trait explained by loci (gene).
# alpha = level of the test.  total size, number of cases, number of ncontrols.
# All you genes have ensembl names. So you can do power calculation for particular gene.
# Instead of name you can provide number of variants in a gene with MAF<0.01.
# Right now, the power calculation is based on data from Exome Aggregation Consortium, so it is exome based calculations. It may underestimate number of variants in a gene (Based on write-up).
# Power calculations can be done under two assumptions: 1) there is no relationship between MAF and % of variations explained by a variant; 2) there is no relationship between MAF and effect size (logOR). 3) we have also developed method where one specify relationship
# between MAF and genetic effect beta. 

# Under assumption #1
# No Gene Specification
source('PowerCalc_Rare.r')

Total = 10000
ncase = 5000
ncont = 5000
alpha= 0.0001

EV  =0.01 # 1%
EV <- seq(0.0005,0.006,(0.006-0.0005)/10)
PC <- 0.9
PRC <- 0

result <- get_AproxL(EV,alpha,Total,ncase,ncont,PC=PC,PRC = PRC) #case control

result <- get_AproxQ(EV,alpha,Total,QT='QT') #QT

result <- get_AproxQ(EV,alpha,Total,QT='QT',PC=0.1) # specify % of causal 


#
# Clapha
#
result <- get_AproxQ(EV,TEST='Calpha',alpha,Total,QT='QT',PC=0.1)
#
# Hotelling
#
result <- get_AproxQ(EV,TEST='Hotelling',alpha,Total,QT='QT',PC=0.1)


#
# For burden test
# Output for only one relationship
#
resulta <- get_AproxL(0.005,alpha,Total,ncase,ncont,QT='CC',PC=0.1)



mytable <- data.frame(
  Name = c("meanPower",
           "sd"),
  Value = c(
            round(get_Aprox(EV,alpha,Total,ncase,ncont)[1],
            get_Aprox(EV,alpha,Total,ncase,ncont)[2])
)
)

# Gene AKT1:  ENSG00000142208.


get_Aprox(EV,alpha,Total,ncase,ncont,"ENSG00000142208")
 
 
 #  You specify, gene size

 get_Aprox(EV,alpha,Total,ncase,ncont,JJ = 20)
 
 
# Under assumption #2 change get_Aprox to get_Aprox_BPind


