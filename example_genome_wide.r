# Input parameters
# K - number of genes
# m - number of discovries
# GEV - genome-wide EV
# grid - gread search (quick, adequate, complete)
# Total - total sample size
# ncase - numebr of cases (C-C study)
# ncont - numebr of controls (C-C study)
# alphaT  - level of the test, default value
# pc - proportion of causal (optopnal)
# type of test SKAT, Calpha, Hotelling
# Outcome: QT ,CC
# jj - numbrt variants per gene (optional)

K = 500
m = 1
GEV = 0.4
grid=20 # effects speed number of models  20 quick, 50 adequate, 100 complete  
epr = GEV/K
Total = 10000
ncase =5000
ncont = 5000
alphaT= 0.05/20000
pc = NA
test = 'SKAT'
qt = 'CC'
jj = NA

source('calc_whole.r')

Obj = calcGenomeLevel(K,m,grid,epr,alphaT,Total,CASE=ncase,CONTROL=ncont,PC=pc,TEST = test,QT=qt,nameEsseble=NULL,JJ=jj)


#
# Output to print: ModelS1_E, ModelS2_E, ModelS3_E - expected values minimum and maximum for three genetic models S1-E ind MAF, S2 - Beta ind MAF, S3 - log(MAF) 
#                  PrbMax1_m, PrbMax2_m, PrbMax3_m  - maximum probability of m-discoveries    
#                  PrbMin1_m, PrbMin2_m, PrbMin3_m  - minimum probability of m-discoveries  