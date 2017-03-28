Total = 10000
ncase = 5000
ncont = 5000
alpha= 0.1

EV  =0.01 # 1%

result <- get_Aprox(EV,alpha,Total,CASE=ncase,CONTROL=ncont) #case control

PowerEindP <- result[[3]]
PowerBindP <- result[[4]]
PowerRelBandP <- result[[5]]