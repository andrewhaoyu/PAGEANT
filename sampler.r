#
# Sampler
# 

dist = function(T){
V = as.vector(T)
L = length(V)
Total = sum(V)
P = V/Total
S = 0
D=NULL
for (i in 1:L){
S = S+P[i]
D = c(D,S)
}
return (D)
}

sampleD=function(dist,name,J){
p = sort(runif(J))
nn= NULL
kkk = 0
for (i in 1:J){
kkk0 = sum(dist<p[i])
kkk = kkk + kkk0
nn = c(nn,name[kkk+1])
dist = dist[-(1:kkk0)]
#cat(length(dist),'\n')
}
return (nn)
}