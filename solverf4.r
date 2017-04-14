source('trgamma.r')
library('nleqslv')
solve = function(V){
#cat(epr,pr,'\n')
al = exp(V[1])
beta = exp(V[2])
if(al>1){
return (c(1000,1000))
}
e = al/beta
E = abs(mean(e-epr))/epr
P = abs(1-pgamma(0.00001,al,beta) - pr)/(pr)
if (is.na(P)){
return (c(1000,1000))
}
E[is.na(E)]=1000
P[is.na(P)]=1000
return (c(E,P))
}



get_solution = function(epr,pr,start,end){

epr<<-epr
pr<<-pr
xx=NULL
xj = NULL
xi = NULL
yy = NULL
x = seq(log(end),log(start),length.out=10)
y = log(exp(x)/(epr))
for (j in 1:10){
for (i in 1:10){
#cat(i,j,'\n')
res3 = try(a<-nleqslv(c(x[i],y[j]),solve,control=list(cndtol=10^(-100),btol=.001,ftol=10^(-25),xtol=10^(-25),maxit=1000))$x)
if(inherits(res3, "try-error")){
cat('Error Start','\n')
}else{
bb=solve(a)
yy = cbind(yy,a)
#cat(sum(abs(bb)),'\n')
xx = c(xx,sum(abs(bb)))
if(sum(abs(bb))< -0.0000001){break}
}
}
}

kkk = !is.na(xx)
xx = xx[kkk]
yy = yy[,kkk]

if (min(xx)>0.01){
return (list(error=min(xx),est=c(NA,NA)))
}

kx = which(xx==min(xx))
#cat(yy[,kx],'\n')
yy = yy[,kx[1]]
#cat(min(xx),'\n')
return (list(error=min(xx),est=exp(yy)))
}





