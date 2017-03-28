#
# Proposed Estimation
#

# Power calculation for SKAT under assumption: EVs are independent of MAFs
# Parameters:
# 	MAFCausal = MAF for J SNPs in a locus (tests)
# 	EVC = a coefficient of explained variation of a locus 
# 	level = level of significance
# 	SampleSize = sample size 
#   
#
# Output:
#  power = estimated power

calcPower_SC = function(MAFCausal,SC,EVC,SampleSize,aa,bb,level){
    n = SampleSize
	pj = MAFCausal
	EV = 0
    J = length(MAFCausal)
	lambda = n*dbeta(pj,aa,bb)^2*pj*(1-pj)

	c1 = sum(lambda^1)*(1+1/J*EV) 
	c2 = sum(lambda^2)*(1+2/J*EV) 
	c3 = sum(lambda^3)*(1+3/J*EV) 
	c4 = sum(lambda^4)*(1+4/J*EV)

	s1 = c3/c2^(3/2)
	s2 = c4/c2^2
    #save(c1,c2,c3,c4,pj,J,file='temp')
    if ((s1^2)>s2){
	a = 1/(s1-sqrt(s1^2-s2))
	delta = s1*a^3 - a^2
	l = a^2 - 2*delta
	#cat(l,'AL',EV,'\n')
	}else{
	l =c2^3/c3^2
	delta=0
	#cat(l,'A',EV,'\n')
	}

	ct = qchisq(1-level,df=l,nc=delta)

	mux = l+delta
	muq = c1
	sigmax = sqrt(2*(l+2*delta)) 
	sigmaq = sqrt(2*c2)

	ctt = sigmaq/sigmax*(ct - mux) + muq

	EV = EVC
	c1 = sum(lambda)*1+sum(lambda[1:SC])/SC*EV 
	c2 = sum(lambda^2)*1+2*sum(lambda[1:SC]^2)/SC*EV
	c3 =  sum(lambda^3)*1+3*sum(lambda[1:SC]^3)/SC*EV
	c4 =   sum(lambda^4)*1+4*sum(lambda[1:SC]^4)/SC*EV

	s1 = c3/c2^(3/2)
	s2 = c4/c2^2

	if ((s1^2)>s2){
	a = 1/(s1-sqrt(s1^2-s2))
	delta = s1*a^3 - a^2
	l = a^2 - 2*delta
	}else{
	l = c2^3/c3^2
	delta=0
	}

	mux = l+delta
	muq = c1
	sigmax = sqrt(2*(l+2*delta)) 
	sigmaq = sqrt(2*c2)
	ctt = sigmax/sigmaq*(ctt - muq) + mux
	powerCT = 1-pchisq(ctt,df=l,nc=delta)
	
return (powerCT)
}



# Power calculation for SKAT under assumption: betajs (genetic effects) are independent of MAFs
# Parameters:
# 	MAFCausal = MAF for J SNPs in a locus (tests)
# 	EVC = a coefficient of explained variation of a locus 
# 	level = level of significance
# 	SampleSize = sample size 
#   
#
# Output:
#  power = estimated power

calcPowerBPind_SC = function(MAFCausal,SC,EVC,SampleSize,aa,bb,level){
    n = SampleSize
	pj = MAFCausal
	EV = 0
    J = length(MAFCausal)
	lambda = n*dbeta(pj,aa,bb)^2*pj*(1-pj)
    Ep = sum(pj[1:SC])
	El1 = sum(lambda^1)
	El1p = sum((lambda^1*pj)[1:SC])
	El2 = sum(lambda^2)
	El2p = sum((lambda^2*pj)[1:SC])
	El3 = sum(lambda^3)
	El3p = sum((lambda^3*pj)[1:SC])
	El4 = sum(lambda^4)
	El4p = sum((lambda^4*pj)[1:SC])
	
	a1 = El1p/(Ep)
	a2 = El2p/(Ep)
	a3 = El3p/(Ep)
	a4 = El4p/(Ep)
	
	c1 = sum(lambda^1)*(1+1/J*EV) 
	c2 = sum(lambda^2)*(1+2/J*EV) 
	c3 = sum(lambda^3)*(1+3/J*EV) 
	c4 = sum(lambda^4)*(1+4/J*EV)

	s1 = c3/c2^(3/2)
	s2 = c4/c2^2
    #save(c1,c2,c3,c4,pj,J,file='temp')
    if ((s1^2)>s2){
	a = 1/(s1-sqrt(s1^2-s2))
	delta = s1*a^3 - a^2
	l = a^2 - 2*delta
	#cat(l,'AL',EV,'\n')
	}else{
	l =c2^3/c3^2
	delta=0
	#cat(l,'A',EV,'\n')
	}

	ct = qchisq(1-level,df=l,nc=delta)

	mux = l+delta
	muq = c1
	sigmax = sqrt(2*(l+2*delta)) 
	sigmaq = sqrt(2*c2)

	ctt = sigmaq/sigmax*(ct - mux) + muq

	EV = EVC
	c1 = sum(lambda)*1+1/1*EV*a1
	c2 = sum(lambda^2)*1+2/1*EV*a2
	c3 =  sum(lambda^3)*1+3/1*EV*a3
	c4 =   sum(lambda^4)*1+4/1*EV*a4

	s1 = c3/c2^(3/2)
	s2 = c4/c2^2

	if ((s1^2)>s2){
	a = 1/(s1-sqrt(s1^2-s2))
	delta = s1*a^3 - a^2
	l = a^2 - 2*delta
	}else{
	l = c2^3/c3^2
	delta=0
	}

	mux = l+delta
	muq = c1
	sigmax = sqrt(2*(l+2*delta)) 
	sigmaq = sqrt(2*c2)
	ctt = sigmax/sigmaq*(ctt - muq) + mux
	powerCT = 1-pchisq(ctt,df=l,nc=delta)
	
return (powerCT)
}


# Power calculation for SKAT under assumption: genetic effects are related to MAF through log10 function.
# Parameters:
# 	MAFCausal = MAF for J SNPs in a locus (tests)
# 	EVC = a coefficient of explained variation of a locus 
# 	level = level of significance
# 	SampleSize = sample size 
#   
#
# Output:
#  power = estimated power


calcPowerD_SC = function(MAFCausal,SC,EVC,SampleSize,gm,aa,bb,level){
n = SampleSize
	pj = MAFCausal
	EV = 0
    J = length(MAFCausal)
	lambda = n*dbeta(pj,aa,bb)^2*pj*(1-pj)
    Ep = sum(pj[1:SC])
	El1 = sum(lambda[1:SC]^1)
	El1p = sum((lambda^1*pj)[1:SC])
	El2 = sum(lambda[1:SC]^2)
	El2p = sum((lambda^2*pj)[1:SC])
	El3 = sum(lambda[1:SC]^3)
	El3p = sum((lambda^3*pj)[1:SC])
	El4 = sum(lambda[1:SC]^4)
	El4p = sum((lambda^4*pj)[1:SC])
	
	a1 = El1p - El1*Ep/SC
	a2 = El2p - El2*Ep/SC
	a3 = El3p - El3*Ep/SC
	a4 = El4p - El4*Ep/SC

	b1 = 1*El1/sum(lambda^1)
	b2 = 1*El2/sum(lambda^2)
	b3 = 1*El3/sum(lambda^3)
	b4 = 1*El4/sum(lambda^4)
	
	c1 = sum(lambda^1)*(1+1/J*EV) 
	c2 = sum(lambda^2)*(1+2/J*EV) 
	c3 = sum(lambda^3)*(1+3/J*EV) 
	c4 = sum(lambda^4)*(1+4/J*EV)

	s1 = c3/c2^(3/2)
	s2 = c4/c2^2
    #save(c1,c2,c3,c4,pj,J,file='temp')
    if ((s1^2)>s2){
	a = 1/(s1-sqrt(s1^2-s2))
	delta = s1*a^3 - a^2
	l = a^2 - 2*delta
	#cat(l,'AL',EV,'\n')
	}else{
	l =c2^3/c3^2
	delta=0
	#cat(l,'A',EV,'\n')
	}

	ct = qchisq(1-level,df=l,nc=delta)

	mux = l+delta
	muq = c1
	sigmax = sqrt(2*(l+2*delta)) 
	sigmaq = sqrt(2*c2)

	ctt = sigmaq/sigmax*(ct - mux) + muq

	EV = EVC
	c1 = sum(lambda)*(1+b1/SC*EV + 1*gm*a1/sum(lambda^1)) 
	c2 = sum(lambda^2)*(1+2*b2/SC*EV + 2*gm*a2/sum(lambda^2)) 
	c3 =  sum(lambda^3)*(1+3*b3/SC*EV + 3*gm*a3/sum(lambda^3))
	c4 =   sum(lambda^4)*(1+4*b4/SC*EV + 4*gm*a4/sum(lambda^4)) 

	s1 = c3/c2^(3/2)
	s2 = c4/c2^2

	if ((s1^2)>s2){
	a = 1/(s1-sqrt(s1^2-s2))
	delta = s1*a^3 - a^2
	l = a^2 - 2*delta
	}else{
	l = c2^3/c3^2
	delta=0
	}

	mux = l+delta
	muq = c1
	sigmax = sqrt(2*(l+2*delta)) 
	sigmaq = sqrt(2*c2)
	ctt = sigmax/sigmaq*(ctt - muq) + mux
	powerCT = 1-pchisq(ctt,df=l,nc=delta)
	
return (powerCT)
}

