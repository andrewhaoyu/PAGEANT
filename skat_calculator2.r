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

calcPower = function(MAFCausal,EVC,SampleSize,aa,bb,level){
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
	c1 = sum(lambda)*(1+1/J*EV) 
	c2 = sum(lambda^2)*(1+2/J*EV) 
	c3 =  sum(lambda^3)*(1+3/J*EV)
	c4 =   sum(lambda^4)*(1+4/J*EV) 

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

# Power calculation for SKAT
# Parameters:
# 	MAFCausal = MAFsof SNPs in a locus
# 	EVC =  vector of  coefficients of explained variations for all SNPs in a locus
# 	level = level of significance
# 	SampleSize = sample size 
#
# Output:
#  power = theoretical power


calcPowerTher = function(MAFCausal,EVC,SampleSize,aa,bb,level){
    n = SampleSize
	pj = MAFCausal
	EV = 0
    J = length(MAFCausal)
	lambda = n*dbeta(pj,aa,bb)^2*pj*(1-pj)*2

	c1 = sum(lambda^1*(1+1/J*EV)) 
	c2 = sum(lambda^2*(1+2/J*EV)) 
	c3 = sum(lambda^3*(1+3/J*EV)) 
	c4 = sum(lambda^4*(1+4/J*EV))

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

	EV = EVC*n
	c1 = sum(lambda^1*(1+1*EV)) 
	c2 = sum(lambda^2*(1+2*EV)) 
	c3 = sum(lambda^3*(1+3*EV)) 
	c4 = sum(lambda^4*(1+4*EV))

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

calcPowerBPind = function(MAFCausal,EVC,SampleSize,aa,bb,level){
    n = SampleSize
	pj = MAFCausal
	EV = 0
    J = length(MAFCausal)
	lambda = n*dbeta(pj,aa,bb)^2*pj*(1-pj)*2
    Ep = sum(pj)
	El1 = sum(lambda^1)
	El1p = sum(lambda^1*pj)
	El2 = sum(lambda^2)
	El2p = sum(lambda^2*pj)
	El3 = sum(lambda^3)
	El3p = sum(lambda^3*pj)
	El4 = sum(lambda^4)
	El4p = sum(lambda^4*pj)
	
	a1 = El1p/(El1*Ep)
	a2 = El2p/(El2*Ep)
	a3 = El3p/(El3*Ep)
	a4 = El4p/(El4*Ep)
	
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
	c1 = sum(lambda)*(1+1/1*EV*a1) 
	c2 = sum(lambda^2)*(1+2/1*EV*a2) 
	c3 =  sum(lambda^3)*(1+3/1*EV*a3)
	c4 =   sum(lambda^4)*(1+4/1*EV*a4) 

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


calcPowerD = function(MAFCausal,EVC,SampleSize,gm,aa,bb,level){
n = SampleSize
	pj = MAFCausal
	EV = 0
    J = length(MAFCausal)
	lambda = n*dbeta(pj,aa,bb)^2*pj*(1-pj)
    Ep = sum(pj)
	El1 = sum(lambda^1)
	El1p = sum(lambda^1*pj)
	El2 = sum(lambda^2)
	El2p = sum(lambda^2*pj)
	El3 = sum(lambda^3)
	El3p = sum(lambda^3*pj)
	El4 = sum(lambda^4)
	El4p = sum(lambda^4*pj)
	
	a1 = El1p - El1*Ep/J
	a2 = El2p - El2*Ep/J
	a3 = El3p - El3*Ep/J
	a4 = El4p - El4*Ep/J
	
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
	c1 = sum(lambda)*(1+1/J*EV + 1*gm*a1/sum(lambda^1)) 
	c2 = sum(lambda^2)*(1+2/J*EV + 2*gm*a2/sum(lambda^2)) 
	c3 =  sum(lambda^3)*(1+3/J*EV + 3*gm*a3/sum(lambda^3))
	c4 =   sum(lambda^4)*(1+4/J*EV + 4*gm*a4/sum(lambda^4)) 

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

