from newtonraph import newton_raph
import math
import pylab as plt
import numpy as np
def solve(x,T,uL,uR,pL,pR,rhoL,rhoR,gamma,TOL):
    aL,aR=math.sqrt(gamma*pL/rhoL),math.sqrt(gamma*pR/rhoR)
    if not 2*aL/(gamma-1)+2*aR/(gamma-1)>uR-uL: #check that the pressure positivity condition holds
        raise ValueError('Vacuum created by initial data')
    AL,AR=2/((gamma+1)*rhoL),2/((gamma+1)*rhoR)
    BL,BR=(gamma-1)*pL/(gamma+1),(gamma-1)*pR/(gamma+1)
    c,d=(gamma-1)/(2*gamma),(gamma-1)/(gamma+1)
    pTR=((aL+aR-(0.5*(gamma-1)*(uR-uL)))/((aL/(pL**c))+(aR/(pR**c))))**(1/c)
    pPV=0.5*(pL+pR)-0.125*(uR-uL)*(rhoL+rhoR)*(aL+aR)
    def gL(p):
        return math.sqrt(AL/(p+BL))
    def gR(p):
        return math.sqrt(AR/(p+BR))
    p1=max(TOL,pPV)
    pTS=(gL(p1)*pL+gR(p1)*pR-(uR-uL))/(gL(p1)+gR(p1))
    p0,p1,p2,p3=pTR,max(TOL,pPV),max(TOL,pTS),(pL+pR)/2
    #Define the determining function f and its first derivative fderiv
    def fL(p):
        if p>pL:
            return (p-pL)*math.sqrt(AL/(p+BL))
        else:
            return (2*aL/(gamma-1))*((p/pL)**c-1)
    def fR(p):
        if p>pR:
            return (p-pR)*math.sqrt(AR/(p+BR))
        else:
            return (2*aR/(gamma-1))*((p/pR)**c-1)
    def f(p):
        return fL(p)+fR(p)+uR-uL
    def fLderiv(p):
        """Derivative of fL above"""
        if p>pL:
            return math.sqrt(AL/(BL+p))*(1-(p-pL)/(2*(BL+p)))
        else:
            return (1/(rhoL*aL))*(p/pL)**(-(gamma+1)/(2*gamma))
    def fRderiv(p):
        """Derivative of fR above"""
        if p>pR:
            return math.sqrt(AR/(BR+p))*(1-(p-pR)/(2*(BR+p)))
        else:
            return (1/(rhoL*aR))*(p/pR)**(-(gamma+1)/(2*gamma))
    def fderiv(p):
        return fLderiv(p)+fRderiv(p)
    pstar=newton_raph(f,fderiv,p0,TOL)[0]
    ustar=0.5*(uL+uR)+0.5*(fR(pstar)-fL(pstar))
    leftshock=pstar>pL
    rightshock=pstar>pR
    def density(x,t=T):
        S=x/t
        if S<ustar:
            if not leftshock: #left rarefaction
                rhostarL,astarL=rhoL*(pstar/pL)**(1/gamma),aL*(pstar/pL)**c
                SHL,STL=uL-aL,ustar-astarL
                def rdensityL(s): #function to compute the density within the rarefaction
                    return rhoL*(2/(gamma+1)+(gamma-1)*(uL-s/t)/((gamma+1)*aL))**(2/(gamma-1))
                if S<SHL:
                    return rhoL
                elif S>STL:
                    return rhostarL
                else:
                    return rdensityL(x)
            else: #left shock
                rhostarL=rhoL*(((pstar/pL)+d)/((d*pstar/pL)+1))
                SL=uL-aL*math.sqrt(((gamma+1)/(2*gamma)*(pstar/pL)+c))
                if S<SL:
                    return rhoL
                else:
                    return rhostarL
        else:
            if rightshock: #right shock
                rhostarR=rhoR*(((pstar/pR)+d)/(d*(pstar/pR)+1))
                SR=uR+aR*math.sqrt((gamma+1)/(2*gamma)*(pstar/pR)+c)
                if S<SR:
                    return rhostarR
                else:
                    return rhoR
            else: #right rarefaction
                rhostarR,astarR=rhoR*(pstar/pR)**(1/gamma),aR*(pstar/pR)**c
                SHR,STR=uR+aR,ustar+astarR
                def rdensityR(s): #function to compute the density within the rarefaction
                    return rhoR*(2/(gamma+1)-(gamma-1)*(uR-s/t)/((gamma+1)*aR))**(2/(gamma-1))
                if S<STR:
                    return rhostarR
                elif S>SHR:
                    return rhoR
                else:
                    return rdensityR(x)
    def velocity(x,t=T):
        S=x/t
        if S<ustar:
            if not leftshock: #left rarefaction
                rhostarL,astarL=rhoL*(pstar/pL)**(1/gamma),aL*(pstar/pL)**c
                SHL,STL=uL-aL,ustar-astarL
                def rvelocityL(s): #function to compute the velocity within the rarefaction
                    return (2/(gamma+1))*(aL+(gamma-1)/2*uL+x/t)
                if S<SHL:
                    return uL
                elif S>STL:
                    return ustar
                else:
                    return rvelocityL(x)
            else: #left shock
                rhostarL=rhoL*(((pstar/pL)+d)/((d*pstar/pL)+1))
                SL=uL-aL*math.sqrt(((gamma+1)/(2*gamma)*(pstar/pL)+c))
                if S<SL:
                    return uL
                else:
                    return ustar
        else:
            if rightshock: #right shock
                rhostarR=rhoR*(((pstar/pR)+d)/(d*(pstar/pR)+1))
                SR=uR+aR*math.sqrt((gamma+1)/(2*gamma)*(pstar/pR)+c)
                if S<SR:
                    return ustar
                else:
                    return uR
            else: #right rarefaction
                rhostarR,astarR=rhoR*(pstar/pR)**(1/gamma),aR*(pstar/pR)**c
                SHR,STR=uR+aR,ustar+astarR
                def rvelocityR(s): #function to compute the velocity within the rarefaction
                    return (2/(gamma+1))*(-aR+(gamma-1)/2*uR+x/t)
                if S<STR:
                    return ustar
                elif S>SHR:
                    return uR
                else:
                    return rvelocityR(x)
    def pressure(x,t=T):
        S=x/t
        if S<ustar:
            if not leftshock: #left rarefaction
                rhostarL,astarL=rhoL*(pstar/pL)**(1/gamma),aL*(pstar/pL)**c
                SHL,STL=uL-aL,ustar-astarL
                def rpressureL(s): #function to compute the pressure within the rarefaction
                    return pL*((2/(gamma+1))+d*(uL-x/t)/aL)**(1/c)
                if S<SHL:
                    return pL
                elif S>STL:
                    return pstar
                else:
                    return rpressureL(x)
            else: #left shock
                rhostarL=rhoL*(((pstar/pL)+d)/((d*pstar/pL)+1))
                SL=uL-aL*math.sqrt(((gamma+1)/(2*gamma)*(pstar/pL)+c))
                if S<SL:
                    return pL
                else:
                    return pstar
        else:
            if rightshock: #right shock
                rhostarR=rhoR*(((pstar/pR)+d)/(d*(pstar/pR)+1))
                SR=uR+aR*math.sqrt((gamma+1)/(2*gamma)*(pstar/pR)+c)
                if S<SR:
                    return pstar
                else:
                    return pR
            else: #right rarefaction
                rhostarR,astarR=rhoR*(pstar/pR)**(1/gamma),aR*(pstar/pR)**c
                SHR,STR=uR+aR,ustar+astarR
                def rpressureR(s): #function to compute the pressure within the rarefaction
                    return pR*((2/(gamma+1))-d*(uR-x/t)/aR)**(1/c)
                if S<STR:
                    return pstar
                elif S>SHR:
                    return pR
                else:
                    return rpressureR(x)
    def speeds():
        if leftshock:
            SL=uL-aL*math.sqrt(((gamma+1)/(2*gamma)*(pstar/pL)+c))
            S1=abs(SL)
        elif not leftshock:
            SHL=uL-aL
            astarL=aL*(pstar/pL)**c
            STL=ustar-astarL
            S1=max(abs(SHL),abs(STL))
        if rightshock:
            SR=uR+aR*math.sqrt((gamma+1)/(2*gamma)*(pstar/pR)+c)
            S2=abs(SR)
        elif not rightshock:
            astarR=aR*(pstar/pR)**c
            STR=ustar+astarR
            SHR=uR+aR
            S2=max(abs(STR),abs(SHR))
        return S1,S2
    return velocity(x),density(x),pressure(x),speeds()[0],speeds()[1]
##def ret(a,b,c,d,e,f,g,h,i,j):
##    return a,b,c,d,e,f,g,h,i,j
##vret=np.vectorize(ret)
vsolve=np.vectorize(solve,otypes=[object])
##u1=np.array([0.75,-2.0],dtype=object)
##u2=np.array([0.0,2.0],dtype=object)
##p1=np.array([1.0,0.4],dtype=object)
##p2=np.array([0.1,0.4],dtype=object)
##rho1=np.array([1.0,1.0],dtype=object)
##rho2=np.array([0.125,1.0],dtype=object)
##T=np.array([0.2,0.15],dtype=object)
##Gamma=np.array([1.4,1.4],dtype=object)
##tol=np.array([10**-6,10**-6],dtype=object)
##x=np.array([0,0],dtype=object)
##
##print(vsolve(x,T,u1,u2,p1,p2,rho1,rho2,Gamma,tol))
##xvals=np.linspace(0,1,200)
##yvals=yvals=[solve(xi-0.3,0.2,0.75,0.0,1.0,0.1,1.0,0.125,1.4,TOL=10**-6)[0] for xi in xvals]
##plt.plot(xvals,yvals)
##plt.show()
