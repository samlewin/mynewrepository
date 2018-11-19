from newtonraph import newton_raph
import math

def pstarexact(initvals, gamma=1.4):
    uL, uR, rhoL, rhoR, pL, pR = (initvals[0], initvals[1], initvals[2],
                                  initvals[3], initvals[4], initvals[5])
    TOL = 10**-6
    aL, aR = math.sqrt(gamma*pL/rhoL), math.sqrt(gamma*pR/rhoR)
    if not 2*aL/(gamma - 1) + 2*aR/(gamma - 1) > uR - uL: #check that the pressure positivity condition holds
        return 0, 0
    AL, AR = 2/((gamma + 1)*rhoL), 2/((gamma + 1)*rhoR)
    BL, BR = (gamma - 1)*pL/(gamma + 1), (gamma - 1)*pR/(gamma + 1)
    c, d = (gamma - 1)/(2*gamma), (gamma - 1)/(gamma + 1)
    pTR = ((aL + aR - (0.5*(gamma - 1)*(uR - uL)))/((aL/(pL**c)) + (aR/(pR**c))))**(1/c)

    def fL(p):
        if p > pL:
            return (p - pL)*math.sqrt(AL/(p + BL))
        else:
            return (2*aL/(gamma - 1))*((p/pL)**c - 1)
    def fR(p):
        if p > pR:
            return (p - pR)*math.sqrt(AR/(p + BR))
        else:
            return (2*aR/(gamma - 1))*((p/pR)**c - 1)
    def f(p):
        return fL(p)+fR(p) + uR - uL
    def fLderiv(p):
        """Derivative of fL above"""
        if p > pL:
            return math.sqrt(AL/(BL + p))*(1 - (p - pL)/(2*(BL + p)))
        else:
            return (1/(rhoL*aL))*(p/pL)**(-(gamma + 1)/(2*gamma))
    def fRderiv(p):
        """Derivative of fR above"""
        if p > pR:
            return math.sqrt(AR/(BR + p))*(1-(p - pR)/(2*(BR + p)))
        else:
            return (1/(rhoL*aR))*(p/pR)**(-(gamma + 1)/(2*gamma))
    def fderiv(p):
        return fLderiv(p)+fRderiv(p)

    try:
        pstar = newton_raph(f, fderiv, pTR, TOL)[0]
        leftshock = int(pstar > pL)
        rightshock = 2*int(pstar > pR)
    except:
        pstar, leftshock, rightshock = 0, 0, 0

    return pstar, leftshock + rightshock
