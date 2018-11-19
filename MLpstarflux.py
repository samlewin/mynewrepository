import math, time
import pylab as plt
import numpy as np


def MLflux(uL, uR, pL, pR, rhoL, rhoR, gamma, pstar):
    #Define initial variables, check pressure positivity condition
    aL, aR = (gamma*pL/rhoL)**0.5, (gamma*pR/rhoR)**0.5
    if not 2*aL/(gamma - 1) + 2*aR/(gamma - 1) > uR - uL :
        raise ValueError('Vacuum created by data')
    rhobar, abar = 0.5*(rhoL + rhoR), 0.5*(aL + aR)
    EL = rhoL*(0.5*uL**2 + pL/((gamma - 1)*rhoL))
    ER = rhoR*(0.5*uR**2 + pR/((gamma - 1)*rhoR))
    FL = [rhoL*uL, rhoL*uL**2 + pL, uL*(EL + pL)]
    FR = [rhoR*uR, rhoR*uR**2 + pR, uR*(ER + pR)]

    #Calculate wave speed estimates
    def qL(p):
        if pstar <= pL:
            return 1
        else:
            return (1 + (gamma + 1)*(pstar/pL - 1)/(2*gamma))**0.5
    def qR(p):
        if pstar <= pR:
            return 1
        else:
            return (1 + (gamma + 1)*(pstar/pR - 1)/(2*gamma))**0.5
    SL, SR = uL - aL*qL(pstar), uR + aR*qR(pstar)
    #Basic outputs
    if 0 <= SL:
        return FL
    elif 0>= SR:
        return FR

    #Outputs with more calculation
    else:
        Sstar = (pR - pL + (rhoL*uL*(SL - uL)) - (rhoR*uR*(SR - uR)))/(rhoL*(SL - uL) - rhoR*(SR - uR))
        #Calculate the functions UK and UstarK
        UL = [rhoL, rhoL*uL, EL]
        UR = [rhoR, rhoR*uR, ER]
        UlistL = [1, Sstar, EL/rhoL + (Sstar - uL)*(Sstar + pL/(rhoL*(SL - uL)))]
        UlistR = [1, Sstar, ER/rhoR + (Sstar - uR)*(Sstar + pL/(rhoR*(SR - uR)))]
        UstarL = [rhoL*((SL - uL)/(SL - Sstar))*Ui for Ui in UlistL]
        UstarR = [rhoR*((SR - uR)/(SR - Sstar))*Ui for Ui in UlistR]

        #Calculate the HLLC flux
        FLstar = [FL[i] + SL*(UstarL[i] - UL[i]) for i in range(3)]
        FRstar = [FR[i] + SR*(UstarR[i] - UR[i]) for i in range(3)]

        # Designate the alternative outputs
        if SL < 0 and 0 <= Sstar:
            return FLstar
        else:
            return FRstar
