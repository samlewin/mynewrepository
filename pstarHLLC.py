def pstarHLLC(uL, uR, pL, pR, rhoL, rhoR, gamma=1.4):
    #Define initial variables, check pressure positivity condition
    aL, aR = (gamma*pL/rhoL)**0.5, (gamma*pR/rhoR)**0.5
    if not 2*aL/(gamma - 1) + 2*aR/(gamma - 1) > uR - uL :
        raise ValueError('Vacuum created by data')
    rhobar, abar = 0.5*(rhoL + rhoR), 0.5*(aL + aR)
    EL = rhoL*(0.5*uL**2 + pL/((gamma - 1)*rhoL))
    ER = rhoR*(0.5*uR**2 + pR/((gamma - 1)*rhoR))
    FL = [rhoL*uL, rhoL*uL**2 + pL, uL*(EL + pL)]
    FR = [rhoR*uR, rhoR*uR**2 + pR, uR*(ER + pR)]
    #Compure the pressure estimate
    Ppvrs = 0.5*(pL + pR) - 0.5*(uR - uL)*rhobar*abar
    pstar = max(0, Ppvrs)
    return pstar
