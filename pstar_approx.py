def pstarapprox(initvals, gamma=1.4):
    uL, uR, rhoL, rhoR, pL, pR = (initvals[0], initvals[1], initvals[2],
                                  initvals[3], initvals[4], initvals[5])
    aL, aR = (gamma*pL/rhoL)**0.5, (gamma*pR/rhoR)**0.5
    if not 2*aL/(gamma - 1) + 2*aR/(gamma - 1) > uR - uL :
        return 0
    rhobar, abar = 0.5*(rhoL + rhoR), 0.5*(aL + aR)
    Ppvrs = 0.5*(pL + pR) - 0.5*(uR - uL)*rhobar*abar

    return max(0, Ppvrs)
