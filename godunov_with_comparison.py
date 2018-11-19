import numpy as np, MLpstarflux
from ML_calc_pstar import ML_calc_pstar
from pstarHLLC import pstarHLLC
from pstarexact import pstarexact
from sklearn.metrics import mean_squared_error

def MLgodunov(rho, u, p, gamma, limit, dx, Ccfl = 0.9):
    """Implements the HLLC Godunov method with steps number of iterations,
    starting values rho, u, p (including the fictitious boundary data) and cell
    length dx. The number of cells is len(u)-2 and the CFL coefficient is Ccfl."""
    currentrho, currentu, currentp = rho, u, p
    currentU1 = currentrho
    currentU2 = currentrho*currentu
    currentU3 = currentrho*(0.5*currentu**2 + currentp/((gamma - 1)*currentrho))
    counter, time, l = 0, 0, len(currentu) - 1

    while time < limit:
        currenta = (gamma*currentp/currentrho)**0.5
        #Calculate the maximum speed and use it to determine the timestep
        speeds = np.vectorize(abs)(currentu)+np.vectorize(abs)(currenta)
        maxspeed = max(speeds)
        if counter <= 5:
            currentdt = 0.2*Ccfl*dx/maxspeed
        else:
            currentdt = Ccfl*dx/maxspeed
        time += currentdt
        counter += 1

        #Calculate the fluxes using the HLLC method
        pstarML = ML_calc_pstar(currentu[:-1], currentu[1:], currentrho[:-1], currentrho[1:], currentp[:-1], currentp[1:])
        pstarHLLC = pstarHLLC(currentu[:-1], currentu[1:], currentrho[:-1], currentrho[1:], currentp[:-1], currentp[1:])
        pstarexact = pstarexact(currentu[:-1], currentu[1:], currentrho[:-1], currentrho[1:], currentp[:-1], currentp[1:])

        rmse1 = mean_squared_error(pstarML, pstarexact)
        rmse2 = mean_squared_error(pstarHLLC, pstarexact)

        print("HLLC error:", rmse2, " ", "ML error:", rmse1)

        result = np.array([MLpstarflux.MLflux(currentu[j], currentu[j+1], currentp[j], currentp[j+1],
                               currentrho[j], currentrho[j+1], gamma, pstar[j]) for j in range(l)])
        fluxes1, fluxes2, fluxes3 = np.array(result[:,0]), np.array(result[:,1]), np.array(result[:,2])

        #Use the fluxes to calculate the new (capital) Us
        newU1 = currentU1[1:-1] + (currentdt/dx)*(fluxes1[:-1] - fluxes1[1:])
        newU2 = currentU2[1:-1] + (currentdt/dx)*(fluxes2[:-1] - fluxes2[1:])
        newU3 = currentU3[1:-1] + (currentdt/dx)*(fluxes3[:-1] - fluxes3[1:])
        #Invert these capital Us and update rho, u, p and
        currentrho = np.concatenate((np.array([newU1[0]]),newU1,np.array([newU1[-1]])))
        currentu = np.concatenate((np.array([newU2[0]/newU1[0]]), newU2/newU1, np.array([newU2[-1]/newU1[-1]])))
        currentp = np.concatenate((np.array([(newU3[0] - (0.5*newU2[0]**2)/newU1[0])*(gamma - 1)]),
                                  (newU3 - (0.5*newU2**2)/newU1)*(gamma - 1),
                                  np.array([(newU3[-1] - (0.5*newU2[-1]**2)/newU1[-1])*(gamma - 1)])))
        currentU1 = currentrho
        currentU2 = currentrho*currentu
        currentU3 = currentrho*(0.5*currentu**2 + currentp/((gamma - 1)*currentrho))
    return currentu[1:len(currentu)-1],currentrho[1:len(currentrho)-1],currentp[1:len(currentp)-1],counter
