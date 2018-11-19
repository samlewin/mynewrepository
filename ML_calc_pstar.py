import MLpstar as ML
import numpy as np

def ML_calc_pstar(uL, uR, rhoL, rhoR, pL, pR):
    inputs = np.stack((uL, uR, rhoL, rhoR, pL, pR), axis=1)
    scaled_inputs = ML.scaler.transform(inputs)
    pstar = ML.forest_reg.predict(scaled_inputs)
    return pstar
