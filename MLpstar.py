import numpy as np
import pandas as pd
from pstar_exact import pstarexact
import math, time

#Generate learning data
start = time.time()

np.random.seed(42)
velocitiesL = 10*(np.random.rand(10000, 1) - 0.5)
np.random.seed(43)
velocitiesR = 10*(np.random.rand(10000, 1) - 0.5)
np.random.seed(44)
densitiesL = 2*(np.random.rand(10000, 1) + 0.125)
np.random.seed(45)
densitiesR = 2*(np.random.rand(10000, 1) + 0.125)
np.random.seed(46)
pressuresL = 2*(np.random.rand(10000, 1))
np.random.seed(47)
pressuresR = 2*(np.random.rand(10000, 1))

initial_vals = np.concatenate((velocitiesL, velocitiesR, densitiesL, densitiesR, pressuresL, pressuresR), axis=1)
exactvals = np.array([pstarexact(initial_vals[i])[0] for i in range(10000)])
exactshock = np.array([pstarexact(initial_vals[i])[1] for i in range(10000)])
outputs = np.stack((exactvals, exactshock),axis=1)
inputs_outputs = np.concatenate((initial_vals, outputs), axis=1)
inputs_outputs_DF = pd.DataFrame(inputs_outputs, columns=["uL", "uR", "rhoL", "rhoR", "pL", "pR"] +
                                 ["pstarexact", "shocks_rarefactions"])
final_vals = inputs_outputs_DF.drop(inputs_outputs_DF[inputs_outputs_DF.pstarexact < 0.0000001].index)

inputs = final_vals[["uL","uR","rhoL", "rhoR", "pL", "pR"]].values

from sklearn.preprocessing import StandardScaler

scaler = StandardScaler()
scaled = scaler.fit_transform(inputs)

#Apply the Random Forest Regressor method
from sklearn.ensemble import RandomForestRegressor
forest_reg = RandomForestRegressor(n_estimators=30, max_features=6)
forest_reg.fit(scaled, final_vals["pstarexact"].values)

end = time.time() - start
print("Training took", end)
