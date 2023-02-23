import numpy as np
import SOAPify
import os
from numpy.linalg import norm
import h5py
from mpl_toolkits import mplot3d
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import sys
from scipy.optimize import curve_fit
from operator import itemgetter

nframes = 1  # nof consecutive frames

# 1. read SOAPified trajectory df from hdf5
os.chdir('Water')

# --------------------------------------------------------#
# with h5py.File('211_600.hdf5', 'r') as f:
#     mask = f["Trajectories/211_600/Trajectory"][0, :, 2] > 5.0
#
# with h5py.File('211_600soap.hdf5', 'r') as f:
#     X = f["SOAP/211_600"][:][:, mask]
# ---------------------------------------------------------#

# LIPIDS
# with h5py.File('martini22_293Ksoap.hdf5', 'r') as f:
#     X = f["/SOAP/martini22_293K"][5000:10001, :, :]

# ICE/WATER
with h5py.File('ice_watersoap.hdf5', 'r') as f:
    X = f["/SOAP/ice_water"][:, :, :]

X = SOAPify.fillSOAPVectorFromdscribe(
        X[:], l_max=8, n_max=8, atomTypes=["H", "O"], atomicSlices={'HH': slice(0, 324, None), 'HO': slice(324, 900, None), 'OH': slice(324, 900, None), 'OO': slice(900, 1224, None)})
#
# # WATER: atomTypes=["H", "O"], atomicSlices={'HH': slice(0, 324, None), 'HO': slice(324, 900, None), 'OH': slice(324, 900, None), 'OO': slice(900, 1224, None)
# # METALS: atomTypes=["Cu"], atomicSlices={'CuCu': slice(0, 324, None)}
## BTAw: atomTypes=["N"], atomicSlices={'NN': slice(0, 324, None)}
## ICO: atomTypes=["Au"], atomicSlices={'AuAu': slice(0, 324, None)}
# LIPIDS: atomTypes=["P"], atomicSlices={'PP': slice(0, 324, None)}
print(X.shape)
np.savez('ice_water.npz', name1=X)  # save dataset for COEXISTENCE


X = SOAPify.normalizeArray(X)
np.savez('X_normalized.npz', name1=X)  # save dataset for COEXISTENCE

# load dataset (if already present)
# X = np.load('X_normalized.npz')['name1']

# # 2. Get norm of SOAP vectors
# os.mkdir('norm_SOAP')
# os.chdir('norm_SOAP')
#
# norm_SOAP = np.zeros((X.shape[0], X.shape[1]))
# for frame in range(0, X.shape[0]):
#     for molecule in range(0, X.shape[1]):
#         SOAP_vector = X[frame, molecule, :]
#         l2norm = norm(SOAP_vector)
#         norm_SOAP[frame, molecule] = l2norm
#
# np.savez('norm_SOAP.npz', name1=norm_SOAP)
# os.chdir('../')
# os.chdir('1frame')


# 3. Get tSOAP (SOAP distance from frame t+1 and frame t-1)
# os.mkdir('dSOAP_1ns')
# os.chdir('dSOAP_1ns')
os.mkdir('dSOAP')
os.chdir('dSOAP')

# DERIVATA IN AVANTI
time_dSOAP = np.zeros((X.shape[0]-nframes, X.shape[1]))
for frame in range(nframes, X.shape[0]):
    for molecule in range(0, X.shape[1]):
        x = X[frame, molecule, :]
        y = X[frame-nframes, molecule, :]
        distance = SOAPify.simpleSOAPdistance(x, y)
        time_dSOAP[frame-nframes, molecule] = distance  # fill the matrix (each molecule for each frame)

np.savez('time_dSOAP.npz', name1=time_dSOAP)
# time_dSOAP = np.load('dSOAP/time_dSOAP.npz')['name1']
os.chdir('../')

# 4. time VARIATION tSOAP
os.mkdir('delta_dSOAP')
os.chdir('delta_dSOAP')

# delta_time_dSOAP = np.zeros((time_dSOAP.shape[0]-1, time_dSOAP.shape[1]))
# for frame in range(nframes, time_dSOAP.shape[0]-1):
#     for molecule in range(0, time_dSOAP.shape[1]):
#         delta = time_dSOAP[frame, molecule] - time_dSOAP[frame-nframes, molecule]
#         delta_time_dSOAP[frame-nframes, molecule] = delta

delta_time_dSOAP = []
for molecule in range(0, time_dSOAP.shape[1]):
    derivative = np.diff(time_dSOAP[:, molecule])
    delta_time_dSOAP.append(derivative)


np.savez('delta_time_dSOAP.npz', name1=delta_time_dSOAP)
# delta_time_dSOAP = np.load('delta_time_dSOAP.npz')['name1']
os.chdir('../')




