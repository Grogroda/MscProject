import numpy as np
from pyctg import ctg4py

OmegaM_fid=0.3
ls_data, ctg_data = ctg4py(OmegaM_fid)
ctg_sigmas=[i*0.1 for i in ctg_data]
print("Fake ctg calculated")

def Likelihood(OmegaM):

    print("Inside Likelihood")

    cl_theo=ctg4py(OmegaM)

    print("ctg_theo calculated")

    theory_array, data_array, sigma_array=np.array(cl_theo), np.array(ctg_data), np.array(ctg_sigmas)

    Xi2    = np.sum(((theory_array-data_array)/sigma_array)**2)
    sigSum = np.sum(np.log(sigma_array))
    logp = - sigSum - Xi2/2

    return logp

Omega_min, Omega_max=0.01, 0.99 #minimum value
n_vals=25

delta_Om=(Omega_max-Omega_min)/n_vals

Omegas, Likes=[],[]
OmegaM=Omega_min
for i in range(n_vals+1):
    print("i=", i)
    Like=Likelihood(OmegaM)
    
    Omegas.append(OmegaM)
    Likes.append(Like)

import matplotlib.pyplot as plt

plt.figure()
plt.plot(Omegas, Likes)
plt.xlabel(r"$\Omega_m$")
plt.ylabel(r'$\mathcal{L}(\Omega_m)$')
plt.savefig('Omegam_profile.png')

