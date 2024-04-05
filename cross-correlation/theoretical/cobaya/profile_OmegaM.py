import numpy as np
from pyctg import ctg4py
import pandas as pd

'''
OmegaM_fid=0.3
ls_data, ctg_data = ctg4py(OmegaM_fid)
ctg_sigmas=[i*0.1 for i in ctg_data]
print("[profile_OmegaM.py] Fake ctg calculated")
'''
arr_names=['ls', 'Dtt', 'err_Dtt', 'Cgg', 'err_Cgg', 'Ctg', 'err_Ctg']
data_band1=pd.read_csv('../../data/Data_2016/errors_data_wmap9QVW_xsc1_jeffrey_dipfix_50000_ttggtg_new.dat', header=None, names=arr_names, sep=' ')

ls1=data_band1['ls']
Ctg1=data_band1["Ctg"]
print("ls=", ls1)
sigma_ctg1=data_band1['err_Ctg']

def Likelihood(OmegaM):

    print("Inside Likelihood")

    cl_theo=ctg4py(OmegaM)

    print("ctg_theo calculated")

    theory_array, data_array, sigma_array=np.array(cl_theo), np.array(Ctg1), np.array(sigma_ctg1)

    Xi2    = np.sum(((theory_array-data_array)/sigma_array)**2)
    sigSum = np.sum(np.log(sigma_array))
    logp = - sigSum - Xi2/2

    return logp

Omega_min, Omega_max=0.15, 0.6 
n_vals=100

delta_Om=(Omega_max-Omega_min)/n_vals

Omegas, Likes=[],[]
OmegaM=Omega_min

print("[profile_OmegaM] About to start profiling loop.")

for i in range(n_vals+1):
    print("i=", i)
    Like=Likelihood(OmegaM)
    
    Omegas.append(OmegaM)
    Likes.append(Like)

    OmegaM+=delta_Om

logpmax=max(Likes)
like_diff=[Likes[i]-logpmax for i in range(len(Likes))]
exp_like=[np.exp(logp) for logp in like_diff]

import matplotlib.pyplot as plt

print("Omegas=", Omegas)
print("Exp(logp-logpmax)=", exp_like)

plt.figure()
plt.plot(Omegas, exp_like)
plt.xlabel(r"$\Omega_m$")
plt.ylabel(r'$\mathcal{L}(\Omega_m)$')
plt.savefig('Omegam_profile.png')

