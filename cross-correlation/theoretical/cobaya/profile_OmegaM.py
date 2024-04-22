import numpy as np
from pyctg import ctg4py, cgg4py
import pandas as pd
import argparse

'''
OmegaM_fid=0.3
ls_data, ctg_data = ctg4py(OmegaM_fid)
ctg_sigmas=[i*0.1 for i in ctg_data]
print("[profile_OmegaM.py] Fake ctg calculated")
'''
arr_names=['ls', 'Dtt', 'err_Dtt', 'Cgg', 'err_Cgg', 'Ctg', 'err_Ctg']
data_band1=pd.read_csv('../../data/Data_2016/errors_data_wmap9QVW_xsc1_jeffrey_dipfix_50000_ttggtg_new.dat', header=None, names=arr_names, sep=' ')
data_band2=pd.read_csv('../../data/Data_2016/errors_data_wmap9QVW_xsc2_jeffrey_dipfix_50000_ttggtg_new.dat', header=None, names=arr_names, sep=' ')
data_band3=pd.read_csv('../../data/Data_2016/errors_data_wmap9QVW_xsc3_jeffrey_dipfix_50000_ttggtg_new.dat', header=None, names=arr_names, sep=' ')
data_band4=pd.read_csv('../../data/Data_2016/errors_data_wmap9QVW_xsc4_jeffrey_dipfix_50000_ttggtg_new.dat', header=None, names=arr_names, sep=' ')

bands=[data_band1, data_band2, data_band3, data_band4]

bands_dict={1:{},
            2:{},
            3:{},
            4:{}}

for band in range(1, 5):
    for col in arr_names:
        bands_dict[band][col]=bands[band-1][col]

def Likelihood(OmegaM, band=1, n=1):

    print("Inside Likelihood")

    ls_ctg, ctg_theo=ctg4py(OmegaM, band, n)
    print("ctg_theo calculated")
    
    ls_cgg, cgg_theo=cgg4py(OmegaM, band, n)
    print("cgg_theo calculated")

    ctg_theory_array, ctg_data_array, ctg_sigma_array=np.array(ctg_theo), np.array(bands_dict[band]['Ctg']), np.array(bands_dict[band]['err_Ctg'])
    
    cgg_theory_array, cgg_data_array, cgg_sigma_array=np.array(cgg_theo), np.array(bands_dict[band]['Cgg']), np.array(bands_dict[band]['err_Cgg'])

    #logp calculation for ctg data:
    ctg_xi2    = np.sum(((ctg_theory_array-ctg_data_array)/ctg_sigma_array)**2)
    ctg_sigSum = np.sum(np.log(ctg_sigma_array))
    ctg_logp   = - ctg_sigSum - ctg_xi2/2

    #logp calculation for cgg data:
    cgg_xi2    = np.sum(((cgg_theory_array-cgg_data_array)/cgg_sigma_array)**2)
    cgg_sigSum = np.sum(np.log(cgg_sigma_array))
    cgg_logp   = - cgg_sigSum - cgg_xi2/2

    #joint logp:
    logp=ctg_logp+cgg_logp

    return logp

Omega_min, Omega_max=0.15, 0.6 
n_vals=100

delta_Om=(Omega_max-Omega_min)/n_vals

Omegas, Likes=[],[]
OmegaM=Omega_min

band, n=1,1
parser=argparse.ArgumentParser()
parser.add_argument('-b', '--band')
parser.add_argument('-n','--nprocess')
args=parser.parse_args()

if args.band!=None:
    band=int(args.band)
print("band=", band)
if args.nprocess!=None:
    n=int(args.nprocess)
print("n=", n)

print("[profile_OmegaM] About to start profiling loop.")

for i in range(n_vals+1):
    print("i=", i)
    Like=Likelihood(OmegaM, band, n)
    
    Omegas.append(OmegaM)
    Likes.append(Like)

    OmegaM+=delta_Om

logpmax=max(Likes)
like_diff=[Likes[i]-logpmax for i in range(len(Likes))]
exp_like=[np.exp(logp) for logp in like_diff]

import matplotlib.pyplot as plt

print("Band=", band)
print("Omegas=", Omegas)
print("Exp(logp-logpmax)=", exp_like)

data=pd.DataFrame({"Omega":Omegas, "Like":exp_like})
data.to_csv("like_profiles/ProfileData_band{0}_Nmc{1}.dat".format(band, '2e5'), sep=' ')

plt.figure()
plt.plot(Omegas, exp_like)
plt.xlabel(r"$\Omega_m$")
plt.ylabel(r'$\mathcal{L}(\Omega_m)$')
plt.savefig('like_profiles/OmegamProfile_band{0}_Nmc{1}.png'.format(band, '2e5'))

