import numpy as np
from pyctg import ctg4py, cgg4py
import pandas as pd
import argparse

arr_names=['ls', 'Dtt', 'err_Dtt', 'Cgg', 'err_Cgg', 'Ctg', 'err_Ctg']
data_band1=pd.read_csv('../../data/Data_2016/errors_data_wmap9QVW_xsc1_jeffrey_dipfix_50000_ttggtg_new.dat', header=None, names=arr_names, sep=' ')
data_band2=pd.read_csv('../../data/Data_2016/errors_data_wmap9QVW_xsc2_jeffrey_dipfix_50000_ttggtg_new.dat', header=None, names=arr_names, sep=' ')
data_band3=pd.read_csv('../../data/Data_2016/errors_data_wmap9QVW_xsc3_jeffrey_dipfix_50000_ttggtg_new.dat', header=None, names=arr_names, sep=' ')
data_band4=pd.read_csv('../../data/Data_2016/errors_data_wmap9QVW_xsc4_jeffrey_dipfix_50000_ttggtg_new.dat', header=None, names=arr_names, sep=' ')

#Not ready to run with band='min'

bands=[data_band1, data_band2, data_band3, data_band4]

bands_dict={1:{},
            2:{},
            3:{},
            4:{}}

for band in range(1, 5):
    for col in arr_names:
        bands_dict[band][col]=bands[band-1][col]

def Likelihood_ctg(OmegaM, band=1, n=1, ncalls=1000000):

    print("Inside Likelihood_ctg")

    ls_ctg, ctg_theo=ctg4py(OmegaM, band=band, nmp=n, ncalls=ncalls)
    print("ctg_theo calculated")
    
    ctg_theory_array, ctg_data_array, ctg_sigma_array=np.array(ctg_theo), np.array(bands_dict[band]['Ctg']), np.array(bands_dict[band]['err_Ctg'])
    
    #logp calculation for ctg data:
    ctg_xi2    = np.sum(((ctg_theory_array-ctg_data_array)/ctg_sigma_array)**2)
    ctg_sigSum = np.sum(np.log(ctg_sigma_array))
    ctg_logp   = - ctg_sigSum - ctg_xi2/2

    return ctg_logp, ctg_xi2

def Likelihood_cgg(OmegaM, band=1, n=1, ncalls=1000000):

    print("Inside Likelihood_cgg")
    
    ls_cgg, cgg_theo=cgg4py(OmegaM, band=band, nmp=n, ncalls=ncalls)
    print("cgg_theo calculated")

    cgg_theory_array, cgg_data_array, cgg_sigma_array=np.array(cgg_theo), np.array(bands_dict[band]['Cgg']), np.array(bands_dict[band]['err_Cgg'])

    #logp calculation for cgg data:
    cgg_xi2    = np.sum(((cgg_theory_array-cgg_data_array)/cgg_sigma_array)**2)
    cgg_sigSum = np.sum(np.log(cgg_sigma_array))
    cgg_logp   = - cgg_sigSum - cgg_xi2/2

    return cgg_logp, cgg_xi2

def Likelihood_both(OmegaM, band=1, n=1, ncalls=1000000):

    cgg_logp, cgg_xi2=Likelihood_cgg(OmegaM, band=band, n=n, ncalls=ncalls)
    ctg_logp, ctg_xi2=Likelihood_ctg(OmegaM, band=band, n=n, ncalls=ncalls)

    return ctg_logp+cgg_logp, ctg_xi2+cgg_xi2

Omega_min, Omega_max=0.05, 0.95 
n_vals=100
band, n, ncalls=1,1,1000000

parser=argparse.ArgumentParser(description="Profile OmegaM Likelihood distribution using 2MASS+WMAP data.")
parser.add_argument('-b', '--band', type=int, help="Band of the 2MASS catalog (1,2,3,4) or 'min' to use minimizer parameters")
parser.add_argument('-n','--nprocess', type=int, help="Number of parallelized processes")
parser.add_argument('-N', '--ncalls', type=int, help="Number of Monte Carlo integration calls")
parser.add_argument('-c', '--correlation', type=str, help="Either 'ctg' or 'cgg' to calculate the likelihood for the corresponding correlation function. If 'both', returns the product of the likelihoods.") 
parser.add_argument('--npoints', type=int, help="Number of points (OmegaM's) to calculate")
parser.add_argument('--Omin', type=float, help="Minimum value of OmegaM to profile")
parser.add_argument('--Omax', type=float, help="Maximum value of OmegaM to profile")
parser.add_argument('-i', '--identifier', type=str, help="If you want to run multiple instances of the same profile separately, use different identifiers to make sure the output files are not being overwritten.")

args=parser.parse_args()

if args.band!=None:
    band=int(args.band)
print("band=", band)
if args.nprocess!=None:
    n=int(args.nprocess)
print("nmp=", n)
if args.ncalls!=None:
    ncalls=int(args.ncalls)
print('ncalls=', ncalls)
if args.npoints!=None:
    n_vals=int(args.npoints)
print('npoints=', n_vals)
if args.Omin!=None:
    Omega_min=float(args.Omin)
if args.Omax!=None:
    Omega_max=float(args.Omax)
print('Omega_min={0} & Omega_max={1}'.format(Omega_min, Omega_max))

delta_Om=(Omega_max-Omega_min)/n_vals

Omegas, Likes, Xi2s=[],[],[]
OmegaM=Omega_min

corr='ctg'

if args.correlation!=None:
    corr=args.correlation

if corr=='ctg':
    Likelihood=Likelihood_ctg
elif corr=='cgg':
    Likelihood=Likelihood_cgg
elif corr=='both':
    Likelihood=Likelihood_both
else:
    raise Exception("-c [--correlation] must be either None, 'ctg', 'cgg' or 'both'. Possibly running with 'ctg'.")

print("[profile_OmegaM] About to start profiling loop for {0}.".format(corr))

for i in range(n_vals): #Open in Omegamax, makes divided runs easier
    print("i=", i)
    Like, Xi2=Likelihood(OmegaM, band, n, ncalls)
    
    Omegas.append(OmegaM)
    Likes.append(Like)
    Xi2s.append(Xi2)

    OmegaM+=delta_Om

logpmax=max(Likes)
print('logpmax=', logpmax)
Xi2max=max(Xi2s)
print('Xi2max=', Xi2max)

like_diff=[Likes[i]-logpmax for i in range(len(Likes))]
Xi2_diff=[Xi2s[i]-Xi2max for i in range(len(Xi2s))]

exp_like=[np.exp(logp) for logp in like_diff]
exp_Xi2=[np.exp(xi2) for xi2 in Xi2_diff]

import matplotlib.pyplot as plt

print("Band=", band)
print("Omegas=", Omegas)
print("Exp(logp-logpmax)=", exp_like)
print("Exp(Xi2-Xi2max)=", exp_Xi2)

ident=1
if args.identifier!=None:
    ident=args.identifier

data=pd.DataFrame({"Omega":Omegas, "P":exp_like, "expXi2":exp_Xi2})
data.to_csv("like_profiles/ProfileData{0}_{1}_band{2}_Nmc{3:.0e}.dat".format(ident, corr, band, ncalls), sep=' ')

'''
plt.figure()
plt.plot(Omegas, exp_like)
plt.xlabel(r"$\Omega_m$")
plt.ylabel(r'$\mathcal{L}(\Omega_m)$')
plt.savefig('like_profiles/OmegamProfile_{0}_band{1}_Nmc{2:.0e}.png'.format(corr, band, ncalls))
'''
