from ctypes import * 
import os
#import camb
from cobaya.model import get_model
import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
#from tqdm import tqdm

# Load the shared library
lib_correlations = CDLL("../src/libcorrelations.so")

# Define the function to calculate Ctg
ctg4py_raw=lib_correlations.ctg4py
cgg4py_raw=lib_correlations.cgg4py

# Define the function signature (argument types and return type)
ctg4py_raw.argtypes = [c_double, c_double, c_int, c_double, c_double, c_double, c_double, c_double, c_int, c_int, np.ctypeslib.ndpointer(dtype=np.float64, flags="C_CONTIGUOUS"), np.ctypeslib.ndpointer(dtype=np.float64, flags="C_CONTIGUOUS"), c_int]
ctg4py_raw.restype = c_double

cgg4py_raw.argtypes = [c_double, c_double, c_int, c_double, c_double, c_double, c_double, c_double, c_int, c_int, np.ctypeslib.ndpointer(dtype=np.float64, flags="C_CONTIGUOUS"), np.ctypeslib.ndpointer(dtype=np.float64, flags="C_CONTIGUOUS"), c_int]
cgg4py_raw.restype = c_double

lmax = 46 

As=1e-10*np.e**(3.044)
params = {'ombh2':0.02237, 'omch2':0.12, 'H0':67, 'omk':0., 'tau':0.0544,
        'As': As, 'ns':0.9649}

info = {
        'params':params,
        'likelihood':{'one':None},
        'theory':{'camb':None}
        }

'''
import camb

pars=camb.CAMBparams()
pars.set_cosmology(H0=67, ombh2=0.02237, omch2=0.12, omk=0, tau=0.0544)
pars.InitPower.set_params(As=(1e-10)*np.e**(3.044), ns=0.9649)
pars.set_for_lmax(lmax)
pars.set_matter_power(redshifts=[0.], kmax=2.0, nonlinear=True)

results=camb.get_results(pars)
khs, zs, pks=results.get_matter_power_spectrum(minkh=1e-4, maxkh=2, npoints=200)

matplotlib.rcParams.update({'font.size': 15})

plt.figure(figsize=(8,6))
plt.plot(khs, pks[0], label="Direct CAMB")
#plt.title("Matter Power Spectrum (Halofit)") 
plt.xlabel(r"$k/h [Mpc^{-1} h^{-1}]$")
plt.xscale("log")
plt.yscale("log")
plt.ylabel(r"$P(k/h) [h^{-3} Mpc^3]$")
#plt.savefig("Pk_directCAMB.png")
'''

def ctg4py(OmegaM):

    print("[pyctg.py] Inside ctg4py")
    
    #First calculate matter PS with Cobaya's CAMB wrapper:
    h = 0.67
    omb=0.02237/(h**2)
    omch2=(OmegaM-omb)*h**2

    params['omch2'] = omch2
    info['params']=params
    
    model = get_model(info)
    model.add_requirements({"Pk_grid":{"z":0, "k_max":2}})
    model.logposterior({})
    karr, zarr, pkarr = model.provider.get_Pk_grid()

    kh=karr/h
    pkh=pkarr[0]*h**3

    #pkarr is a 2D array representing P(k,z), each subarray is P(k) for a z requested in Pk_grid.

    '''
    #plt.figure()
    plt.plot(karr/h, pkarr[0]*h**3, label="Cobaya wrapper")
    #plt.title("P(k) with CAMB indirectly used by Cobaya")
    #plt.xlabel("k/h")
    #plt.ylabel("P(k/h)")
    #plt.xscale("log")
    #plt.yscale("log")
    plt.legend()
    plt.savefig("Pk_comparison.png")
    '''
   
    nks = karr.size 

    #Then calculate ctg:
    OmegaL = 1 - OmegaM
    z0 = 0.043
    beta = 1.825 
    lbda = 1.524
    bg = 1.37
    mode = 1
    ncalls = 200000
    #fname  = c_char_p(pkfname.encode("ascii"))
    ls=[]
    ctg = []
    for l in range(2, round(lmax)): #around 2-3 minutes for the whole spectrum
        print('ctg for l=', l)
        ls.append(l)
        print("About to start ctg")
        cl= ctg4py_raw(OmegaL, OmegaM, l, z0, beta, lbda, h, bg, mode, ncalls, kh, pkh, nks)
        print("ctg calculated")
        ctg.append(cl)

    return ls, ctg

def cgg4py(OmegaM):

    #First calculate matter PS with CAMB:
    h = 0.67
    omb=0.02237/(h**2)
    omch2=(OmegaM-omb)*h**2

    params['omch2'] = omch2
    info['params']=params
    
    model = get_model(info)
    model.add_requirements({"Pk_grid":{"z":0, "k_max":2.0}})
    model.logposterior({})
    karr, zarr, pkarr = model.provider.get_Pk_grid()
   
    kh=karr/h
    pkh=pkarr[0]*h**3
    #karr, pkarr = np.array(kh, dtype=np.float64), np.array(pk, dtype=np.float64)
    nks = karr.size 

    #Then calculate cgg:
    OmegaL = 1 - OmegaM
    z0 = 0.043
    beta = 1.825
    lbda = 1.524
    bg = 1.37
    mode = 1 #mode and ncalls don't change the results, bur are still needed
    ncalls = 200000
    ls=[]
    cgg = []
    for l in range(2, round(lmax)): #about 12 seconds per point, ~6mins for 54 points
        #print("cgg for l=", l)
        ls.append(l)
        cl= cgg4py_raw(OmegaL, OmegaM, l, z0, beta, lbda, h, bg, mode, ncalls, kh, pkh, nks)
        cgg.append(cl)

    return ls, cgg

###Testing

if __name__=='__main__':
    #import matplotlib.pyplot as plt

    OmegaM=(0.02237+0.12)/(0.67**2)
    print("OmegaM=", OmegaM)
    
    ls,ctg=ctg4py(OmegaM)
    print('ls=', ls)
    print('ctg({0})={1}'.format(OmegaM, ctg))

    ls,cgg=cgg4py(OmegaM)
    print('ls=', ls)
    print('cgg({0})={1}'.format(OmegaM, ctg))

    matplotlib.rcParams.update({'font.size': 15})

    plt.figure()
    plt.plot(ls,ctg)
    plt.xlabel(r'$\ell$')
    plt.ylabel(r'$C^{tg}$')
    plt.xscale('log')
    plt.savefig("pyctg_full_test.png")

    plt.figure()
    plt.plot(ls, cgg)
    plt.xlabel(r'$\ell$')
    plt.ylabel(r'$C^{gg}$')
    plt.xscale('log')
    plt.yscale('log')
    plt.savefig("pycgg_full_test.png")

    plt.figure(figsize=(14,6))
    #Double plot for dissertation:

    plt.subplot(121, box_aspect=0.75)
    plt.plot(ls, cgg)
    plt.xlabel(r'$\ell$')
    plt.ylabel(r'$C^{gg}$')
    plt.xscale('log')
    plt.yscale('log')

    plt.subplot(122, box_aspect=0.75)
    plt.plot(ls,ctg)
    plt.xlabel(r'$\ell$')
    plt.ylabel(r'$C^{tg} [\mu K]$')
    plt.xscale('log')

    plt.savefig("Correlations_DoublePlot.png")
