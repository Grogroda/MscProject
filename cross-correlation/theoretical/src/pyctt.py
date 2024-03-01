import numpy as np
import matplotlib.pyplot as plt
from ctypes import * 
from tqdm import tqdm

# Load the shared library
lib_correlations = CDLL("../src/libcorrelations.so")

# Define the function to calculate Ctg
ctt4py=lib_correlations.ctt

# Define the function signature (argument types and return type)
ctt4py.argtypes = [c_double, c_double, c_int, c_double, c_double]
ctt4py.restype = c_double

if __name__=='__main__':

    lmax = 96

    wM=0.31
    wL=1-wM
    h, bg=0.67, 1.

    ls,ctt=[],[]

    for l in tqdm(range(2,lmax+1)):
        ls.append(l)
        ctt_l=ctt4py(wL, wM, l, h, bg)
        ctt.append(ctt_l)

    print("ls={}".format(ls))
    print("ctt={}".format(ctt))

    Dtt=[ls[i]*(ls[i]+1)/(2*np.pi)*ctt[i] for i in range(len(ctt))]

    plt.figure()
    plt.plot(ls, Dtt)
    plt.xscale('log')
    plt.xlabel(r'$\ell$')
    plt.ylabel(r'$\ell(\ell+1)/2\pi C_\ell^{tt}$')
    plt.show()
    
