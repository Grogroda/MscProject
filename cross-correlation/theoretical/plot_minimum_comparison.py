import sys
from tqdm import tqdm
import numpy as np
from matplotlib import pyplot as plt
import matplotlib
import pandas as pd
import scipy.special as scp

matplotlib.rcParams.update({'font.size':18})

LCDM_correlations=pd.read_csv('./tables/final_table.dat', sep=' ', header=None, names=['ls', 'ctt', 'cgg', 'ctg'])
min_correlations=pd.read_csv('./cobaya/ctg_bandmin.dat', sep=' ', header=None, names=['ind', 'ls', 'ctg'])

ls=LCDM_correlations['ls']
ctg_2mass=LCDM_correlations['ctg']

pars_min={'lbda':4.9401, 'z0':0.1508, 'beta':3.088}
pars_2mass={'lbda':1.524, 'z0':0.043, 'beta':1.800}

def selection(z, z0, lbda, beta):
    k=(beta/scp.gamma(lbda))*(z/z0)**(beta*lbda-1)
    exp=np.exp(-(z/z0)**beta)/z0

    return k*exp

def plot_selections(zmin=0, zmax=0.5, npoints=10000, name='select_compare.png'):
    dz=(zmax-zmin)/npoints
    zs_2mass, selection_2mass=[], []
    zs_min, selection_min=[],[]

    z=zmin
    for i in tqdm(range(npoints), desc='Selection functions'):
        zs_2mass.append(z)
        zs_min.append(z)

        f_2mass=selection(z, pars_2mass['z0'], pars_2mass['lbda'], pars_2mass['beta'])
        f_min=selection(z, pars_min['z0'], pars_min['lbda'], pars_min['beta'])

        selection_2mass.append(f_2mass)
        selection_min.append(f_min)

        z+=dz

    plt.figure(figsize=(9,6))
    plt.plot(zs_2mass, selection_2mass, label='2MASS')
    plt.plot(zs_min, selection_min, label='Minimum')
    plt.xlabel(r'$z$')
    plt.ylabel(r'$\frac{dN}{dz}$', fontsize=24)
    plt.legend()
    plt.tight_layout()
    plt.savefig(name)

    return None

ls_min=min_correlations['ls']
ctg_min=min_correlations['ctg']

Dtg_2mass=[ls[i]*(ls[i]+1)*ctg_2mass[i]/(2*np.pi) for i in range(len(ls))]
Dtg_min=[ls_min[i]*(ls_min[i]+1)*ctg_min[i]/(2*np.pi) for i in range(len(ls_min))]

def plot_ctg():
    
    plt.figure(figsize=(9,6))
    plt.plot(ls[2:], Dtg_2mass[2:], label='2MASS')
    plt.plot(ls_min, Dtg_min, label='Minimum')
    plt.xlabel(r'$\ell$')
    plt.ylabel(r'$\ell(\ell+1)/2\pi$ $C^{tg}$')
    plt.xscale('log')
    plt.legend()
    plt.tight_layout()
    plt.savefig('minimum_Ctg_compare.png')

    return None

if __name__=='__main__':
    plot_selections()
    plot_ctg()
