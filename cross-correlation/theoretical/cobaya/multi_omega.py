from pyctg import ctg4py, cgg4py
import pandas as pd
import numpy as np
import camb
from tqdm import tqdm
from matplotlib import pyplot as plt

omegas=[0.2,0.3,0.45,0.6,0.7]

plots={} #{omegaMi:{'ctg':[ls, ctgs], 'cgg':[ls, cggs]}} where ls and ctgs are lists containing multipolesl and ctg(l) to plot

fig_ctg, ax_ctg=plt.subplots()
ax_ctg.set_xlabel(r'$\ell$')
ax_ctg.set_ylabel(r'$C_{tg}$')
ax_ctg.set_xscale('log')

fig_cgg, ax_cgg=plt.subplots()
ax_cgg.set_xlabel(r'$\ell$')
ax_cgg.set_ylabel(r'$C_{gg}$')
ax_cgg.set_xscale('log')
ax_cgg.set_yscale('log')

for i in tqdm(range(len(omegas))):
    
    #First calculate matter PS with CAMB:

    H0=67
    h=H0/100
    omb=0.02237/(h**2)
    omch2=(omegas[i]-omb)*h**2 #=0.12 for LCDM
    print("OmegaM={0} & OmegaC*h^2={1}".format(omegas[i], omch2))
    
    pars=camb.CAMBparams()
    pars.set_cosmology(H0=67, ombh2=0.02237, omch2=omch2, omk=0, tau=0.0544)
    pars.InitPower.set_params(As=(1e-10)*np.e**(3.044), ns=0.9649)
    pars.set_for_lmax(128)
    pars.set_matter_power(redshifts=[0.], kmax=2.0, nonlinear=True)

    results=camb.get_results(pars)
    power=results.get_cmb_power_spectra(pars, CMB_unit='muK')
    kh, z, pk=results.get_matter_power_spectrum(minkh=1e-4, maxkh=2, npoints=200)

    data_Cl=power['total']
    data_cltt=data_Cl[:,0]
    ls=np.arange(data_Cl.shape[0])
    dataF_ctt=pd.DataFrame({"l":ls, "Cltt":data_cltt})
    dataF_ctt.to_csv('Cltt_Om={}_multiomegas.dat'.format(omegas[i]), sep=' ', header=False, index=False)

    data_matter=pd.DataFrame({'k':kh, 'P(k)':pk[0]}, index=None)
    data_matter.to_csv('../tables/pk_3dmatter.dat', sep=' ', header=False, index=False)

    #Then calculate ctg and cgg using that spectrum

    omegaM=omegas[i]
    plots[omegaM]={'ctg':[], 'cgg':[]}
    
    ls, ctg=ctg4py(omegaM)
    ls, cgg=cgg4py(omegaM)

    ax_ctg.plot(ls, ctg, label=r'$\Omega_M={}$'.format(omegaM))
    ax_cgg.plot(ls, cgg, label=r'$\Omega_M={}$'.format(omegaM))
    
    plots[omegaM]['ctg'].append(ls)
    plots[omegaM]['ctg'].append(ctg)
    plots[omegaM]['cgg'].append(ls)
    plots[omegaM]['cgg'].append(cgg)

ax_ctg.legend()
ax_cgg.legend()
fig_ctg.savefig('MultiOmega_ctg.png')
fig_cgg.savefig("MultiOmega_cgg.png")

#print(plots)
