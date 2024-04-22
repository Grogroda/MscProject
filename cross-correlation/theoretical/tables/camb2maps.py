import sys
import camb
import healpy as hp
import numpy as np
import pandas as pd

sys.path.append('../cobaya/')
import pyctg 

pars=camb.CAMBparams()
pars.set_cosmology(H0=None, ombh2=0.02237, omch2=0.12, omk=0, tau=0.0544, cosmomc_theta=1.04092/100)
pars.InitPower.set_params(As=(1e-10)*np.e**(3.044), ns=0.9649)
pars.set_for_lmax(128)
pars.set_matter_power(redshifts=[0.], kmax=2.0, nonlinear=True)

results=camb.get_results(pars)
power=results.get_cmb_power_spectra(pars, CMB_unit='muK', raw_cl=True)
kh, z, pk=results.get_matter_power_spectrum(minkh=1e-4, maxkh=2, npoints=200)

data_Cl=power['total']
ctt=data_Cl[:,0]
ls_ctt=np.arange(data_Cl.shape[0])
print("Ctt computed")

OmegaM=0.3135
ls_ctg, ctg=pyctg.ctg4py(OmegaM)
print("Ctg computed")
ls_cgg, cgg=pyctg.cgg4py(OmegaM)
print("Cgg computed")

temp_map=hp.sphtfunc.synfast(ctt, int(len(ctt)-1)/2, lmax=len(ctt)-1)
hp.write_cl("CttMap.fits", temp_map, overwrite=True)
hp.mollview(temp_map, title='Temperature Map', remove_dip=True, norm='hist', unit=r'$\mu K$', cmap='RdYlBu_r')
hp.graticule()
plt.savefig('CMB_TempMap.png')

galaxy_map=hp.sphtfunc.synfast(cgg, int(len(cgg)-1)/2, lmax=len(cgg)-1)
hp.write_cl("CggMap.fits", galaxy_map, overwrite=True)
hp.mollview(temp_map, title='Contraste de Galáxias', remove_dip=True, norm='hist', cmap='RdYlBu_r')
hp.graticule()
plt.savefig('Galaxy_Map.png')

ctt_anafast=hp.sphtfunc.anafast(temp_map, lmax=len(ctt)-1)
ls_anafast=np.arange(len(ctt_anafast))

plt.figure()
plt.plot(ls_ctt[2:], ctt[2:], label="Theoretical")
plt.plot(ls_anafast[2:], ctt_anafast[2:], label="Anafast")
plt.xlabel(r'$\ell$')
plt.ylabel(r'$C_{\ell}^{tt}$')
plt.xscale('log')
plt.savefig('Ctt_healpix.png')

ctg_anafast=hp.sphtfunc.anafast(temp_map, map2=galaxy_map, lmax=len(ctg)-1)
ls_anafast=np.arange(len(ctg_anafast))

plt.figure()
plt.plot(ls_ctg[2:], ctg[2:], label="Theoretical")
plt.plot(ls_anafast[2:], ctg_anafast[2:], label="Anafast")
plt.xlabel(r'$\ell$')
plt.ylabel(r'$C_{\ell}^{tg}$')
plt.xscale('log')
plt.savefig('Ctg_healpix.png')
