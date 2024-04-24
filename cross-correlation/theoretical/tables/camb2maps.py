import sys
from tqdm import tqdm
import camb
import healpy as hp
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
import matplotlib

matplotlib.rcParams.update({"font.size":15})

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
print("Ctt computed (size={}):".format(len(ctt)))
print(ctt)

OmegaM=0.3135
ls_ctg, ctg=pyctg.ctg4py(OmegaM)
print("Ctg computed")
ls_cgg, cgg=pyctg.cgg4py(OmegaM)
print("Cgg computed")

#Remember that, although we set lmax=128, camb still calculated some further multipoles, so we have to force the use of only 0 to 128

ls_tt, ctt=ls_ctt[:129], ctt[:129]
ls_gg, cgg=ls_cgg[:129], cgg[:129]
ls_tg, ctg=ls_ctg[:129], ctg[:129]

monodip=[0,0]
cgg=monodip+cgg
ctg=monodip+ctg

#I have to create a "final_table.dat" file with: ls | ctt | cgg | ctg

print('len(ls_tt)={0}\nctt={1}\ncgg={2}\nctg={3}'.format(len(ls_tt), ctt, cgg, ctg))
final=pd.DataFrame({'ls': ls_tt, 'ctt': ctt, 'cgg':cgg, 'ctg':ctg})
final.to_csv('final_table.dat', sep=' ', header=None)

# Test Zone:
temp_map=hp.sphtfunc.synfast(ctt, int((len(ctt)-1)/2), lmax=len(ctt)-1)
galaxy_map=hp.sphtfunc.synfast(cgg, int((len(cgg)-1)/2), lmax=len(cgg)-1)

#hp.write_map("CttMap.fits", temp_map, overwrite=True)
hp.mollview(temp_map, title='Temperature Map', remove_dip=True, unit=r'$\mu K$', cmap='RdYlBu_r')
hp.graticule()
plt.savefig('CMB_TempMap.png')

#hp.write_map("CggMap.fits", galaxy_map, overwrite=True)
hp.mollview(temp_map, title='Contraste de Galáxias', remove_dip=True, cmap='RdYlBu_r')
hp.graticule()
plt.savefig('Galaxy_Map.png')

Dtt=[ls_tt[i]*(ls_tt[i]+1)*ctt[i]/(2*np.pi) for i in range(len(ctt))]
ctt_anafast=hp.sphtfunc.anafast(temp_map, lmax=len(ctt)-1)
ls_anafast=np.arange(len(ctt_anafast))
Dtt_anafast=[ls_anafast[i]*(ls_anafast[i]+1)*ctt_anafast[i]/(2*np.pi) for i in range(len(ls_anafast))]

ls=ls_tt

plt.figure()
plt.plot(ls[2:], Dtt[2:], label="Theoretical")
plt.plot(ls_anafast[2:], Dtt_anafast[2:], label="Anafast")
plt.xlabel(r'$\ell$')
plt.ylabel(r'$C_{\ell}^{tt}$')
plt.xscale('log')
plt.savefig('Ctt_healpix.png')

ctg_anafast=hp.sphtfunc.anafast(temp_map, map2=galaxy_map, lmax=len(ctg)-1)
ls_anafast=np.arange(len(ctg_anafast))

plt.figure()
plt.plot(ls[2:], ctg[2:], label="Theoretical")
plt.plot(ls_anafast[2:], ctg_anafast[2:], label="Anafast")
plt.xlabel(r'$\ell$')
plt.ylabel(r'$C_{\ell}^{tg}$')
plt.xscale('log')
plt.savefig('Ctg_healpix.png')

#Actual code:

N=100

for i in tqdm(range(N), desc='Synthesizing maps'):
    temp_map=hp.sphtfunc.synfast(ctt, int((len(ctt)-1)/2), lmax=len(ctt)-1)
    galaxy_map=hp.sphtfunc.synfast(cgg, int((len(cgg)-1)/2), lmax=len(cgg)-1)

    almt=hp.sphtfunc.map2alm(temp_map, lmax=128)
    almg=hp.sphtfunc.map2alm(galaxy_map, lmax=128)

    hp.fitsfunc.write_alm('alms/cmb_alm{0}.fits'.format(i+1), almt, overwrite=True)
    hp.fitsfunc.write_alm('alms/galaxy_alm{0}.fits'.format(i+1), almg, overwrite=True)
