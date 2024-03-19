import camb
import matplotlib
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

pars=camb.CAMBparams()
pars.set_cosmology(H0=None, ombh2=0.02237, omch2=0.12, omk=0, tau=0.0544, cosmomc_theta=1.04092/100)
pars.InitPower.set_params(As=(1e-10)*np.e**(3.044), ns=0.9649)
pars.set_for_lmax(int(input("\nlmax=")))
pars.set_matter_power(redshifts=[0.], kmax=2.0, nonlinear=True)

results=camb.get_results(pars)
power=results.get_cmb_power_spectra(pars, CMB_unit='muK')
kh, z, pk=results.get_matter_power_spectrum(minkh=1e-4, maxkh=2, npoints=200)

data_Cl=power['total']
data_cltt=data_Cl[:,0]
Dtt_full=[l*(l+1)*data_cltt[l]/(2*np.pi) for l in range(len(data_cltt))]
ls=np.arange(data_Cl.shape[0])

data_tt_ISW=pd.read_csv("ctt_ISW_LCDM.dat", sep=' ', header=None, names=["l", "ctt"])
ls_isw=data_tt_ISW['l']
ctt_isw=data_tt_ISW['ctt']
dtt_isw=[l*(l+1)*ctt_isw[l]/(2*np.pi) for l in range(len(ctt_isw))]

matplotlib.rcParams.update({'font.size':15})

plt.figure(figsize=(8,6))
plt.plot(ls[2:], data_cltt[2:], label="Full CMB")
plt.plot(ls_isw[2:], ctt_isw[2:], label="Late-ISW contribution")
plt.xscale('log')
plt.yscale('log')
plt.xlabel(r'$\ell$')
plt.ylabel(r'$\ell(\ell+1)C_\ell^{tt}/2\pi$')
plt.legend()
plt.savefig("Ctt_comparison.png")
