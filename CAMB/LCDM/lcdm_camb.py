import camb
import numpy as np
import pandas as pd

approx=input("Linear (L) or halofit (H)? ")
if approx=='H':
    nlinear=True
else:
    nlinear=False

pars=camb.CAMBparams()
pars.set_cosmology(H0=None, ombh2=0.02237, omch2=0.12, omk=0, tau=0.0544, cosmomc_theta=1.04092/100)
pars.InitPower.set_params(As=(1e-10)*np.e**(3.044), ns=0.9649)
pars.set_for_lmax(int(input("\nlmax=")))
pars.set_matter_power(redshifts=[0.], kmax=2.0, nonlinear=nlinear)

results=camb.get_results(pars)
power=results.get_cmb_power_spectra(pars, CMB_unit='muK')
kh, z, pk=results.get_matter_power_spectrum(minkh=1e-4, maxkh=2, npoints=200)

data_Cl=power['total']
data_cltt=data_Cl[:,0]
ls=np.arange(data_Cl.shape[0])
dataF_ctt=pd.DataFrame({"l":ls, "Cltt":data_cltt})
dataF_ctt.to_csv('Cltt_camb.dat', sep=' ', header=False, index=False)

data_matter=pd.DataFrame({'k':kh, 'P(k)':pk[0]}, index=None)
data_matter.to_csv('matterPS.dat', sep=' ', header=False, index=False)

# add plot option later
