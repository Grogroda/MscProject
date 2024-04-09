import camb 
import numpy as np
from matplotlib import pyplot as plt

h2=0.6736**2
n=[-5,-4,-3,-2,-1,0,1,2,3,4,5]
Om=0.315
ombh2=0.02237
omch2=0.12
ob_list, oc_list=[ombh2+0.003*i for i in n],[Om-(ombh2+0.003*i) for i in n]

s8_list1=[]
pars=camb.CAMBparams()
plt.figure()
plt.title(r'Dependência de $\sigma_8$ com $\Omega_b$ para $\Omega_m=0.315$ fixo')
plt.xlabel(r'$\Omega_b=\Omega_m-\Omega_c$')
plt.ylabel(r'$\sigma_8$')
for omb in ob_list:
	pars.set_cosmology(H0=67.36, ombh2=omb, omch2=0.12, omk=0, tau=0.0544)
	pars.InitPower.set_params(As=(1e-10)*np.e**(3.044), ns=0.9649)
	pars.set_matter_power(redshifts=[0])
	
	results=camb.get_results(pars)
	powers=results.get_cmb_power_spectra(pars, CMB_unit='muK')
	sigma8=results.get_sigma8()[0]
	s8_list1.append(sigma8)
	
	print('Ob={0} feito!'.format(omb))
	
plt.scatter(ob_list,s8_list1)
plt.savefig('s8_omfixo.png')
