import camb 
import numpy as np
from matplotlib import pyplot as plt

h2=0.6736**2
n=[-5,-4,-3,-2,-1,0,1,2,3,4,5]
Om_list=[0.315+0.003*i for i in n]
ombh2, omch2=0.02237, 0.12
ob_list, oc_list=[Om*h2-omch2 for Om in Om_list],[Om*h2-ombh2 for Om in Om_list]

s8_list1=[]
pars1=camb.CAMBparams()
plt.figure(1)
plt.title(r'Dependência de $\sigma_8$ com $\Omega_m$, variando apenas $\Omega_b$')
plt.xlabel(r'$\Omega_m$')
plt.ylabel(r'$\sigma_8$')
for omb in ob_list:
	pars1.set_cosmology(H0=67.36, ombh2=omb, omch2=0.12, omk=0, tau=0.0544)
	pars1.InitPower.set_params(As=(1e-10)*np.e**(3.044), ns=0.9649)
	pars1.set_matter_power(redshifts=[0])
	
	results=camb.get_results(pars1)
	powers=results.get_cmb_power_spectra(pars1, CMB_unit='muK')
	sigma8=results.get_sigma8()[0]
	s8_list1.append(sigma8)
	
	print('Ob={0} feito!'.format(omb))
	
plt.scatter(Om_list,s8_list1)
plt.savefig('s8_obvar.png')

s8_list2=[]
pars2=camb.CAMBparams()
plt.figure(2)
plt.title(r'Dependência de $\sigma_8$ com $\Omega_m$, variando apenas $\Omega_c$')
plt.xlabel(r'$\Omega_m$')
plt.ylabel(r'$\sigma_8$')
for omc in oc_list:
	pars2.set_cosmology(H0=67.36, ombh2=0.02237, omch2=omc, omk=0, tau=0.0544)
	pars2.InitPower.set_params(As=(1e-10)*np.e**(3.044), ns=0.9649)
	pars2.set_matter_power(redshifts=[0])
	
	results=camb.get_results(pars2)
	powers=results.get_cmb_power_spectra(pars2, CMB_unit='muK')
	sigma8=results.get_sigma8()[0]
	s8_list2.append(sigma8)
	
	print('Oc={0} feito!'.format(omc))
	
plt.scatter(Om_list,s8_list2)
plt.savefig('s8_ocvar.png')
