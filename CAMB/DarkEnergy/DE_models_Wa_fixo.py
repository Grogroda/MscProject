from matplotlib import pyplot as plt
import camb
import numpy as np
import pandas as pd

pars_std=camb.CAMBparams()

#Usando as constantes definidas na dissertação do Esteban:

w0=-1

#Para cs2 constante:
wa_list=[0, 1/3, 2/3, 1]
cs2_list=[1e-2, 1e-1, 1]

c=1

data=pd.read_csv('COM_PowerSpect_CMB-TT-binned_R301.csv')

for wa in wa_list:

	plt.figure(c)
	plt.title(r'Espectro TT para $w(a)=w+(1-a)w_a$ com $wa={0:.2f}$'.format(wa))
	plt.figure(c+1)
	plt.title('Escala log x log')
	
	#Para cada wa, criarei um dataframe pandas no formato | l | ClTT001 | ClTT01 | ClTT1 | onde o numero depois do ClTT indica o valor de cs2
	#Esse dataframe vai ser armazenado na lista wa_sims
	
	sim=pd.DataFrame({'l':[], 'ClTT001':[], 'ClTT01':[], 'ClTT1':[]}, index=None)

	for cs2 in cs2_list:
		pars=camb.CAMBparams()
		pars.set_cosmology(H0=None, ombh2=0.02237, omch2=0.12, omk=0, tau=0.0544, cosmomc_theta=1.04092/100)
		pars.InitPower.set_params(As=(1e-10)*np.e**(3.044), ns=0.9649)
		pars.DarkEnergy.set_params(w=w0, wa=wa, cs2=cs2)
		pars.set_for_lmax(1700)
		pars.set_matter_power(redshifts=[0])

		results=camb.get_results(pars)
		powers=results.get_cmb_power_spectra(pars, CMB_unit='muK')
		sigma8=results.get_sigma8()[0]
		H0=results.hubble_parameter(0)

		Cl_Tot=powers['total']
		ls=np.arange(Cl_Tot.shape[0])
		ClTT=Cl_Tot[:,0]
		
		plt.figure(c)
		plt.plot(ls[2:], ClTT[2:], label=r'$c_s^2={0:.2f}$, $\sigma_8={1:.4f}$'.format(cs2, sigma8))
		plt.figure(c+1)
		plt.plot(ls[2:], ClTT[2:], label=r'$c_s^2={0:.2f}$'.format(cs2))
		
		sim['l']=ls
		
		if cs2==1e-2:
			sim['ClTT001']=ClTT
			
		if cs2==1e-1:
			sim['ClTT01']=ClTT
			
		if cs2==1:
			sim['ClTT1']=ClTT
	
	sim.to_csv(r'ClTT_for_Wa{0:.2f}.csv'.format(wa))
			
	plt.figure(c)
	plt.xlabel(r'$\ell$')
	plt.ylabel(r'$\ell(\ell+1)C_l/2\pi \; (\mu K)^2$')
	plt.xlim([0,max(ls)])
	plt.errorbar(data['l'], data['Dl'], fmt='.',label='Dados Planck', yerr=[data['-dDl'], data['+dDl']])
	plt.legend()
	#plt.show()
	plt.savefig('DE_wa_fixo{0:.2f}.png'.format(wa))
	
	plt.figure(c+1)
	plt.xlabel(r'$\ell$')
	plt.ylabel(r'$\ell(\ell+1)C_l/2\pi \; (\mu K)^2$')
	plt.errorbar(data['l'], data['Dl'], fmt='.',label='Dados Planck', yerr=[data['-dDl'], data['+dDl']])
	plt.xscale('log')
	plt.yscale('log')
	plt.xlim([1,max(ls)])
	plt.legend()
	#plt.show()
	plt.savefig('DE_wa_fixo{0:.2f}log.png'.format(wa))
	
	print('wa=', wa, ' terminado!')
	c+=2
