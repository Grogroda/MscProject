from matplotlib import pyplot as plt
import matplotlib
import camb
import numpy as np
import pandas as pd

matplotlib.rcParams.update({'font.size':15})
pars_std=camb.CAMBparams()

#Usando as constantes definidas na dissertação do Esteban:

w0=-1

#Para cs2 constante:
wa_list=[0, 1/3, 2/3, 1]
cs2_list=[1e-2, 1e-1, 1]

c=1

data=pd.read_csv('COM_PowerSpect_CMB-TT-binned_R301.csv')

for cs2 in cs2_list:


	fig1=plt.figure(c)
	#plt.title(r'Espectro TT para $w(a)=w+(1-a)w_a$ com $c_s^2={0}$'.format(cs2))
    #figure(c)=full linear spectrum
	plt.figure(c+1)
	#plt.title('Escala log x log')
    #figure(c+1)=full log/log spectrum
    '''
	plt.figure(c+2)
	#plt.title('Diferença relativa para o modelo \u039BCDM')
    #figure(c+2)=residue
    '''
	
	#Para cada cs2, criarei um dataframe pandas no formato | l | ClTT0 | ClTT033 | ClTT067 | ClTT1 onde o numero depois do ClTT indica o valor de wa
	#Esse dataframe vai ser armazenado na lista wa_sims
	
	sim=pd.DataFrame({'l':[], 'ClTT0':[], 'ClTT033':[], 'ClTT067':[], 'ClTT1':[]}, index=None)

	for wa in wa_list:
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
		
		sim['l']=ls
        
        frame11=fig1.add_axes((.1,.3,.8,.6))
        plt.plot(ls[2:], ClTT[2:], label=r'$w_a={0:.2f}$'.format(wa))
		
        frame21=fig2.add_axes((.1,.3,.8,.6))
        plt.plot(ls[2:], ClTT[2:], label=r'$w_a={0:.2f}$'.format(wa))

        #code will probably break at some point after this. I realized I didn't need it for now so I stopped working on it

		if wa==0:
			LCDM=ClTT
			sim['ClTT0']=ClTT
			
		if wa==1/3:
			sim['ClTT033']=ClTT
			
		if wa==2/3:
			sim['ClTT067']=ClTT
			
		if wa==1:
			sim['ClTT1']=ClTT

		dif=[(ClTT[i]-LCDM[i])/LCDM[i] for i in range(len(ClTT))]
        '''
		plt.figure(c+2)
		plt.plot(ls, dif, label=r'$w_a={0:.2f}$'.format(wa))
        '''
	
	sim.to_csv(r'ClTT_for_cs2{0:.2f}.csv'.format(cs2))
			
    #plt.figure(c)
	plt.ylabel(r'$\ell(\ell+1)C_l/2\pi \; (\mu K)^2$')
	plt.xlim([0,max(ls)])
	plt.errorbar(data['l'], data['Dl'], fmt='.',label='Planck Data', yerr=[data['-dDl'], data['+dDl']])
	plt.legend()

    frame2=fig1.add_axes((.1,.1,.8,.2))
	plt.xlabel(r'$\ell$')
    plt.plot(ls, dif)
	#plt.show()
	plt.savefig('DE_cs_fixo{0}.png'.format(cs2))
	
    plt.figure(c+1)
    frame1=fig1.add_axes((.1,.3,.8,.6))
	plt.xlabel(r'$\ell$')
	plt.ylabel(r'$\ell(\ell+1)C_l/2\pi \; (\mu K)^2$')
	plt.errorbar(data['l'], data['Dl'], fmt='.',label='Dados Planck', yerr=[data['-dDl'], data['+dDl']])
	plt.xscale('log')
	plt.yscale('log')
	plt.xlim([1,max(ls)])
	plt.legend()
	#plt.show()
	plt.savefig('DE_cs_fixo{0}log.png'.format(cs2))
	
	plt.figure(c+2)
	plt.xlabel(r'$\ell$')
	plt.ylabel('(ClTT-\u039BCDM)/\u039BCDM')
	plt.xlim([0,max(ls)])
	plt.legend()
	#plt.show()
	plt.savefig('DE_cs_fixo{0}_dif.png'.format(cs2))
	
	print('cs2=', cs2, ' terminado!')
	c+=3

