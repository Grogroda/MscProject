from correlations_theory import ctg 
import numpy as np
import random

#pkfname = "/home/arthur/Documentos/Mestrado/Projeto/cross-correlation/theoretical/tables/pk_3dmatter.dat"

def correlation_like(_self=None):
    
    #Descobrir como flexibilizar: Deixar o usuário escolher se quer ctg ou cgg
    cl_theo=_self.provider.get_ctg()
    print("theo=", cl_theo)
    #Por enquanto usando dados fictícios:
    cl_data=[C*(1+random.uniform(-0.2,0.2)) for C in cl_theo] #importar os dados
    print("data=", cl_data)
    cl_sigmas=[C/10 for C in cl_data]
    print("sigmas=", cl_sigmas)

    theory_array, data_array, sigma_array=np.array(cl_theo), np.array(cl_data), np.array(cl_sigmas)
    
    sum1=sum(2*np.pi*sigma_array**2)
    sum2=sum(((theory_array-data_array)/sigma_array)**2)

    logp=-0.5*(np.log(sum1)+sum2)

    print('log(p)=',logp)

    return logp    

if __name__=='__main__':
    #Testing area

    import matplotlib.pyplot as plt

    info = {
            #Dict to sample the parameters
            'params':{
                #Fixed first
                'lmax':32,
                #'surv_z0':0.15,
                #'surv_beta':3.1,
                #'surv_lbda':4.9,
                'bg':1,
                #'fname':pkfname,
                'ombh2':0.02237,
                'omch2':0.12,
                #Now sampled
                'H0':{'ref':70.0, 'prior':{'min':30.0, 'max':100.0}},
                'omega_de':{'ref':0.69, 'prior':{'min':0.0, 'max':1.0}}
                },
            'likelihood':{
                'my_like':{'external':correlation_like, 'requires':{'ctg':None}},
                'planck_2018_lowl.TT': None},
            'theory':{
                'correlations_theory.ctg':None,
                'camb': {'extra_args': {'bbn_predictor': 'PArthENoPE_880.2_standard.dat',
                    'halofit_version': 'mead',
                    'lens_potential_accuracy': 1,
                    'nnu': 3.044,
                    'num_massive_neutrinos': 1,
                    'theta_H0_range': [20, 100]}}}}

    from cobaya.model import get_model
    model=get_model(info)

    from cobaya.run import run
    info['sampler'] = {'mcmc':{"Rminus1_stop": 0.1, "max_tries": 9000}}
    updated_info, sampler=run(info)

    from getdist.mcsamples import MCSamplesFromCobaya
    import getdist.plots as gdplt
    import matplotlib.pyplot as plt

    gdsamples = MCSamplesFromCobaya(updated_info, sampler.products()["sample"])
    gdplot = gdplt.get_subplot_plotter(width_inch=5)
    gdplot.triangle_plot(gdsamples, ["H0", "omega_de"], filled=True)
    gdplot = gdplt.get_subplot(width_inch=5)
    plt.savefig("my_samples.png")
