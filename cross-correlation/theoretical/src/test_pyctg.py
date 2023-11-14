from correlations_theory import ctg 
import numpy as np
import random
from pyctg import ctg4py
import ctypes

OmegaL, Omegam, h = 0.7, 0.3, 0.67
bg, mode, ncalls, lmax = 1.0, 1, 100000, 32 
beta, lbda, z0 = 3.09, 4.94, 0.15
fname="/home/arthur/Documentos/Mestrado/Projeto/cross-correlation/theoretical/tables/pk_3dmatter.dat"
fname=ctypes.c_char_p(fname.encode("ascii"))
print("ctg(lmax)=", ctg4py(OmegaL, Omegam, lmax, z0, beta, lbda, h, bg, mode, ncalls, fname))

def correlation_like(_self=None):
    print("inside correlation!")
    
    #Descobrir como flexibilizar: Deixar o usuário escolher se quer ctg ou cgg
    cl_theo=_self.provider.get_ctg()
    print("theo=", cl_theo)
    #Por enquanto usando dados fictícios:
    cl_sigmas=[1/10 for i in range(2,lmax)]
    cl_data=[(ctg4py(OmegaL, Omegam, l, z0, beta, lbda, h, bg, mode, ncalls, fname)+random.gauss(0., 0.1)) for l in range(2,lmax)] #importar os dados
    print("data=", cl_data)
    print("sigmas=", cl_sigmas)

    theory_array, data_array, sigma_array=np.array(cl_theo), np.array(cl_data), np.array(cl_sigmas)
    
    sum1=sum(2*np.pi*sigma_array**2)
    sum2=sum(((theory_array-data_array)/sigma_array)**2)

    logp=-0.5*(np.log(sum1)+sum2)

    print('log(p)=',logp)

    return logp    

if __name__=='__main__':
    #Testing area

    info={
         'params':{
             #Fixed 
             'lmax':32,
             'bg':1.,
             'surv_z0':{'ref':1., 'prior':{'min':0.0001, 'max':1.0}},#0.15,
             'surv_beta':{'ref':1., 'prior':{'min':1., 'max':4.}},#3.1,
             'surv_lbda':{'ref':1., 'prior':{'min':1., 'max':6.}}#4.9,
             },
         'likelihood':{
             'my_like':{'external':correlation_like, 'requires':{'ctg':None}}
             },
         'theory':{'test_pyctg_theory.ctg':None}
         }

    from cobaya.model import get_model
    model=get_model(info)

    from cobaya.run import run
    info['sampler'] = {'mcmc':{'Rminus1_stop':0.1, 'max_tries':5000}}
    updated_info, sampler=run(info)

    from getdist.mcsamples import MCSamplesFromCobaya
    import getdist.plots as gdplt
    import matplotlib.pyplot as plt

    gdsamples = MCSamplesFromCobaya(updated_info, sampler.products()["sample"])
    gdplot = gdplt.get_subplot_plotter(width_inch=5)
    gdplot.triangle_plot(gdsamples, ['surv_z0', 'surv_beta', 'surv_lbda'], filled=True)
    gdplot = gdplt.get_subplot(width_inch=5)
    plt.savefig("my_samples.png")
