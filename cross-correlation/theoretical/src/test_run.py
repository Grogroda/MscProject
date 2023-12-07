from correlations_theory import ctg 
import numpy as np
import random
from pyctg import ctg4py
import ctypes

#OmegaL, Omegam, h = 0.7, 0.3, 0.67
#bg, mode, ncalls, lmax = 1.0, 1, 100000, 32 
#beta, lbda, z0 = 3.09, 4.94, 0.15
#fname="../tables/pk_3dmatter.dat"
#fname=ctypes.c_char_p(fname.encode("ascii"))
#print("ctg(lmax)=", ctg4py(OmegaL, Omegam, lmax, z0, beta, lbda, h, bg, mode, ncalls, fname))


#Fiducial model to generate data with 10% precision. Has to be outside the likelihood function, otherwise will be recomputed
# for each point. 
OmegaM_fid = 0.3
cl_data  = ctg4py(OmegaM_fid)
cl_sigmas= [i*0.1 for i in cl_data]

def correlation_like(_self=None):    
    #Descobrir como flexibilizar: Deixar o usuário escolher se quer ctg ou cgg
    cl_theo=_self.provider.get_ctg()
#    print(cl_theo,cl_data,cl_sigmas)

    theory_array, data_array, sigma_array=np.array(cl_theo), np.array(cl_data), np.array(cl_sigmas)

# Old form...   
#    sum1=sum(2*np.pi*sigma_array**2)
#    sum2=sum(((theory_array-data_array)/sigma_array)**2)
#    logp=-0.5*(np.log(sum1)+sum2)

    Xi2    = np.sum(((theory_array-data_array)/sigma_array)**2)
    sigSum = np.sum(np.log(sigma_array))
    # Constant not included.
    logp = - sigSum - Xi2/2

    return logp    

if __name__=='__main__':
    #Testing area

    info={
         'params':{
             'OmegaM':{'ref':0.3, 'prior':{'min':0.01, 'max':0.99}}
             },
         'likelihood':{
             'my_like':{'external':correlation_like, 'requires':{'ctg':None}}
             },
         'theory':{'test_pyctg_theory.ctg':None}
         }

    from cobaya.model import get_model
    model=get_model(info)

    from cobaya.run import run
    info['sampler'] = {'mcmc':{'Rminus1_stop':0.001, 'max_tries':10000}}
    updated_info, sampler=run(info)

    from getdist.mcsamples import MCSamplesFromCobaya
    import getdist.plots as gdplt
    import matplotlib.pyplot as plt

    gdsamples = MCSamplesFromCobaya(updated_info, sampler.products()["sample"])
    gdplot = gdplt.get_subplot_plotter(width_inch=5)
    gdplot.triangle_plot(gdsamples, ['OmegaM'], filled=True)
#    gdplot = gdplt.get_subplot(width_inch=5)
    plt.savefig("my_samples.png")
