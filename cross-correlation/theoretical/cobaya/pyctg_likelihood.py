import numpy as np
import random
from pyctg import ctg4py
from pyctg import cgg4py
import ctypes
import matplotlib.pyplot as plt

#Fiducial model to generate data with 10% precision. Has to be outside the likelihood function, otherwise will be recomputed
# for each point. 
OmegaM_fid = 0.3
ls, ctg_data  = ctg4py(OmegaM_fid)
ctg_sigmas= [i*0.1 for i in ctg_data]

ls, cgg_data = cgg4py(OmegaM_fid)
cgg_sigmas = [i*0.1 for i in cgg_data]

#cl_2nd = ctg4py(0.42)
#ls=[i for i in range(len(cl_data))]

test_plots=True

def ctg_like(_self=None):    
    #Descobrir como flexibilizar: Deixar o usuário escolher se quer ctg ou cgg
    cl_theo=_self.provider.get_ctg()
#    print(cl_theo,cl_data,cl_sigmas)

    theory_array, data_array, sigma_array=np.array(cl_theo), np.array(ctg_data), np.array(ctg_sigmas)

    #Chi2 calculation:
    Xi2    = np.sum(((theory_array-data_array)/sigma_array)**2)
    sigSum = np.sum(np.log(sigma_array))
    # Constant not included.
    logp = - sigSum - Xi2/2

    return logp    

def cgg_like(_self=None):
    #Descobrir como flexibilizar: Deixar o usuário escolher se quer ctg ou cgg
    cl_theo=_self.provider.get_cgg()
#    print(cl_theo,cl_data,cl_sigmas)

    theory_array, data_array, sigma_array=np.array(cl_theo), np.array(cgg_data), np.array(cgg_sigmas)

    #Chi2 calculation:
    Xi2    = np.sum(((theory_array-data_array)/sigma_array)**2)
    sigSum = np.sum(np.log(sigma_array))
    # Constant not included.
    logp = - sigSum - Xi2/2

    return logp

if __name__=='__main__':
    #Testing area

    if test_plots:
        plt.figure()
        plt.plot(ls, ctg_data, label='Fake data (OmegaM={})'.format(OmegaM_fid))
        plt.ylabel(r'$C^{tg}$')
        plt.xscale('log')
        plt.xlabel(r'$\ell$')
        plt.savefig('Ctg_RunPlot.png')

        plt.figure()
        plt.plot(ls, cgg_data, label='Fake data (OmegaM={})'.format(OmegaM_fid))
        plt.ylabel(r'$C^{gg}$')
        plt.xscale('log')
        plt.yscale('log')
        plt.xlabel(r'$\ell$')
        plt.savefig('Cgg_RunPlot.png')

    info={
         'params':{
             'OmegaM':{'ref':0.3, 'prior':{'min':0.01, 'max':0.99}}
             },
         'likelihood':{
             'ctg_like':{'external':ctg_like, 'requires':{'ctg':None}},
             'cgg_like':{'external':cgg_like, 'requires':{'cgg':None}}
             },
         'theory':{'pyctg_theory.ctg':None, 'pyctg_theory.cgg':None}
         }

    from cobaya.model import get_model
    model=get_model(info)

    from cobaya.run import run
    info['sampler'] = {'mcmc':{'Rminus1_stop':0.001, 'max_tries':10000}}
    info['output'] = "scriptrun_test/ctg_cgg_lmax54_ncalls5e5"
    updated_info, sampler=run(info, resume=True, no_mpi=False)

    from getdist.mcsamples import MCSamplesFromCobaya
    import getdist.plots as gdplt

    gdsamples = MCSamplesFromCobaya(updated_info, sampler.products()["sample"])
    gdplot = gdplt.get_subplot_plotter(width_inch=5)
    gdplot.triangle_plot(gdsamples, ['OmegaM'], filled=True)
#    gdplot = gdplt.get_subplot(width_inch=5)
    plt.savefig("scriptrun_test/ctg_cgg_lmax54_ncalls5e5.png")
