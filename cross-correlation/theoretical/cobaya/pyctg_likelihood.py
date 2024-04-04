# -*- coding: utf-8 -*-
import numpy as np
import random
from pyctg import ctg4py
from pyctg import cgg4py
import ctypes
import matplotlib.pyplot as plt

#Fiducial model to generate data with 10% precision. Has to be outside the likelihood function, otherwise will be recomputed for each point
OmegaM_fid = 0.3
ls_ctg, ctg_data  = ctg4py(OmegaM_fid)
ctg_sigmas= [i*0.1 for i in ctg_data]
print("fake ctg calculated")

'''
ls, cgg_data = cgg4py(OmegaM_fid)
cgg_sigmas = [i*0.1 for i in cgg_data]
'''
ls_cgg, cgg_data, cgg_sigmas=[1],[1],[1]

#cl_2nd = ctg4py(0.42)
#ls=[i for i in range(len(cl_data))]

test_plots=True

def ctg_like(_self=None):    
    #Descobrir como flexibilizar: Deixar o usuário escolher se quer ctg ou cgg
    cl_theo=_self.provider.get_ctg()
#    print(cl_theo,cl_data,cl_sigmas)
    print("cl_theo=", cl_theo)

    theory_array, data_array, sigma_array=np.array(cl_theo), np.array(ctg_data), np.array(ctg_sigmas)

    #Chi2 calculation:
    Xi2    = np.sum(((theory_array-data_array)/sigma_array)**2)
    print("Xi2=", Xi2)
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

print("likelihoods defined")

def estimate_proposal(min_max, n=25, plot=True):
    '''
    min_max (dict)=dictionary with parameter name (key) associated with a size 2 array (tuple or list) indication par_min and par_max
    This function returns a diagonal matrix with proposal variances for the ith parameter in the element [i][i]
    To do this, the function calculates the likelihood of the parameters for values of that parameter between par_min and par_max (without varying the others) and calculates the variance of the distribution. If plot=True, it plots the distribution so the user can analyse it if necessary. 
    n (int)=number of points to for which L is being calculated for each parameter.
    '''

    expected_pars=['OmegaM']
    #The output matrix is ordered according to expected_pars, not min_max keys order. Make this very clear for the user
    default_values={'OmegaM':0.31} #use middle value as default instead?
    pars={'OmegaM':0.31}

    proposal_matrix=[]
    npars=len(expected_pars)

    if len(list(min_max.keys()))!=len(expected_pars):
        raise Exception("Expected {0} parameters, {1} were passed. List of expected parameters: {2}".format(len(expected_pars),len(list(min_max.keys())), expected_pars))
    else:
        counter=0
        for par in expected_pars:
            #checks if the user passed a wrong parameter:
            if par not in min_max:
                raise Exception("Parameter {0} not expected. List of expected parameters: {1}".format(par, expected_pars))

            proposal_matrix.append([0 for i in range(npars)]) #creates line with npars zeroes
            #Sets the current parameter to its minimum value:
            pars[par]=min_max[par][0]
            delta_par=(min_max[par][1]-min_max[par][0])/n
            
            #Sets all parameter values
            OmegaM=pars['OmegaM']

            xs, ys=[],[]

            for i in range(n+1):
                print("i=", i)
                ctg_theo=ctg4py(OmegaM)

                theory_array, data_array, sigma_array=np.array(ctg_theo), np.array(cgg_data), np.array(cgg_sigmas)

                #Chi2 calculation:
                Xi2    = np.sum(((theory_array-data_array)/sigma_array)**2)
                sigSum = np.sum(np.log(sigma_array))

                xs.append(pars[par])
                ys.append(-Xi2-sigSum)

                pars[par]+=delta_par

            pars[par]=default_values[par]

            plt.figure()
            plt.plot(xs, ys)
            plt.xaxis(par)
            plt.yaxis(r'$\mathcal{L}=log(p)$')
            plt.savefig('LProfile_{}.png'.format(par))

    return None #como retornar a matriz "exata"?

if __name__=='__main__':
    #Testing area

    '''
    print("Inside testing area")
    estimate_proposal({'OmegaM':(0.01, 0.99)})
    '''

    if test_plots:
        plt.figure()
        plt.plot(ls_ctg, ctg_data, label='Fake data (OmegaM={})'.format(OmegaM_fid))
        plt.ylabel(r'$C^{tg}$')
        plt.xscale('log')
        plt.xlabel(r'$\ell$')
        plt.savefig('Ctg_RunPlot.png')

        plt.figure()
        plt.plot(ls_cgg, cgg_data, label='Fake data (OmegaM={})'.format(OmegaM_fid))
        plt.ylabel(r'$C^{gg}$')
        plt.xscale('log')
        plt.yscale('log')
        plt.xlabel(r'$\ell$')
        plt.savefig('Cgg_RunPlot.png')

        print("Test plots completed")

    info={
         'params':{
             'OmegaM':{'ref':0.3, 'prior':{'min':0.05, 'max':0.80}}
             },
         'likelihood':{
             'ctg_like':{'external':ctg_like, 'requires':{'ctg':None}}
             },
         'theory':{'pyctg_theory.ctg':None}#, 'pyctg_theory.cgg':None}
         }

    print("info dict defined")

    from cobaya.model import get_model

    model=get_model(info)

    print("Model created")

    from cobaya.run import run
    info['sampler'] = {'mcmc':{'Rminus1_stop':0.001, 'max_tries':10000}}
    info['output'] = "scriptrun_test/ctg_Pkarray_lmax54_N2e5"
    print("Run about to start")
    updated_info, sampler=run(info, resume=True, no_mpi=False)

    from getdist.mcsamples import MCSamplesFromCobaya
    import getdist.plots as gdplt

    gdsamples = MCSamplesFromCobaya(updated_info, sampler.products()["sample"])
    gdplot = gdplt.get_subplot_plotter(width_inch=5)
    gdplot.triangle_plot(gdsamples, ['OmegaM'], filled=True)
#    gdplot = gdplt.get_subplot(width_inch=5)
    plt.savefig("scriptrun_test/ctg_Pkarray_lmax54_N2e5.png")
