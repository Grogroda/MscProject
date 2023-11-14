from cobaya.theory import Theory
import numpy as np
from pyctg import ctg4py

'''
In this code I'll define theory classes for the calculation and integration of Ctg and Cgg into Cobaya's sampler.
'''

pkfname = "/home/arthur/Documentos/Mestrado/Projeto/cross-correlation/theoretical/tables/pk_3dmatter.dat"
OmegaL, Omegam=0.7, 0.3

def ctg4py2(l, surv_z0, surv_beta, surv_lbda, h, bg, mode, ncalls):

    return ctg4py(OmegaL, Omegam, l, surv_z0, surv_beta, surv_lbda, h, bg, mode, ncalls, pkfname)

class ctg(Theory):

    params={"lmax":None, "surv_z0":None, "surv_beta":None, "surv_lbda":None, "bg":None, "h":0.67, "mode":1, "ncalls":1000000}

    def initialize(self):
        '''
        Called from __init__ to initialize
        '''

    def initialize_with_provider(self, provider):
        '''
        Initialization after other components are initialized, using Provider class instance, which is used to return any dependencies
        '''

        self.provider=provider

    def get_can_provide(self):
        '''
        Outputs a list of parameters that can be calculated within this class
        '''

        return ['ctg']

    def calculate(self, state, want_derived=False, **params_values_dict):
        '''
        This function is used to calculate and store the results in the "state" dictionary
        '''

        lmax, bg, mode, ncalls = self.provider.get_param('lmax'), self.provider.get_param('bg'), self.provider.get_param('mode'), self.provider.get_param('ncalls')
        h = self.provider.get_param('h')
        z0, beta, lbda=self.provider.get_param('surv_z0'), self.provider.get_param('surv_beta'), self.provider.get_param('surv_lbda')

        ctg=[]

        for l in range(2, round(lmax)):
            cl=ctg4py2(l, z0, beta, lbda, h, bg, mode, ncalls)
            ctg.append(cl)

        state['ctg']=cl

    def get_ctg(self):
        '''
        This function should return the current value of ctg stored in the "state" dictionary
        '''

        return self.current_state['ctg']

