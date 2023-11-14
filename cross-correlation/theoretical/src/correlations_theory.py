from cobaya.theory import Theory
import numpy as np
from pyctg import ctg4py

'''
In this code I'll define theory classes for the calculation and integration of Ctg and Cgg into Cobaya's sampler.
'''

pkfname = "/home/arthur/Documentos/Mestrado/Projeto/cross-correlation/theoretical/tables/pk_3dmatter.dat"

def ctg4py2(OmegaL, Omegam, lmax, surv_z0, surv_beta, surv_lbda, h, bg, mode, ncalls):

    return ctg4py(OmegaL, Omegam, lmax, surv_z0, surv_beta, surv_lbda, h, bg, mode, ncalls, pkfname)


class ctg(Theory):

    params={"lmax":None, "surv_z0":None, "surv_beta":None, "surv_lbda":None, "bg":None, "mode":1, "ncalls":1000000}

    def initialize(self):
        '''
        Called from __init__ to initialize
        '''

        #self.pk3d=np.loadtxt("/home/arthur/Documentos/Mestrado/Projeto/cross-correlation/theoretical/tables/pk_3dmatter.dat")

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

    def get_requirements(self):
        '''
        Returns a dictionary of quantities that are needed by this class, these can be calculated by another theory class or just be samples parameters
        '''

        return {'omega_de':None, 'ombh2':None, 'omch2':None, 'H0':None} # omega_de=OmegaL em LCDM (?)

    def calculate(self, state, want_derived=False, **params_values_dict):
        '''
        This function is used to calculate and store the results in the "state" dictionary
        '''

        lmax, bg, mode, ncalls = self.provider.get_param('lmax'), self.provider.get_param('bg'), self.provider.get_param('mode'), self.provider.get_param('ncalls')
        h=self.provider.get_param('H0')/100
        OmegaL, Omegam=self.provider.get_param('omega_de'), (self.provider.get_param('ombh2')+self.provider.get_param('omch2'))/h**2
        z0, beta, lbda=self.provider.get_param('surv_z0'), self.provider.get_param('surv_beta'), self.provider.get_param('surv_lbda')

        ctg=[]

        for l in range(round(lmax)):
            cl=ctg4py2(OmegaL, Omegam, lmax, z0, beta, lbda, h, bg, mode, ncalls)
            ctg.append(cl)

        state['ctg']=cl

    def get_ctg(self):
        '''
        This function should return the current value of ctg stored in the "state" dictionary
        '''

        return self.current_state['ctg']

#depois que ctg estiver testado, eu preparo cgg

"""
class cgg(Theory):

    params={"lmax":None, "surv_z0":None, "surv_beta":None, "surv_lbda":None, "bg":None, "mode":1, "ncalls":1000000, "fname":pkfname}

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

        return ['cgg']

    def get_requirements(self):
        '''
        Returns a dictionary of quantities that are needed by this class, these can be calculated by another theory class or just be samples parameters
        '''

        return {}

    def calculate(self, state, want_derived=False, **params_values_dict):
        '''
        This function is used to calculate and store the results in the "state" dictionary
        '''

        cgg=0

    def get_cgg(self):
        '''
        This function should return the current value of ctg stored in the "state" dictionary
        '''

        return self.current_state['cgg']
"""
