from cobaya.theory import Theory
import numpy as np
from pyctg import ctg4py

'''
In this code I'll define theory classes for the calculation and integration of Ctg and Cgg into Cobaya's sampler.
'''

class ctg(Theory):

    params={"OmegaM":None}

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
        OmegaM = self.provider.get_param('OmegaM')
        cl     = ctg4py(OmegaM)
        state['ctg']=cl

    def get_ctg(self):
        '''
        This function should return the current value of ctg stored in the "state" dictionary
        '''

        return self.current_state['ctg']

