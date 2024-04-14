import numpy as np
from matplotlib import pyplot as plt
import pandas as pd

xi_data=pd.read_csv('xi_list_10M.dat', header=None, names=['OmegaM', 'logP'], sep=' ')

Omegas=xi_data['OmegaM']

minus_logp=[(-1)*logp for logp in xi_data['logP']]

max_minus_logp=max(minus_logp)
rel_logp=[logp-max_minus_logp for logp in minus_logp]
P=[np.exp(logp) for logp in rel_logp]

plt.figure()
plt.plot(Omegas, P)
#plt.show()
plt.xlabel(r'$\Omega_m$')
plt.ylabel(r'$\mathcal{L}(\Omega_m)$')
plt.savefig('OmegaM_prof_data.png')
