import numpy as np
from matplotlib import pyplot as plt
import pandas as pd

import matplotlib

matplotlib.rcParams.update({'font.size':15})

xi_data=pd.read_csv('xi_list_band_best.dat', header=None, names=['OmegaM', 'logP'], sep=' ')

Omegas=xi_data['OmegaM']

minus_logp=[(-1)*logp for logp in xi_data['logP']]

max_minus_logp=max(minus_logp)
rel_logp=[logp-max_minus_logp for logp in minus_logp]
P=[np.exp(logp) for logp in rel_logp]

plt.figure()
plt.plot(Omegas, P)
plt.vlines(0.3135, 0, 2, linestyle='dashed', colors='k', label='Planck best-fit')
#plt.show()
plt.ylim(0.84, 1.01)
plt.xlabel(r'$\Omega_m$')
plt.ylabel(r'$\mathcal{L}/\mathcal{L}_\text{max}$')
plt.legend()
plt.tight_layout()
plt.savefig('OmegaM_prof_BestBand.png')
