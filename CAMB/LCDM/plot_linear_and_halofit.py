import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

import matplotlib
matplotlib.rcParams.update({'font.size':16})

PS_linear=pd.read_csv('matterPS_linear_lmax128.dat', sep=' ', header=None, names=['kh', 'Pk'])
PS_halofit=pd.read_csv('matterPS_halofit_lmax128.dat', sep=' ', header=None, names=['kh', 'Pk'])

kh_lin, pk_lin=PS_linear['kh'], PS_linear['Pk']
kh_halo, pk_halo=PS_halofit['kh'], PS_halofit['Pk']

plt.figure()
plt.plot(kh_lin, pk_lin, label='Linear')
plt.plot(kh_halo, pk_halo, label='Halofit')
plt.xscale('log')
plt.yscale('log')
plt.xlabel(r'$k/h \text{ } [\text{h}^{-1}\text{Mpc}^{-1}]$')
plt.ylabel(r'$P(k/h) \text{ }[\text{h}^{-3}\text{Mpc}^{-3}]$')
plt.legend()
plt.tight_layout()
plt.savefig('Pk_comparison.png')
