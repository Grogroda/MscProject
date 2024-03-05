import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

ctt_data=pd.read_csv("ctt_LCDM.dat", sep=' ', names=['l', 'cltt'])

#print(ctt_data)
ls=ctt_data['l']
Dtt=[ls[i]*(ls[i]+1)/(2*np.pi)*ctt_data['cltt'][i] for i in range(len(ls))]

plt.figure()
plt.plot(ls[2:], Dtt[2:])
plt.xscale('log')
plt.yscale('log')
plt.xlabel(r'$\ell$')
plt.ylabel(r'$\ell(\ell+1)/2\pi C_\ell^{tt}$')
plt.savefig('ctt_LCDM.png')
