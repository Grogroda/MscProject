import numpy as np
from matplotlib import pyplot as plt
import pandas as pd
import matplotlib

matplotlib.rcParams.update({'font.size':15})

full_data={1:{}, 2:{}, 3:{}, 4:{}} #{band:{'Omegas':[], raw_logp:[], raw_xi2:[]}}

for band in list(full_data.keys()):
    print('Band: ', band)
    Omegas, raw_logp, raw_xi2=[],[],[]

    for i in [1,2,3,4]:
        data=pd.read_csv('ProfileData{0}_both_band{1}_Nmc2e+07.dat'.format(i, band), sep=' ', index_col=0)
        Omegas+=list(data['Omega'])
        raw_logp+=list(data['raw_logp'])
        raw_xi2+=list(data['raw_xi2'])

    full_data[band]['Omegas']=Omegas
    full_data[band]['raw_logp']=raw_logp
    full_data[band]['raw_xi2']=raw_xi2

    max_logp=max(raw_logp)
    diff_logp=[raw_logp[i]-max_logp for i in range(len(raw_logp))]
    exp_pdiff=[np.exp(logp) for logp in diff_logp]

    plt.figure()
    plt.title('Band {}'.format(band))
    plt.xlabel(r'$\Omega_m$')
    plt.ylabel(r'$P(\Omega_m)$')
    plt.plot(Omegas, exp_pdiff)
    plt.savefig('profile_band{0}_Nmc2e7.png'.format(band))

print('Data dict:\n', full_data)

'''
max_logp=max(raw_logp1)
diff_logp=[raw_logp1[i]-max_logp for i in range(len(raw_logp1))]
exp_pdiff=[np.exp(logp) for logp in diff_logp]

plt.figure()
plt.plot(Omegas1, exp_pdiff)
plt.show()
'''  
