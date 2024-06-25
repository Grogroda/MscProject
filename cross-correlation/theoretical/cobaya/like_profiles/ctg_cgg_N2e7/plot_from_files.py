import numpy as np
from matplotlib import pyplot as plt
import pandas as pd
import matplotlib

matplotlib.rcParams.update({'font.size':15})

full_data={1:{}, 2:{}, 3:{}, 4:{}} #{band:{'Omegas':[], raw_logp:[], raw_xi2:[]}}

plt.figure()
plt.xlabel(r'$\Omega_m$')
plt.ylabel(r'$\mathcal{L}/\mathcal{L}_\text{max}$')
plt.tight_layout()

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

    plt.plot(Omegas, exp_pdiff, label='Band {}'.format(band))


print('Data dict:\n', full_data)

sum_table={'Omegas':[]} #table containing the sums of raw_logp and raw_xi2

sum_table['Omegas']=full_data[1]['Omegas']
sum_table['raw_logp']=[0 for i in range(len(sum_table['Omegas']))]
sum_table['raw_xi2']=[0 for i in range(len(sum_table['Omegas']))]

for i in range(len(sum_table['raw_logp'])):
    print('i=', i)
    for band in [1,2,3,4]:
        print('Band=', band)
        sum_table['raw_logp'][i]+=full_data[band]['raw_logp'][i]
        sum_table['raw_xi2'][i]+=full_data[band]['raw_xi2'][i]

max_logp=max(sum_table['raw_logp'])
diff_logp=[sum_table['raw_logp'][i]-max_logp for i in range(len(sum_table['raw_logp']))]
exp_pdiff=[np.exp(logp) for logp in diff_logp]

plt.plot(sum_table['Omegas'], exp_pdiff, label='All bands', linestyle='dashed')
planck_bf, error=0.3135, 0.0073
plt.vlines(planck_bf, -0.05, 2, linestyle='dashdot', colors='k', label='Planck best-fit')
plt.axvspan(planck_bf-error, planck_bf+error, alpha=0.5, color='gray')
plt.ylim(-0.02,1.02)
plt.xlim(min(sum_table['Omegas']), max(sum_table['Omegas']))

plt.legend()
plt.savefig('profile_allbands_Nmc2e7.png')
