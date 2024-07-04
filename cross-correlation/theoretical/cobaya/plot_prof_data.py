import numpy as np
from matplotlib import pyplot as plt
import pandas as pd

import matplotlib

matplotlib.rcParams.update({'font.size':15})

xi_data=pd.read_csv('xi_list_bestband_Om0_188_merged.dat', header=None, names=['OmegaM', 'logP'], sep=' ')
fid_best_band=pd.read_csv('xi_list_band_best_new.dat', header=None, names=['OmegaM', 'logP'], sep=' ')
band4_data=pd.read_csv('like_profiles/ctg_N2e7/full_data4.dat', sep=' ', index_col=0)

Omegas_best=xi_data['OmegaM']
minus_logp=[(-1)*logp for logp in xi_data['logP']]

Omegas_fid=fid_best_band['OmegaM']
logp_fid=[(-1)*logp for logp in fid_best_band['logP']]

Omegas_band4=band4_data['Omegas']
raw_logp=band4_data['raw_logp']

max_minus_logp=max(minus_logp)
rel_logp=[logp-max_minus_logp for logp in minus_logp]
P=[np.exp(logp) for logp in rel_logp]

max_logp_fid=max(logp_fid)
diff_logp_fid=[logp_fid[i]-max_logp_fid for i in range(len(logp_fid))]
P_fid=[np.exp(logp) for logp in diff_logp_fid]

max_logp4=max(raw_logp)
diff_logp=[raw_logp[i]-max_logp4 for i in range(len(raw_logp))]
exp_diff=[np.exp(logp) for logp in diff_logp]

plt.figure()
plt.plot(Omegas_best, P, label=r'Optimized $\Omega_m=0.188$', linestyle="dashed")
plt.plot(Omegas_fid, P_fid, label="Optimized Fiducial")
plt.plot(Omegas_band4, exp_diff, label='Band 4')
planck_bf, error=0.3135, 0.0073
plt.vlines(planck_bf, 0, 2, linestyle='dashdot', colors='k', label='Planck best-fit')
plt.axvspan(planck_bf-error, planck_bf+error, alpha=0.5, color='gray')
#plt.show()
plt.xlim(0.05378, 1.02)
plt.ylim(0.5, 1.03)
plt.xlabel(r'$\Omega_m$')
plt.ylabel(r'$\mathcal{L}/\mathcal{L}_\text{max}$')
plt.legend(fontsize=12)
plt.tight_layout()
plt.savefig('OmegaM_prof_BestBand_Om0188.png')
