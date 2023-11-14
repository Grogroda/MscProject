import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

name="prior_lcdm_halofit_xsc_fixed_monodip_and_above51_fittedbgl50_miliK_cross_band"

band1=pd.read_csv(name+"1.dat", header=None, sep=' ', names=['l', 'ctt', 'ctt_include', 'cgg', 'cgg_include', 'ctg', 'ctg_include'])
band2=pd.read_csv(name+"2.dat", header=None, sep=' ', names=['l', 'ctt', 'ctt_include', 'cgg', 'cgg_include', 'ctg', 'ctg_include'])
band3=pd.read_csv(name+"3.dat", header=None, sep=' ', names=['l', 'ctt', 'ctt_include', 'cgg', 'cgg_include', 'ctg', 'ctg_include'])
band4=pd.read_csv(name+"4.dat", header=None, sep=' ', names=['l', 'ctt', 'ctt_include', 'cgg', 'cgg_include', 'ctg', 'ctg_include'])

bands=[band1, band2, band3, band4]

for i in range(len(bands)):
    ls=bands[i]['l']
    ctt=bands[i]['ctt']
    cgg=bands[i]['cgg']
    ctg=bands[i]['ctg']

    kctt=[ls[k]*(ls[k]+1)/(2*np.pi)*ctt[k] for k in range(len(ctt))]

    plt.figure()
    plt.title('Ctt for band {}'.format(i+1))
    plt.plot(ls[2:51], kctt[2:51], label='Sampled')
    plt.plot(ls[51:], kctt[51:], label='Theoretical')
    plt.legend()
    plt.xlabel(r'$\ell$')
    plt.ylabel(r'$\ell(\ell+1)/2\pi C_\ell^{tt} [mK^2]$')
    plt.xscale('log')
    plt.savefig('ctt_band{}.png'.format(i+1), bbox_inches = 'tight')
    plt.close()

    plt.figure()
    plt.title('Cgg for band {}'.format(i+1))
    plt.plot(ls[2:51], cgg[2:51], label='Sampled')
    plt.plot(ls[51:], cgg[51:], label='Theoretical')
    plt.legend()
    plt.xlabel(r'$\ell$')
    plt.ylabel(r'$C_\ell^{gg}$')
    plt.xscale('log')
    plt.yscale('log')
    plt.savefig('cgg_band{}.png'.format(i+1), bbox_inches = 'tight')
    plt.close()

    plt.figure()
    plt.title('Ctg for band {}'.format(i+1))
    plt.plot(ls[2:51], ctg[2:51], label='Sampled')
    plt.plot(ls[51:], ctg[51:], label='Theoretical')
    plt.legend()
    plt.xlabel(r'$\ell$')
    plt.ylabel(r'$C_\ell^{tg} [mK]$')
    plt.xscale('log')
    plt.savefig('ctg_band{}.png'.format(i+1), bbox_inches = 'tight')
    plt.close()

