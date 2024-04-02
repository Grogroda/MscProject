import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

#dict={n_band:{col_name:[values], ...}, ...}
band_dict={}

arr_names=['ls', 'Dtt', 'err_Dtt', 'Cgg', 'err_Cgg', 'Ctg', 'err_Ctg']

band=1
while band<=4:
    dict_band=pd.read_csv('errors_data_wmap9QVW_xsc{}_jeffrey_dipfix_50000_ttggtg_new.dat'.format(band), header=None, names=arr_names, sep=' ')
    band_dict[band]=dict_band
    band+=1

#print(band_dict)

import matplotlib

matplotlib.rcParams.update({'font.size': 15})

fig, axs = plt.subplots(4,3, figsize=(20, 16))

band=1
while band<=4:
    print("Starting band {}".format(band))
    #Dtt plot for each band:
    axs[band-1,0].errorbar(band_dict[band]['ls'], band_dict[band]['Dtt'], yerr=band_dict[band]['err_Dtt'], label="Best-fit+total error", fmt='ko')
    axs[band-1,0].set_ylabel(r'$(\ell+1)\ell/2\pi C_\ell^{tt}$ $[\mu K^2]$', fontsize=20)
    axs[band-1,0].set_xlabel(r'$\ell$', fontsize=20)
    axs[band-1,0].set_xscale('log')
    axs[band-1,0].set_ylim(0, 4000)
    axs[band-1,0].legend()
    print("Ctt finished")

    #Cgg plot for each band:
    axs[band-1,1].errorbar(band_dict[band]['ls'], band_dict[band]['Cgg'], yerr=band_dict[band]['err_Cgg'], fmt='ko')
    axs[band-1,1].set_xlabel(r'$\ell$', fontsize=20)
    axs[band-1,1].set_ylabel(r'$C_\ell^{gg}$', fontsize=20)
    axs[band-1,1].set_ylim(1e-4, 1e-1)
    axs[band-1,1].set_xscale('log')
    axs[band-1,1].set_yscale('log')
    print('Cgg finished')

    #Ctg plot for each band:
    axs[band-1,2].errorbar(band_dict[band]['ls'], band_dict[band]['Ctg'], yerr=band_dict[band]['err_Ctg'], fmt='ko')
    axs[band-1,2].set_xlabel(r'$\ell$', fontsize=20)
    axs[band-1,2].set_ylabel(r'$C_\ell^{tg}$ $[\mu K]$', fontsize=20)
    axs[band-1,2].set_ylim(-0.8, 0.8)
    axs[band-1,2].set_xscale('log')
    #axs[band-1,2].set_yscale('log')
    axs[band-1,2].text(32, 0.6, 'Band {}'.format(band))
    print("Band {} finished, moving to next".format(band))
    band+=1

fig.tight_layout()
plt.savefig("Full_Data_Plot.png")

