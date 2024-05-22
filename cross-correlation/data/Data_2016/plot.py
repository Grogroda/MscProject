import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

Ctt_etc=pd.read_csv('../../theoretical/tables/final_table.dat', sep=' ', header=None, names=['ls', 'ctt', 'cgg', 'ctg'])
ls_ctt=Ctt_etc['ls'][2:]
ctt_theo=Ctt_etc['ctt'][2:]

ls_ctt, ctt_theo=list(ls_ctt), list(ctt_theo)

Dtt_theo=[ls_ctt[i]*(ls_ctt[i]+1)*ctt_theo[i]/(2*np.pi) for i in range(len(ls_ctt))]

#dict={n_band:{col_name:[values], ...}, ...}
band_dict={}
arr_names=['ls', 'Dtt', 'err_Dtt', 'Cgg', 'err_Cgg', 'Ctg', 'err_Ctg']
theo_dict={}

band=1
while band<=4:
    dict_band=pd.read_csv('errors_data_wmap9QVW_xsc{}_jeffrey_dipfix_50000_ttggtg_new.dat'.format(band), header=None, names=arr_names, sep=' ')
    band_dict[band]=dict_band

    band_corrs=pd.read_csv('../../theoretical/tables/ctg_cgg_band{}.dat'.format(band), sep=' ', header=None, index_col=0, names=['ls', 'ctg', 'cgg'])

    print('len(ctg_theo)=', len(band_corrs['ctg']))
    theo_dict[band]=band_corrs

    band+=1

print(theo_dict)

import matplotlib

matplotlib.rcParams.update({'font.size': 15})

fig, axs = plt.subplots(4,3, figsize=(20, 16))

band=1
while band<=4:
    print("Starting band {}".format(band))
    #Dtt plot for each band:
    axs[band-1,0].errorbar(band_dict[band]['ls'], band_dict[band]['Dtt'], yerr=band_dict[band]['err_Dtt'], label="Best-fit+total error", fmt='ko')
    axs[band-1,0].plot(ls_ctt, Dtt_theo, label=r'$\Lambda$CDM fiducial spectrum')
    axs[band-1,0].set_ylabel(r'$(\ell+1)\ell/2\pi C_\ell^{tt}$ $[\mu K^2]$', fontsize=22)
    if band==4:#xlabel only in bottom plots
        axs[band-1,0].set_xlabel(r'$\ell$', fontsize=25)
    axs[band-1,0].set_xscale('log')
    axs[band-1,0].set_xlim(1.8, 90)
    axs[band-1,0].set_ylim(0, 4000)
    if band==1:#legend only in the top-left plot 
        axs[band-1,0].legend()
    print("Ctt finished")

    #Cgg plot for each band:
    axs[band-1,1].errorbar(band_dict[band]['ls'], band_dict[band]['Cgg'], yerr=band_dict[band]['err_Cgg'], fmt='ko')
    axs[band-1,1].plot(theo_dict[band]['ls'], theo_dict[band]['cgg']) 
    if band==4:
        axs[band-1,1].set_xlabel(r'$\ell$', fontsize=25)
    axs[band-1,1].set_ylabel(r'$C_\ell^{gg}$', fontsize=22)
    axs[band-1,1].set_xlim(1.8, 90)
    axs[band-1,1].set_ylim(5e-5, 3e-2)
    axs[band-1,1].set_xscale('log')
    axs[band-1,1].set_yscale('log')
    print('Cgg finished')

    #Ctg plot for each band:
    axs[band-1,2].errorbar(band_dict[band]['ls'], band_dict[band]['Ctg'], yerr=band_dict[band]['err_Ctg'], fmt='ko')
    axs[band-1,2].plot(theo_dict[band]['ls'], theo_dict[band]['ctg']) 
    if band==4:
        axs[band-1,2].set_xlabel(r'$\ell$', fontsize=25)
    axs[band-1,2].set_ylabel(r'$C_\ell^{tg}$ $[\mu K]$', fontsize=22)
    axs[band-1,2].set_xlim(1.8, 90)
    axs[band-1,2].set_ylim(-0.8,0.8)
    axs[band-1,2].set_xscale('log')
    axs[band-1,2].text(32, 0.6, 'Band {}'.format(band), fontsize=26)
    print("Band {} finished, moving to next".format(band))
    band+=1

fig.tight_layout()
plt.savefig("Full_Data_Plot.png")

