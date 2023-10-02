import numpy as np
import pandas as pd
from matplotlib import pyplot as plt

#Reminder:
#0=Ctt
#1=Cgg
#2=Ctg

best_fit=pd.read_csv('max_like.csv', names=['#', 'l', 'Cl', '-log(Lmax)'])
camb_path='../../../../CAMB/LCDM/'

ctt_theo=pd.read_csv('{}cttLCDM_camb.dat'.format(camb_path), sep=' ', header=None, names=['l', 'Cl'])
cgg_theo=pd.read_csv('{}cgg2MASSband1_camb.dat'.format(camb_path), sep=' ', header=None, names=['l', 'Cl'])
ctg_theo=pd.read_csv('{}ctg2MASSband1_camb.dat'.format(camb_path), sep=' ', header=None, names=['l', 'Cl'])

#print('ctt=\n', ctt_theo)
#print('cgg=\n', cgg_theo)
#print('ctg=\n', ctg_theo)

ls=[]
ctt=[]
cgg=[]
ctg=[]
for i in range(len(best_fit["#"])):
    if best_fit['#'][i]==0:
       ls.append(best_fit['l'][i])
       ctt.append(best_fit["Cl"][i])
    elif best_fit["#"][i]==1:
        cgg.append(best_fit["Cl"][i])
    elif best_fit["#"][i]==2:
        ctg.append(best_fit["Cl"][i])

#print('ls=', ls)
#print('ctg_ls=\n', ctg_theo['l'])

kctt=[]
kctt_theo=[]
for i in range(len(ctt)):
    l=ls[i]
    kcl=l*(l+1)/(2*np.pi)*ctt[i]
    kcl_theo=l*(l+1)/(2*np.pi)*ctt_theo['Cl'][i]
    kctt.append(kcl)
    kctt_theo.append(kcl_theo)

######## Begin figures ########

fig_tt, axtt=plt.subplots()
fig_gg, axgg=plt.subplots()
fig_tg, axtg=plt.subplots()

axes=[axtt, axgg, axtg]

axtt.plot(ls, kctt_theo, label='CAMB')
axtt.scatter(ls, kctt, s=5, label='Best-fit of samples', c='tab:orange')
axtt.set_ylabel(r'$\ell(\ell+1)C_l^{tt}/2\pi$')

axgg.plot(ls, cgg_theo['Cl'][2:], label='CAMB')
axgg.scatter(ls, cgg, s=5, label='Best-fit of samples', c='tab:orange')
axgg.set_ylabel(r'$C_l^{gg}$')
axgg.set_yscale('log')

axtg.plot(ls, ctg_theo['Cl'][2:], label='CAMB')
axtg.scatter(ls, ctg, s=5, label='Best-fit of samples', c='tab:orange')
axtg.set_xlabel(r'$C_l^{tg}$')

for ax in axes:
    ax.set_xlabel(r'$\ell$')
    ax.set_xscale('log')
    ax.legend()

fig_tt.savefig('cltt.png')
fig_gg.savefig('clgg.png')
fig_tg.savefig('cltg.png')
