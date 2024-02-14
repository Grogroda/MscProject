import numpy as np
import pandas as pd
from matplotlib import pyplot as plt

#Reminder:
#0=Ctt
#1=Cgg
#2=Ctg

best_fit=pd.read_csv('max_likeQ.dat', sep=' ', names=['#', 'l', 'Cl', 'min(-log(L))'])
camb_path='../../../../CAMB/LCDM/'

print(best_fit)

ctt_theo=pd.read_csv('{}cttLCDM_camb.dat'.format(camb_path), sep=' ', header=None, names=['l', 'Cl'])
cgg_theo=pd.read_csv('{}cgg2MASSband1_camb.dat'.format(camb_path), sep=' ', header=None, names=['l', 'Cl'])
ctg_theo=pd.read_csv('{}ctg2MASSband1_camb.dat'.format(camb_path), sep=' ', header=None, names=['l', 'Cl'])

# reminder: ctt output on CAMB is usually l(l+1)/2*pi Cltt, the same goes for the best-fit output of powertools

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

######## Begin figures ########

#Ideia: Criar uma opção que o usuário pode pegar 3 figuras separadas ou uma figura junta
fig, axes=plt.subplots(1,3, figsize=(20,5))
'''
fig_tt, axtt=plt.subplots()
fig_gg, axgg=plt.subplots()
fig_tg, axtg=plt.subplots()

axes=[axtt, axgg, axtg]
'''

lmin, lfinal=2, 50 #samples lists start at l=lmin. Starting from l=lmax, the values are fixed, not sampled.
final_ind=lfinal-lmin #index of the last sampled value of the list
lmax=96 #lmax for theoretical data

axes[0].plot(ls, ctt_theo['Cl'][2:lmax+1], label='CAMB')
axes[0].scatter(ls[:final_ind], ctt[:final_ind], s=5, label='Best-fit of samples', c='tab:orange')
axes[0].set_ylabel(r'$\ell(\ell+1)C_l^{tt}/2\pi$', fontsize=14)

axes[1].plot(ls, cgg_theo['Cl'][2:], label='CAMB')
axes[1].scatter(ls[:final_ind], cgg[:final_ind], s=5, label='Best-fit of samples', c='tab:orange')
axes[1].set_ylabel(r'$C_l^{gg}$', fontsize=14)
axes[1].set_yscale('log')

axes[2].plot(ls, ctg_theo['Cl'][2:], label='CAMB')
axes[2].scatter(ls[:final_ind], ctg[:final_ind], s=5, label='Best-fit of samples', c='tab:orange')
axes[2].set_ylabel(r'$C_l^{tg}$', fontsize=14)

for ax in axes:
    ax.set_xlabel(r'$\ell$', fontsize=17)
    ax.set_xscale('log')
    ax.legend(fontsize=14)

fig.savefig('cl_tripleQ.png', bbox_inches='tight')
'''
fig_gg.savefig('clgg.png')
fig_tg.savefig('cltg.png')
'''
