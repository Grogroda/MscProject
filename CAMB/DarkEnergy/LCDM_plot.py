import pandas as pd
from matplotlib import pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable

data_full=pd.read_csv("COM_PowerSpect_CMB-TT-full_R3.01.csv", header=[0])
data_binned=pd.read_csv('COM_PowerSpect_CMB-TT-binned_R301.csv', header=[0])

simswa0=pd.read_csv("ClTT_for_Wa0.00.csv", header=[0])

plt.figure()
#plt.title(r'Espectro de Potências do CMB para o modelo $\Lambda$CDM', loc='left')
plt.xlabel(r'$\ell$')
	
axMain=plt.subplot(111)
axMain.plot(simswa0['l'][30:], simswa0['ClTT001'][30:])
axMain.errorbar(data_binned['l'], data_binned['Dl'], yerr=[data_binned['-dDl'], data_binned['+dDl']], fmt='.', color='tab:red')
axMain.set_xscale('linear')
axMain.set_xlim((simswa0['l'][30], 1900))
axMain.spines['left'].set_visible(False)
axMain.yaxis.set_ticks_position('right')
axMain.yaxis.set_visible(False)
	
divider=make_axes_locatable(axMain)
axLog = divider.append_axes("left", size=2, pad=0, sharey=axMain)
axLog.plot(simswa0['l'][2:31], simswa0['ClTT001'][2:31], label='Simulação (CAMB)')
axLog.errorbar(data_full['l'][:29], data_full['Dl'][:29], yerr=[data_full['-dDl'][:29], data_full['+dDl'][:29]], fmt='.', color='tab:red', label='Dados (Planck)')
axLog.set_xscale('log')
axLog.set_xlim((1.5, simswa0['l'][30]))
axLog.spines['right'].set(ls=':')
	
plt.ylabel(r'$\ell(\ell+1)C_\ell^{TT}/2\pi \; (\mu K)^2$')

plt.setp(axLog.get_xticklabels(), visible=True)

plt.legend()
plt.savefig('full_doublscale_LCDM.png')
