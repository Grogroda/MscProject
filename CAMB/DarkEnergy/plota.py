'''
Esse código será usado apenas para fazer os plots das simulações e dos dados
'''

import pandas as pd
from matplotlib import pyplot as plt
import matplotlib
from mpl_toolkits.axes_grid1 import make_axes_locatable

matplotlib.rcParams.update({'font.size':15, 'figure.figsize':[10, 7.5]})

data_full=pd.read_csv("COM_PowerSpect_CMB-TT-full_R3.01.csv", header=[0])
data_binned=pd.read_csv('COM_PowerSpect_CMB-TT-binned_R301.csv', header=[0])

simwa0=pd.read_csv("ClTT_for_Wa0.00.csv", header=[0])
simwa033=pd.read_csv("ClTT_for_Wa0.33.csv", header=[0])
simwa067=pd.read_csv("ClTT_for_Wa0.67.csv", header=[0])
simwa1=pd.read_csv("ClTT_for_Wa1.00.csv", header=[0])

simswa=[simwa0, simwa033, simwa067, simwa1]

for i in range(4):

	plt.figure(i)
	plt.xlabel(r'$\ell$',fontsize=20)
	plt.ylim((0, 6200))
	
	axMain=plt.subplot(111)
	axMain.errorbar(data_binned['l'], data_binned['Dl'], yerr=[data_binned['-dDl'], data_binned['+dDl']], fmt='k.', label='Planck Data')
	axMain.plot(simswa[i]['l'][30:], simswa[i]['ClTT001'][30:], label=r'$c_s^2=0.01$')
	axMain.plot(simswa[i]['l'][30:], simswa[i]['ClTT01'][30:], label=r'$c_s^2=0.1$')
	axMain.plot(simswa[i]['l'][30:], simswa[i]['ClTT1'][30:], label=r'$c_s^2=1$')
	axMain.set_xscale('linear')
	axMain.set_xlim((simswa[i]['l'][30], 1750))
	axMain.spines['left'].set_visible(False)
	axMain.yaxis.set_ticks_position('right')
	axMain.yaxis.set_visible(False)
	plt.legend()
	
	divider=make_axes_locatable(axMain)
	axLog = divider.append_axes("left", size=2, pad=0, sharey=axMain)
	axLog.errorbar(data_full['l'][:29], data_full['Dl'][:29], yerr=[data_full['-dDl'][:29], data_full['+dDl'][:29]], fmt='k.', label='Dados para baixo l')
	axLog.plot(simswa[i]['l'][2:31], simswa[i]['ClTT001'][2:31], label=r'$c_s^2=0.01$')
	axLog.plot(simswa[i]['l'][2:31], simswa[i]['ClTT01'][2:31], label=r'$c_s^2=0.1$')
	axLog.plot(simswa[i]['l'][2:31], simswa[i]['ClTT1'][2:31], label=r'$c_s^2=1$')
	axLog.set_xscale('log')
	axLog.set_xlim((1.5, simswa[i]['l'][30]))
	axLog.spines['right'].set(ls=':')
	plt.text(2, 5800, r'$w_a={:.2f}$'.format(i*1/3))
	
	plt.ylabel(r'$D_\ell^{TT}$ $[\mu K^2]$')
	
	plt.setp(axLog.get_xticklabels(), visible=True)
	
	plt.savefig('full_doublscale_Wafixo{:.2f}.png'.format(i*1/3))

simcs001=pd.read_csv('ClTT_for_cs20.01.csv')
simcs01=pd.read_csv('ClTT_for_cs20.10.csv')
simcs1=pd.read_csv('ClTT_for_cs21.00.csv')

simscs=[simcs001, simcs01, simcs1]
	
for i in range(3):

	plt.figure(5+i)
	plt.xlabel(r'$\ell$', fontsize=20)
	plt.ylim((0, 6200))
	
	axMain=plt.subplot(111)
	axMain.errorbar(data_binned['l'], data_binned['Dl'], yerr=[data_binned['-dDl'], data_binned['+dDl']], fmt='k.', label='Planck Data')
	axMain.plot(simscs[i]['l'][30:], simscs[i]['ClTT0'][30:], label=r'$w_a=0$')
	axMain.plot(simscs[i]['l'][30:], simscs[i]['ClTT033'][30:], label=r'$w_a=0.33$')
	axMain.plot(simscs[i]['l'][30:], simscs[i]['ClTT067'][30:], label=r'$w_a=0.67$')
	axMain.plot(simscs[i]['l'][30:], simscs[i]['ClTT1'][30:], label=r'$w_a=1$')
	axMain.set_xscale('linear')
	axMain.set_xlim((simscs[i]['l'][30], 1750))
	axMain.spines['left'].set_visible(False)
	axMain.yaxis.set_ticks_position('right')
	axMain.yaxis.set_visible(False)
	plt.legend()
	
	divider=make_axes_locatable(axMain)
	axLog = divider.append_axes("left", size=2, pad=0, sharey=axMain)
	axLog.errorbar(data_full['l'][:29], data_full['Dl'][:29], yerr=[data_full['-dDl'][:29], data_full['+dDl'][:29]], fmt='k.', label='Dados para baixo l')
	axLog.plot(simscs[i]['l'][2:31], simscs[i]['ClTT0'][2:31], label=r'$w_a=0$')
	axLog.plot(simscs[i]['l'][2:31], simscs[i]['ClTT033'][2:31], label=r'$w_a=0.33$')
	axLog.plot(simscs[i]['l'][2:31], simscs[i]['ClTT067'][2:31], label=r'$w_a=0.67$')
	axLog.plot(simscs[i]['l'][2:31], simscs[i]['ClTT1'][2:31], label=r'$w_a=1$')
	axLog.set_xscale('log')
	axLog.set_xlim((1.5, simscs[i]['l'][30]))
	axLog.spines['right'].set(ls=':')
	plt.text(2, 5800, r'$c_s^2={:.2f}$'.format(10**(i-2)))
	
	plt.ylabel(r'$D_\ell^{TT}$ $[\mu K^2]$')
	
	plt.setp(axLog.get_xticklabels(), visible=True)
	
	plt.savefig('full_doublscale_Cs2fixo{:.2f}.png'.format(10**(i-2)))

