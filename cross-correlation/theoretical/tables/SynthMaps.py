from matplotlib import pyplot as plt
import healpy as hp
import numpy as np
from tqdm import tqdm
import pandas as pd

def MultiSynthesize(N, direct='alms', cmb_sample=True, galaxy_sample=True):

#Synthesizes N maps from the .fits files in the 'direct' directory (must be in the same directory as the running script). If X_sample=True, saves a figure of one of the synthsized maps.

	ctts, cggs, ctgs=[],[],[]
	#list of lists. Each "sublist" is the corresponding list of Cltg.

	for i in tqdm(range(N)):
	    cmb_alm=hp.read_alm('{0}/cmb_alm{1}.fits'.format(direct,i+1))
	    galaxy_alm=hp.read_alm('{0}/galaxy_alm{1}.fits'.format(direct,i+1))
	    
	    cmb_alm=cmb_alm.astype('complex')
	    galaxy_alm=galaxy_alm.astype('complex')

	    mapCMB=hp.alm2map(cmb_alm, nside=64, lmax=128)
	    mapgalaxy=hp.alm2map(galaxy_alm, nside=64, lmax=128)
	    
	    ctt=hp.anafast(mapCMB, lmax=128)
	    cgg=hp.anafast(mapgalaxy, lmax=128)
	    ctg=hp.anafast(mapCMB, map2=mapgalaxy, lmax=128)
	    
	    ctts.append(ctt)
	    cggs.append(cgg)
	    ctgs.append(ctg)
	    
	if cmb_sample==True:
		hp.mollview(mapCMB, title="Um dos mapas sintetizados (CMB)", unit=r'$\mu K$', cmap='RdYlBu_r')
		hp.graticule()
		plt.savefig('CMB_sintetizado.png')
		
	if galaxy_sample==True:
		hp.mollview(mapgalaxy, title='Um dos mapas sintetizados (galáxias)', cmap='RdYlBu_r')
		hp.graticule()
		plt.savefig('galaxias_sintetizado.png')
	    
	return ctts, cggs, ctgs

def Avg_Plots(ctts, cggs, ctgs, plot_ctt=False, plot_cgg=False, plot_ctg=True, plot_theory=True):

#Plots the average of the selected spectra. The Ctg plot is always binned, in the future I might work on giving the user the option to choose. If plot_theory=True, plots the theoretical curves with the data. PS.: Binning not implemented yet.

	theoretical_spectra=pd.read_csv('final_table.dat', sep=' ', header=None, names=['l', 'ctt', 'cgg', 'ctg'])
	ls_th, ctt_th=theoretical_spectra['l'], theoretical_spectra['ctt']
	cgg_th, ctg_th=theoretical_spectra['cgg'], theoretical_spectra['ctg']

	avg_ctt, avg_cgg,avg_ctg=[],[],[] #list of average(Cltg) for each fixed l (0<=l<=128)

	for i in range(len(ctgs[0])):
	    numerator=sum([ctgs[j][i] for j in range(len(ctgs))])
	    avgi=numerator/len(ctgs)
	    avg_ctg.append(avgi)

	for i in range(len(ctts[0])):
	    numerator=sum([ctts[j][i] for j in range(len(ctts))])
	    avgi=numerator/len(ctts)
	    avg_ctt.append(avgi)

	for i in range(len(cggs[0])):
	    numerator=sum([cggs[j][i] for j in range(len(cggs))])
	    avgi=numerator/len(cggs)
	    avg_cgg.append(avgi)

	#print(len(avg_ctg))
	#ls_bins, ctg_binned=log_bin_data(avg_ctg)
	#Dtg=[ls_bins[i]*(ls_bins[i]+1)*ctg_binned[i]/(2*np.pi) for i in range(len(ls_bins))]

	ls_tg=np.arange(len(avg_ctg))
	Dtg=[ls_tg[i]*(ls_tg[i]+1)*avg_ctg[i]/(2*np.pi) for i in range(len(ls_tg))]
	Dtg_th=[ls_th[i]*(ls_th[i]+1)*ctg_th[i]/(2*np.pi) for i in range(len(ls_th))]

	ls_tt=np.arange(len(avg_ctt))
	Dtt=[ls_tt[i]*(ls_tt[i]+1)*avg_ctt[i]/(2*np.pi) for i in range(len(ls_tt))]
	Dtt_th=[ls_th[i]*(ls_th[i]+1)*ctt_th[i]/(2*np.pi) for i in range(len(ls_th))]

	ls_gg=np.arange(len(avg_cgg))
	
	if plot_ctg:
		plt.figure('plotCtg')
		plt.title('Média das correlações cruzadas de {} mapas.'.format(len(ctgs)))
		plt.xlabel(r'$\ell$')
		plt.ylabel(r'$\ell(\ell+1)C_{\ell}^{tg}/2\pi$')
		plt.xscale('log')
		
		if plot_theory:
			plt.plot(ls_th[2:], Dtg_th[2:], label='Espectro teórico')
			
		plt.plot(ls_tg[2:], Dtg[2:], label='Mapas sintetizados')
		plt.legend()
		plt.savefig('media_plot_ctg_for_N{}.png'.format(len(ctgs)))
		
	if plot_ctt:
		plt.figure('plotCtt')
		plt.title('Média das autocorrelações de {} mapas'.format(len(ctts)))
		plt.xlabel(r'$\ell$')
		plt.ylabel(r'$\ell(\ell+1)C_{\ell}^{tt}/2\pi$')
		plt.xscale('log')
		
		if plot_theory:
			plt.plot(ls_th[2:], Dtt_th[2:], label='Espectro teórico')
			
		plt.plot(ls_tt[2:], Dtt[2:], label='Mapas sintetizados')
		plt.legend()
		plt.savefig('media_plot_ctt_for_N{}.png'.format(len(ctts)))
		
	if plot_cgg:
		plt.figure('plotCgg')
		plt.title('Média das autocorrelações de {} mapas'.format(len(cggs)))
		plt.xlabel(r'$\ell$')
		plt.ylabel(r'$C_{\ell}^{gg}$')
		plt.xscale('log')
		plt.yscale('log')
		
		if plot_theory:
			plt.plot(ls_th[2:], cgg_th[2:], label='Espectro teórico')
			
		plt.plot(ls_gg[2:], avg_cgg[2:], label='Mapas sintetizados')
		plt.legend()
		plt.savefig('media_plot_cgg_for_N{}.png'.format(len(cggs)))
		
def main():
	N=int(input('Quantos mapas a serem sintetizados? (máximo=1e4)\n'))
	ctts, cggs, ctgs=MultiSynthesize(N, cmb_sample=False, galaxy_sample=False)
	Avg_Plots(ctts,cggs,ctgs,plot_ctt=True, plot_cgg=True)
	
if __name__=='__main__':
	main()
	
