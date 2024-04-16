import numpy as np
import healpy as hp
from matplotlib import pyplot as plt
import pandas as pd
from tqdm import tqdm #this module is optional, it shows a progression bar
import SynthMaps

N=int(input('Quantos mapas a serem sintetizados? (máximo=1e4)\n'))

ctts, cggs, ctgs=SynthMaps.MultiSynthesize(N,direct='alms_nocorrelation', cmb_sample=False, galaxy_sample=False) #Synthetizes N maps

ctg_dict={} #dictionary format=l:[all Cltg simulated for that l] 
gauss_df=pd.DataFrame({'l':[0,1], 'mean':[0,0], 'variance':[0,0]}) #this table stores the mean and variance for every l across the synthesized maps

nbins=int(input('Número de bins do histograma:\n'))

for l in tqdm(range(2,len(ctgs[0])), desc='Generating Histograms'): #you can optionally remove the tqdm function to not show a progression bar 
    cls=[ctgs[i][l] for i in range(len(ctgs))]
    ctg_dict[l]=cls
    #hist=np.histogram(cls, bins=nbins)
    #hists_dict[l]=hist
    limits=max(abs(min(cls)), abs(max(cls)))
    
    counts, bin_edges=np.histogram(cls, bins=nbins, range=((-1)*limits, limits), density=True)
    counts=np.append(counts,None) #this last term can't be used, it's here just to have the same length of bin_edges to fit in the dataframe
    hist_df=pd.DataFrame({'bins':bin_edges, 'freqs':counts})
    
    hist_df.to_csv('new_hists/histogram_l{1}.txt'.format(nbins,l), sep=' ', index=False, header=False)

    mean, var=np.mean(cls), np.var(cls)
    gauss_df=gauss_df.append({'l':l, 'mean':mean, 'variance':var}, ignore_index=True)
   # plt.figure()
   # plt.title('N={0}, l={1}'.format(N,l))
   # plt.xlabel(r'$C_l^{tg}$')
   # plt.ylabel('Frequência')
   # plt.hist(bin_edges[:-1], bin_edges, weights=counts[:-1], rwidth=0.9, label='Mean={0:.4f}\nVariance={1:.4f}'.format(mean, var))
   # plt.legend()
   # plt.savefig('histograms/hist_N{0}l{1}.png'.format(N,l))

gauss_df.to_csv('hists/gauss_fits.txt', sep=' ', index=False)
