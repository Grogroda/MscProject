import numpy as np
from matplotlib import pyplot as plt
import argparse
import pandas as pd

def gauss(x, mu, var):
    k=1/(np.sqrt(2*np.pi*var))
    exp=np.exp(-1/(2*var)*(x-mu)**2)

    return k*exp

l=2

parser=argparse.ArgumentParser()
parser.add_argument('-l', '--multipole')
args=parser.parse_args()

if args.multipole!=None:
    l=int(args.multipole)
print("l=", l)

import matplotlib
matplotlib.rcParams.update({'font.size':20})

hist=pd.read_csv('hists/histogram_l{}.txt'.format(l), header=None, sep=' ')
bins, weights=hist[0], hist[1][:-1]
print(bins)
print(weights)

#Gauss fit points to plot:
gauss_fits=pd.read_csv('hists/gauss_fits.txt', header=None, names=['ls', 'mu', 'var'], sep=' ')

xmin, xmax=min(bins), max(bins)
Npoints=1000
dx=(xmax-xmin)/Npoints
xs, ys=[],[]
x=xmin
for i in range(Npoints):
    xs.append(x)
    mu=gauss_fits['mu'][l]
    var=gauss_fits['var'][l]
    y=gauss(x, mu, var)
    ys.append(y)

    x+=dx

plt.figure(figsize=(10,8))
plt.hist(bins[:-1], bins, density=True,weights=weights, label='Normalized\ndistribution')
plt.plot(xs, ys, label='Gaussian Fit')
#plt.title(r'$\ell={}$'.format(l))
plt.xlabel(r'$C^{tg}_{\ell}$', fontsize=25)
#plt.ylabel(r'$f_{\ell}(C^{tg}_{\ell})$', fontsize=25)
plt.legend(loc='upper left')
plt.tight_layout()
plt.savefig('hist_l{}.png'.format(l))
