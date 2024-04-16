import numpy as np
from matplotlib import pyplot as plt
import argparse
import pandas as pd

l=2

parser=argparse.ArgumentParser()
parser.add_argument('-l', '--multipole')
args=parser.parse_args()

if args.multipole!=None:
    l=args.multipole
print("l=", l)

import matplotlib
matplotlib.rcParams.update({'font.size':17})

hist=pd.read_csv('hists/histogram_l{}.txt'.format(l), header=None, sep=' ')
bins, weights=hist[0], hist[1][:-1]
print(bins)
print(weights)

plt.figure(figsize=(10,8))
plt.hist(bins[:-1], bins, weights=weights)
plt.title(r'$\ell={}$'.format(l))
plt.xlabel(r'$C^{tg}$')
plt.ylabel(r'$P(C^{tg})$')
plt.savefig('hist_l{}.png'.format(l))
