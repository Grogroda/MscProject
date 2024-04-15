import numpy as np
import matplotlib.pyplot as plt
import scipy.special as scp
from tqdm import tqdm

band_pars={1:{'z0':0.043, 'beta':1.825, 'lbda':1.524},
           2:{'z0':0.054, 'beta':1.800, 'lbda':1.600},
           3:{'z0':0.067, 'beta':1.765, 'lbda':1.636},
           4:{'z0':0.084, 'beta':1.723, 'lbda':1.684}
           }

def selection(z, z0, lbda, beta):

    k=(beta/scp.gamma(lbda))*(z/z0)**(beta*lbda-1)
    exp=np.exp(-(z/z0)**beta)/z0

    return k*exp

npoints=10000
zmin, zmax=0, 0.5
dz=(zmax-zmin)/npoints

bands={1:{'zs':[], 'f':[]}, 2:{'zs':[], 'f':[]}, 3:{'zs':[], 'f':[]}, 4:{'zs':[], 'f':[]}}

for band in range(1,5):
    z=zmin
    z0=band_pars[band]['z0']
    beta=band_pars[band]['beta']
    lbda=band_pars[band]['lbda']

    for i in tqdm(range(npoints), desc="Band {}".format(band)):
        bands[band]['zs'].append(z)
        f=selection(z, z0, beta, lbda)
        bands[band]['f'].append(f)

        z+=dz

import matplotlib
matplotlib.rcParams.update({'font.size':15})

plt.figure(figsize=[8,5])

for band in range(1,5):
    plt.plot(bands[band]['zs'], bands[band]['f'], label='Band {}'.format(band))

plt.xlim((0,0.3))
plt.ylim((0,20))
plt.xlabel(r'$z$')
plt.ylabel(r'$\frac{dN}{dz}$')
plt.gca().set_aspect(1/80)
plt.legend()
plt.savefig('selection_2MASS.png')
