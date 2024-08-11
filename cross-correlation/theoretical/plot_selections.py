import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import scipy.special as scp

band_pars={1:{'z0':0.043,'beta':1.825,'lambda':1.524, 'bg50':1.32, 'bg96':1.37},
            2:{'z0':0.054,'beta':1.800,'lambda':1.600, 'bg50':1.34, 'bg96':1.35},
            3:{'z0':0.067,'beta':1.765,'lambda':1.636, 'bg50':1.29, 'bg96':1.34},
            4:{'z0':0.084,'beta':1.723,'lambda':1.684, 'bg50':1.28, 'bg96':1.29},
            'min':{'z0':0.1508,'beta':3.088,'lambda':4.9401, 'bg50':1, 'bg96':1}}

def dNdz(z, z0, lbda, beta):
    k=(beta/scp.gamma(lbda))*(z/z0)**(beta*lbda-1)
    exp=np.exp(-(z/z0)**beta)/z0

    return k*exp

zmin, zmax=0, 0.5
npoints=1000
dz=(zmax-zmin)/npoints

zs=[zmin+i*dz for i in range(npoints)]

fs_2mass={1:[], 2:[], 3:[], 4:[]}

for band in [1,2,3,4]:
    z0, beta, lbda=band_pars[band]['z0'],band_pars[band]['beta'],band_pars[band]['lambda']
    for z in zs:
        fz=dNdz(z,z0,beta,lbda)
        fs_2mass[band].append(fz)

z0_min, beta_min, lbda_min=band_pars['min']['z0'],band_pars['min']['beta'],band_pars['min']['lambda']
fmin=[dNdz(z,z0_min,lbda_min,beta_min) for z in zs]

matplotlib.rcParams.update({'font.size':21})

plt.figure(figsize=(9,6))

for band in [1,2,3,4]:
    plt.plot(zs, fs_2mass[band], label="Band {}".format(band))

plt.plot(zs, fmin, label="Optimal band")
plt.xlabel('z')
plt.ylabel(r'$\frac{dN}{dz}$')
plt.legend(fontsize=18)
plt.tight_layout()
plt.savefig('Compare_Selections.png')
#plt.show()
