import numpy as np
import matplotlib.pyplot as plt
import scipy.special as scp

def selection(z,z0,lbda,beta):
    
    k=(beta/scp.gamma(lbda))*(z/z0)**(beta*lbda-1)
    exp=np.exp(-(z/z0)**beta)/z0

    return k*exp

dNdz_list=[]
zs=[0.0005*(i+1) for i in range(10000)]
z0=float(input('z0='))
lbda=float(input('lambda='))
beta=float(input('beta='))
dNdz=[selection(z,z0,lbda,beta) for z in zs]
dNdz_list.append(dNdz)

plt.figure()
plt.title(r'$z0={0}$, $\beta={1}$ e $\lambda={2}$'.format(z0, beta, lbda))
plt.plot(zs, dNdz)
plt.xlabel('z')
plt.ylabel(r'$\frac{dN}{dz}$')
plt.xlim((0,1))

save=input('Save (sa) or just show (sh)? ')

if save=='sa':
    name=input('Image file name (with image extension, ex: image.png): ')
    plt.savefig(name)

elif save=='sh':
    plt.show()

else:
    print('What do you want from me?')
