import numpy as np
import matplotlib.pyplot as plt
import scipy.special as scp

def selection(z,z0,lbda,beta):
    
    k=(beta/scp.gamma(lbda))*(z/z0)**(beta*lbda-1)
    exp=np.exp(-(z/z0)**beta)/z0

    return k*exp
    
dNdz_list=[]
zs=[0.0005*(i+1) for i in range(10000)]
z0, lbda, beta=0.16, 2.2, 1.8
dNdz=[selection(z,z0,lbda,beta) for z in zs]
dNdz_list.append(dNdz)
   
plt.figure()  
plt.title(r'$z0={0}$, $\beta={1}$ e $\lambda={2}$'.format(z0, beta, lbda))
plt.plot(zs, dNdz)
plt.xlabel('z')
plt.ylabel(r'$\frac{dN}{dz}$')
plt.xlim((0,1))
plt.savefig('test1.png')
     
