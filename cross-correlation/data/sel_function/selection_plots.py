import numpy as np
from matplotlib import pyplot as plt
import scipy.special as scp

#reference values correspond to band 2 of the 2MASS catalog:
z0,lbda,beta=0.054,1.6,1.8
z0s=[0.02, 0.06, 0.1, 0.14, 0.18]
betas=[0.4, 1.0, 1.6, 2.2, 2.8]
lbdas=[0.4, 1.0, 1.6, 2.2, 2.8]

def selection(z,z0,lbda,beta):
    
    k=(beta/scp.gamma(lbda))*(z/z0)**(beta*lbda-1)
    exp=np.exp(-(z/z0)**beta)/z0

    return k*exp

#fixing lbda and beta:

def vary_z0(zs):
    #z0s=[0.02*(i+1) for i in range(10)]
    dNdz_list=[]
    for z0 in z0s:
        dNdz=[selection(z,z0,lbda,beta) for z in zs]
        dNdz_list.append(dNdz)

    return dNdz_list, z0s
    
def vary_lbda(zs):
    #lbdas=[0.4*(i+1) for i in range(10)]
    dNdz_list=[]
    for lbda in lbdas:
        dNdz=[selection(z,z0,lbda,beta) for z in zs]
        dNdz_list.append(dNdz)

    return dNdz_list, lbdas
    
def vary_beta(zs):
    #betas=[0.4*(i+1) for i in range(10)]
    dNdz_list=[]
    for beta in betas:
        dNdz=[selection(z,z0,lbda,beta) for z in zs]
        dNdz_list.append(dNdz)

    return dNdz_list, betas

def main():

    zs=[0.0005*(i+1) for i in range(10000)]
    select4z0s, z0s=vary_z0(zs)
    plt.figure()
    for i in range(len(select4z0s)):
        plt.plot(zs, select4z0s[i], label=r'$z_0$={:.2f}'.format(z0s[i]))
    
    plt.xlabel('z')
    plt.ylabel(r'$\frac{dN}{dz}$')
    plt.xlim((0,1))
    plt.title(r'Fixed parameters: $\lambda={0}$ and $\beta={1}$'.format(lbda,beta))
    plt.legend()
    plt.savefig('z0_vary.png')
    
    select4lbdas, lbdas=vary_lbda(zs)
    plt.figure()
    for i in range(len(select4lbdas)):
        plt.plot(zs, select4lbdas[i], label=r'$\lambda$={:.1f}'.format(lbdas[i]))
    
    plt.xlabel('z')
    plt.ylabel(r'$\frac{dN}{dz}$')
    plt.title(r'Fixed parameters: $z_0={0}$ and $\beta={1}$'.format(z0,beta))
    plt.legend()
    plt.xlim((0,0.4))
    plt.savefig('lbda_vary.png')
    
    select4betas, betas=vary_beta(zs)
    plt.figure()
    for i in range(len(select4z0s)):
        plt.plot(zs, select4betas[i], label=r'$\beta$={:.1f}'.format(betas[i]))
    
    plt.xlabel('z')
    plt.ylabel(r'$\frac{dN}{dz}$')
    plt.title(r'Fixed parameters: $\lambda={0}$ and $z_0={1}$'.format(lbda,z0))
    plt.legend()
    plt.xlim((0,0.6))
    plt.savefig('beta_vary.png')   

    #Triple plot for dissertation:
    
    import matplotlib
    matplotlib.rcParams.update({'font.size':15})

    fig = plt.figure(figsize=[15, 5])
    gs = fig.add_gridspec(1,3,wspace=0)
    axs = gs.subplots(sharey=True)
    axs[0].set_ylabel(r"$\frac{dN}{dz}$")
    axs[0].set_xlabel(r'$z$') #z0 vary
    axs[1].set_xlabel(r'$z$') #beta vary
    axs[2].set_xlabel(r'$z$') #lbda vary

    for i in range(len(select4z0s)):
        axs[0].plot(zs, select4z0s[i], label=r'$z_0$={:.2f}'.format(z0s[i]))
        axs[1].plot(zs, select4betas[i], label=r'$\beta$={:.1f}'.format(betas[i]))
        axs[2].plot(zs, select4lbdas[i], label=r'$\lambda$={:.1f}'.format(lbdas[i]))

    axs[0].set_ylim((-1, 40))
    axs[0].set_xlim((0, 0.67))
    axs[1].set_xlim((0, 0.47))
    axs[2].set_xlim((0,0.4))
    axs[0].legend()
    axs[1].legend()
    axs[2].legend()

    plt.tight_layout()
    plt.savefig('Selection_TriplePlot.png')

if __name__=='__main__':
    main()
