import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import tabular as tb
#import HistSearch as hs

lbda=float(input('Value of lbda: '))

def gauss(x,mean,var): #returns f(x) where f is a normalized gaussian distribution with mean=mean and variance=var

	return 1/(np.sqrt(2*np.pi*var))*np.exp(-(x-mean)**2/(2*var))

ls=[i for i in range(129)]

z0min, z0max, dz0=0.02, 0.28, 0.01
Nz0=round((z0max-z0min)/dz0)+1

betamin, betamax, dbeta=0.6, 4.6, 0.1
Nbeta=round((betamax-betamin)/dbeta)+1

z0s, betas=[z0min+dz0*i for i in range(Nz0)],[betamin+dbeta*j for j in range(Nbeta)]
#print('z0s=', z0s)
#print('betas=', betas)
pdfs=pd.read_csv('../tables/hists/gauss_fits.txt', sep=' ', header=None, names=['l', 'mean', 'variance'])
#print(pdfs) 
logPs=[]
error_points=[]

for beta in betas: #beta vai ser o eixo Y da matriz de cores
    sublist=[] #cria uma nova linha para a matriz de cores
    for z0 in z0s: #z0 vai ser o eixo X da matriz de cores
        logP=0
        Cls=pd.read_csv('../ctg_files/funcs_lbda_{0:.1f}/beta_{2:.1f}.z0_{1:.2f}.dat'.format(lbda,z0, beta), sep=' ', header=None, names=['l', 'Cltg'])
        #print(Cls.to_string())
        #print('len(cltg)=', len(Cls['Cltg']))
        #print('len(ls)=', len(Cls['l']))
    	
        if len(Cls['Cltg'])==129:
            for l in ls[2:]:
                Ctg=Cls['Cltg'][l]
                logp=np.log(gauss(Ctg, pdfs['mean'][l], pdfs['variance'][l]))
                logP=logP+logp
            
        else:
            print("Error in (z0,beta)=({0},{1})".format(z0,beta))
            print('len(l)=', len(Cls['l']))
            logP=1#max(max(logPs))
            error_points.append((z0,beta))

        sublist.append(logP) #vai preenchendo as colunas da matriz de cores

    logPs.append(sublist)

MASS_df=pd.read_csv('../2MASS/band1_bg1.37.dat', sep=' ', header=None, names=['l', 'Cltg'])
ls, ctg=MASS_df['l'], MASS_df['Cltg']
logP_2MASS=0
for l in ls[2:]:
    logp=np.log(gauss(ctg[l], pdfs['mean'][l], pdfs['variance'][l]))
    logP_2MASS+=logp

#logmax=max(max(logPs))
#print(logmax)
Ps=[]

for i in range(len(betas)):
	Ps.append([])
	for j in range(len(z0s)):
		P=np.exp(logPs[i][j]-logP_2MASS)
		Ps[i].append(P)
		
		if P==min(map(min,Ps)):
			minimum=[z0s[j], betas[i]]

#print(logPs)
#print(lbdas)

#imin=Ps.index(min(Ps))
#jmin=Ps[imin].index(min(min(Ps)))

print('Error Points: ', error_points)
print('Minimum Cell (z0,beta)=({0},{1})'.format(minimum[0], minimum[1]))
#print('Ps=', Ps)
#print('Vmin=', min(map(min, Ps)))


tb.tab_export(z0s, betas, Ps, 'map_files/cmap_lbda_{0:.1f}.txt'.format(lbda))

plt.figure()
plt.title(r'Color map of $P=\prod_{\ell=2}^{128} f(C_\ell^{tg})$, where $\lambda=$'+str(lbda))
plt.xlabel(r'$z_0$')
plt.ylabel(r'$\beta$')
plt.pcolormesh(z0s, betas, Ps, shading='nearest', vmin=0.97, vmax=1.04, cmap='RdYlBu_r')
plt.colorbar(label=r'$P/P_1$')
plt.scatter([point[0] for point in error_points], [point[1] for point in error_points], c='black', marker='x')
plt.savefig('plots/lbda{0:.1f}.png'.format(lbda))

