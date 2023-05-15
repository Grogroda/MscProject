import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import tabular as tb
#import HistSearch as hs

beta=float(input('Value of beta: '))

def gauss(x,mean,var): #returns f(x) where f is a normalized gaussian distribution with mean=mean and variance=var

	return 1/(np.sqrt(2*np.pi*var))*np.exp(-(x-mean)**2/(2*var))

z0min, z0max, dz0=0.01, 0.22, 0.01
Nz0=round((z0max-z0min)/dz0)+1

lbdamin, lbdamax, dlbda=0.8, 8.0, 0.2
Nlbda=round((lbdamax-lbdamin)/dlbda)+1
 
ls=[i for i in range(129)]
z0s, lbdas=[z0min+i*dz0 for i in range(Nz0)],[lbdamin+dlbda*i for i in range(Nlbda)]
pdfs=pd.read_csv('../tables/hists/gauss_fits.txt', sep=' ', header=None, names=['l', 'mean', 'variance'])
#print(pdfs) 
logPs=[]
error_points=[]

for lbda in lbdas: #lbda vai ser o eixo Y da matriz de cores
    sublist=[] #cria uma nova linha para a matriz de cores
    for z0 in z0s: #z0 vai ser o eixo X da matriz de cores
    	logP=0
    	Cls=pd.read_csv('../ctg_files/funcs_beta_{0:.1f}/z0_{1:.2f}.lbda_{2:.1f}.dat'.format(beta,z0, lbda), sep=' ', header=None, names=['l', 'Cltg'])
    	#print(Cls.to_string())
    	#print(len(Cls['Cltg']))

    	if len(Cls['Cltg'])==129:
            for l in ls[2:]:
                #print(l)
                Ctg=Cls['Cltg'][l]
                logp=np.log(gauss(Ctg, pdfs['mean'][l], pdfs['variance'][l]))
                logP=logP+logp
            
    	else:
    	    print("Error in (z0,beta)=({0},{1})".format(z0,beta))
    	    print('len(l)=', len(Cls['l']))
            logP=1
            error_points.append((z0,lbda))

    	sublist.append(logP) #vai preenchendo as colunas da matriz de cores

    logPs.append(sublist)

MASS_df=pd.read_csv('../tables/2MASS/band1_bg1.37.dat', sep=' ', header=None, names=['l', 'Cltg'])
ls, ctg=MASS_df['l'], MASS_df['Cltg']
logP_2MASS=0
for l in ls[2:]:
    logp=np.log(gauss(ctg[l], pdfs['mean'][l], pdfs['variance'][l]))
    logP_2MASS+=logp

#logmax=max(max(logPs))
#print(logmax)
Ps=[]

for i in range(len(lbdas)):
	Ps.append([])
	for j in range(len(z0s)):
		P=np.exp(logPs[i][j]-logP_2MASS)
		Ps[i].append(P)
		#print('i=', i, ' & j=', j)
		#print('Ps[i]=', Ps[i])
		#print('Ps=', Ps)
		
		if P==min(map(min,Ps)):
			minimum=[z0s[j], lbdas[i]]

#print(logPs)
#print(lbdas)

print('Error Points:', error_points)
print('Minimum Cell (z0,lbda)=({0},{1})'.format(minimum[0], minimum[1]))

tb.tab_export(z0s, lbdas, Ps, 'map_files/cmap_beta_{0:.1f}.txt'.format(beta))

plt.figure()
plt.title(r'Color map of $P=\prod_{\ell=2}^{128} f(C_\ell^{tg})$, where $\beta=$'+str(beta))
plt.xlabel(r'$z_0$')
plt.ylabel(r'$\lambda$')
plt.pcolormesh(z0s, lbdas, Ps, shading='nearest', vmin=0.97, vmax=1.04, cmap='RdYlBu_r')
plt.colorbar(label=r'$P/P_1$')
plt.scatter([point[0] for point in error_points], [point[1] for point in error_points], c='black', marker='x')
plt.savefig('plots/beta{0:.1f}.png'.format(beta))

