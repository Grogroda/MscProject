import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import tabular as tb
#import HistSearch as hs

z0=float(input('Value of z0: '))

def gauss(x,mean,var): #returns f(x) where f is a normalized gaussian distribution with mean=mean and variance=var

	return 1/(np.sqrt(2*np.pi*var))*np.exp(-(x-mean)**2/(2*var))

ls=[i for i in range(129)]

lbdamin, lbdamax, dlbda=0.8, 8.0, 0.2
Nlbda=round((lbdamax-lbdamin)/dlbda)

betamin, betamax, dbeta=0.6, 3.8, 0.1
Nbeta=round((betamax-betamin)/dbeta)

betas, lbdas=[betamin+dbeta*i for i in range(Nbeta)],[lbdamin+dlbda*i for i in range(Nlbda)]
pdfs=pd.read_csv('../tables/hists/gauss_fits.txt', sep=' ', header=None, names=['l', 'mean', 'variance'])
#print(pdfs) 
logPs=[]
error_points=[]

for lbda in lbdas: #lbda vai ser o eixo Y da matriz de cores
    sublist=[] #cria uma nova linha para a matriz de cores
    for beta in betas: #z0 vai ser o eixo X da matriz de cores
    	logP=0
    	Cls=pd.read_csv('../ctg_files/funcs_z0_{0:.2f}/beta_{1:.1f}.lbda_{2:.1f}.dat'.format(z0,beta, lbda), sep=' ', header=None, names=['l', 'Cltg'])
        #print(Cls.to_string())
        
    	if len(Cls['Cltg'])==129:
            for l in ls[2:]:
                Ctg=Cls['Cltg'][l]
                logp=np.log(gauss(Ctg, pdfs['mean'][l], pdfs['variance'][l]))
                logP=logP+logp
            
    	else:
            print("(beta,lbda)=({0},{1})".format(beta,lbda))
            print('len(l)=', len(Cls['l']))    	
            logP=1#max(max(logPs))
            error_points.append((beta,lbda))

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
	for j in range(len(betas)):
		P=np.exp(logPs[i][j]-logP_2MASS)
		Ps[i].append(P)
		
		if P==min(map(min,Ps)):
			minimum=[betas[j], lbdas[i]]		

#print(logPs)
#print(lbdas)

print('Error Points:', error_points)
print('Minimum Cell (lbda,beta)=({0},{1})'.format(minimum[0], minimum[1]))

tb.tab_export(betas, lbdas, Ps, 'map_files/cmap_z0_{0:.2f}.txt'.format(z0))

plt.figure()
plt.title(r'Color map of $P=\prod_{\ell=2}^{128} f(C_\ell^{tg})$, where $z_0=$'+str(z0))
plt.xlabel(r'$\beta$')
plt.ylabel(r'$\lambda$')
plt.pcolormesh(betas, lbdas, Ps, shading='nearest', vmin=0.97, vmax=1.04, cmap='RdYlBu_r')
plt.colorbar(label=r'$P/P_1$')
plt.scatter([point[0] for point in error_points], [point[1] for point in error_points], c='black', marker='x')
plt.savefig('plots/z0_{0:.2f}.png'.format(z0))
