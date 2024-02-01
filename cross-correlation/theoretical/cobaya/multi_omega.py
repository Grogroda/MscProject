from pyctg import ctg4py
from tqdm import tqdm
from matplotlib import pyplot as plt

omegas=[0.2,0.3,0.45,0.6,0.7]

plots={} #{omegaM1:[ls, ctgs]} where ls and ctgs are lists containing multipolesl and ctg(l) to plot

fig, ax=plt.subplots()
ax.set_xlabel(r'$\ell$')
ax.set_ylabel(r'$C_{tg}$')
ax.set_xscale('log')

for i in tqdm(range(len(omegas))):
    omegaM=omegas[i]
    plots[omegaM]=[]
    
    ls, ctg=ctg4py(omegaM)

    ax.plot(ls, ctg, label=r'$\Omega_M={}$'.format(omegaM))
    
    plots[omegaM].append(ls)
    plots[omegaM].append(ctg)

ax.legend()
plt.savefig('multi_omegas.png')

#print(plots)
