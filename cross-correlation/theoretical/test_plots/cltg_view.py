import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

cont=True
c=0

fig, ax=plt.subplots()

while cont:

    direct=input('Cltg_file (with relative or absolute path): ')
    df=pd.read_csv(direct, sep=' ', header=None, names=['l', 'cltg'])
    ls, ctg=df['l'], df['cltg']
    Dtg=[ls[i]*(ls[i]+1)*ctg[i]/(2*np.pi) for i in range(len(ls))]

    lab=input('Label: ')
    
    ax.plot(ls[2:], Dtg[2:], label=lab)

    c+=1
    check=input('Continue plotting? (Y/n) ')

    if check=='n':
        cont=False

if c>1:
    ax.legend()

ax.set_xlabel(r'$\ell$')
ax.set_ylabel(r'$\ell(\ell+1)C_\ell^{tg}/2\pi$')
ax.set_xscale('log')

save=input('Save (sa) or just show (sh)? ')

if save=='sa':
    file_name=input('Image file name (with image extension, ex: image.png): ')
    plt.savefig(file_name)

elif save=='sh':
    plt.show()

else:
    print("Well, I guess that's it then, bye.")
