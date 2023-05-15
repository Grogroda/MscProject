import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

direct=input('Cltg_file (with relative or absolute directory location): ')
df=pd.read_csv(direct, sep=' ', header=None, names=['l', 'cltg'])
ls, ctg=df['l'], df['cltg']
Dtg=[ls[i]*(ls[i]+1)*ctg[i]/(2*np.pi) for i in range(len(ls))]

#for l in range(1,11):
#    Dtg[l]=Dtg[l]/(0.3*l)

save=input('Save the image (sa) or just show it (sh)? ')

plt.figure()
plt.plot(ls[2:], Dtg[2:])
plt.xlabel(r'$\ell$')
plt.ylabel(r'$\ell(\ell+1)C_\ell^{tg}/2\pi$')
plt.xscale('log')

if save=='sa':
    file_name=input('Image file name (with image extension, ex: image.png): ')
    plt.savefig(file_name)

elif save=='sh':
    plt.show()

else:
    print("Well, I guess that's it then, bye.")
