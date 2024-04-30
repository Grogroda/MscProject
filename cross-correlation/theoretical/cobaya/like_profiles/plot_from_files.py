import numpy as np
from matplotlib import pyplot as plt
import pandas as pd

Omegas1, P1, expXi1=[],[],[]

for i in [1,2,3,4]:
    data=pd.read_csv('ProfileData{}_both_band1_Nmc2e+07.dat'.format(i), sep=' ', index_col=0)
    Omegas1+=list(data['Omega'])
    P1+=list(data['P'])
    expXi1+=list(data['expXi2'])

print('Omegas=', Omegas1)
print('P=', P1)

plt.figure()
plt.plot(Omegas1, P1)
plt.show()
    
