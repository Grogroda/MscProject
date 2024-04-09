from matplotlib import pyplot as plt

w0=-1

#Para cs2 constante:
wa_list=[0,0.33,0.66,1]
z0,zmax,N_it=0,0.2,1000
dz=(zmax-z0)/N_it
z_list=[z0+dz*i for i in range(N_it)]

def w(wa,z):

	return w0+(1-1/(1+z))*wa

plt.figure(1)
plt.title(r'Gráficos para $w(a)=w+(1-1/(1+z))w_a$')
plt.xlabel('z')
plt.ylabel('w(z)')
plt.figure(2)
plt.title('Diferenças Relativas ao LCDM')
plt.xlabel('z')
plt.ylabel('(wLCDM-w)/wLCDM')

for wa in wa_list:
	w_list=[w(wa,z) for z in z_list]
	plt.figure(1)
	plt.plot(z_list, w_list, label=r'$w_a={0}$'.format(wa))
	
	if wa==0:
		w_LCDM=w_list
	
	dif=[(w_LCDM[i]-w_list[i])/w_LCDM[i] for i in range(len(w_list))]
	plt.figure(2)
	plt.plot(z_list,dif, label=r'$w_a={0}$'.format(wa))

plt.figure(1)
plt.legend()
plt.savefig('wz.png')

plt.figure(2)
plt.legend()
plt.savefig('wz_dif.png')
	
