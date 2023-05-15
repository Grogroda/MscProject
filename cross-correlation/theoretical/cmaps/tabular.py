import pandas as pd
import random
import matplotlib.pyplot as plt

'''
This is a module I'm using to implement 2 functionalities specifically: One to export a well formatted table of colors into a text format available for the user to look and save the exact values obtained in the map, and another to do the exact opposite: import the file and automatically convert the file into correctly formatted lists/matrices with the values desired.
'''

def tab_export(X,Y,colors, file_name='my_table.txt'):

	table=pd.DataFrame()
	
	for i in range(len(X)):
		x=X[i]
		col=[]
		for j in range(len(Y)):
			col.append(colors[j][i])

		table[x]=col

	table=table.set_index(pd.Index(Y))
	table.to_csv(file_name, sep=' ')
	
def tab_import(file_name='my_table.txt'):

	table=pd.read_csv(file_name, sep=' ', index_col=0)
	
	Xstr=table.columns.to_list()
	X=[float(x) for x in Xstr]
	Y=table.index.to_list()
	
	colors=[]
	
	for y in Y:
		sublist=[]
		for x in X:
			col=table[str(x)][y]
			sublist.append(col)
			
		colors.append(sublist)

	return X,Y,colors

if __name__=='__main__':

	X=[0.5,0.6,0.7,0.8]
	Y=[1,1.2,1.4,1.6, 1.8]
	colors=[[random.random() for i in range(len(X))] for j in range(len(Y))]

	print('X=', X)
	print('Y=', Y, '\n')
	print('CMap Matrix:')
	for line in range(len(colors)):
		print(colors[line])

	plt.figure()
	plt.title('Color Map before export')
	plt.pcolormesh(X,Y,colors,shading='nearest')
	plt.colorbar()
	plt.show()

	tab_export(X,Y,colors,'test.txt')
	newX, newY, newColors=tab_import('test.txt')
	
	print('\nImported X-axis:')
	print(newX)
	
	print('\nImportex Y-axis:')
	print(newY)
	
	print('\nNew color matrix:')
	for line in range(len(newColors)):
		print(newColors[line])
