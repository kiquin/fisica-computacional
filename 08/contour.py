import scipy.interpolate
import numpy as np
import matplotlib.pyplot as plt
import sys

#Carga el archivo producido por el programa principal
data = np.loadtxt(str(sys.argv[1]))

#Obtiene los valores para cada parametro
alfa = data[:,0]
beta = data[:,1]
km1 = data[:,2]
km2 = data[:,3]
S = data[:,4]
chi2 = data[:,5]

#Creo arrays con titulos y valores para facil manejo
a = ['Alfa','Beta','Km1','Km2','S']
b = [alfa,beta,km1,km2,S]

#Metodo iterativo que grafica los contornos de Chi2 en todas las combinaciones de parametros
for i in range(len(b)):
    for j in range(len(b)):
        if i!=j:
        
        	#Creo arreglos para interpolar
            x = np.linspace(b[i].min(), b[i].max(), 500)
            y = np.linspace(b[j].min(), b[j].max(), 500)
            x, y = np.meshgrid(x, y)
            
            #Interpolo en 2D utilizando el metodo de triangulacion griddata
            z2 = scipy.interpolate.griddata((b[i], b[j]), chi2, (x, y), method='cubic')
            title = "Distribucion de chi2 en el plano " + a[j] + '-' + a[i]
            
            #Grafico la figura, para todos los valores en mi distribucion de chi2, en el rango de valores de los parametros usados
            fig = plt.figure(figsize=(12,9))
            fig.suptitle(title, fontsize=20, horizontalalignment='center')
            plt.xlabel(a[i], size=16)
            plt.ylabel(a[j], size=16)
            plt.imshow(z2, vmin=chi2.min(), vmax=chi2.max(), origin='lower',extent=[b[i].min(), b[i].max(), b[j].min(), b[j].max()])
            plt.colorbar()
            
            #Guardo cada archivo como una imagen PNG
            filename = a[j] + '_vs_' + a[i] + '.png'
            fig.savefig(filename)
