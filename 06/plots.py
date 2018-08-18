from matplotlib import pyplot
import numpy as np

#Utilizando la funcion loadtxt de numpy, se cargan los contenidos del archivo especificado
data1 = np.loadtxt('output_1000.dat')
data2 = np.loadtxt('output_2000.dat')
data3 = np.loadtxt('output_3000.dat')
data4 = np.loadtxt('output_4000.dat')
data5 = np.loadtxt('output_5000.dat')

#Usando pyplot, se dibuja la grafica para las columnas que contienen los valores de x y y, y luego se muestra
pyplot.figure(figsize=(15,3))
pyplot.subplot(151)
pyplot.scatter(data1[:,1], data1[:,2], 1),

pyplot.subplot(152)
pyplot.scatter(data2[:,1], data2[:,2], 1)

pyplot.subplot(153)
pyplot.scatter(data3[:,1], data3[:,2], 1)

pyplot.subplot(154)
pyplot.scatter(data4[:,1], data4[:,2], 1)

pyplot.subplot(155)
pyplot.scatter(data5[:,1], data5[:,2], 1)

pyplot.show()