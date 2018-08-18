import numpy as np
import scipy.integrate as integrate
import matplotlib.pyplot as plt

"""
Es metodo es para realizar la expolarcion sobre los datos usando MCMC

   Input:
   - Archivo que contiene los datos encontrados por Juanita Lara

   Output:
   -Graficas de los parametros 

"""

#Primero, cargamos los datos

data_obs = np.loadtxt("dimer_observations.dat")

t = data_obs[3:,0] 
p = data_obs[3:,1]

#Inicializamos algunos parametros del RK4

h = 0.125
n_points = int((25 + h)/h)
p_temp = np.zeros(n_points)
t_temp = np.zeros(n_points)
p_temp[0] = p[0]
t_temp[0] = t[0]

print p[0], t[0] #Para comprobar que leimos los datos correctamente

#Ploteamos datos iniciales 

fig = plt.figure(figsize = (12,9))
plt.plot(t, p, 'ko')
plt.plot(t, p)
fig.suptitle("Datos iniciales de P(t) y t", fontsize=20, horizontalalignment='center')
plt.xlabel("Tiempo (t)", size=16)
plt.ylabel("Sustrato (P)", size=16)

fig.savefig("Datos_iniciales.png")

"""

    Ahora definimos el runge kutta de 4 orden para integrar la funcion dada en el paper

    Input:
    - Parametros de la funcion

    Output:
    -Nuevos vectores de P y t que contendran los valores resultantes de la integral de la funcion 

"""

def rk4(alpha, beta, km1, km2, s):

    def function_dP_dt(p_tem, t_tem, alpha, beta, km1, km2, s):
        termino_1 = (alpha * (s - p_tem))/(km1 + s - p_tem)
        termino_2 = (beta * p_tem)/(km2 + p_tem)
        return termino_1 - termino_2
    
    for i in range(1, n_points):
  
        k1 = function_dP_dt(p_temp[i-1], t_temp[i-1], alpha, beta, km1, km2, s)
    
        #first step
        t1 = t_temp[i-1] + (h/2.0)
        p1 = p_temp[i-1] + (h/2.0)*k1
        k2 = function_dP_dt(p1, t1, alpha, beta, km1, km2, s)
    
        #second step
        t2 = t_temp[i-1] + (h/2.0)
        p2 = p_temp[i-1] + (h/2.0)*k2
        k3 = function_dP_dt(p2, t2, alpha, beta, km1, km2, s)
        
        #third step
        t3 = t_temp[i-1] + h
        p3 = p_temp[i-1] + h*k3
        k4 = function_dP_dt(p3, t3, alpha, beta, km1, km2, s)
        
        #fourth step
        average_k = (1.0/6.0)*(k1 + 2.0*k2 + 2.0*k3 + k4)
        
        t_temp[i] = t_temp[i-1] + h
        p_temp[i] = p_temp[i-1] + h*average_k
        
    return t_temp, p_temp

"""

   Ahora definimos el likelihood que contendra el chi cuadrado para cada P

   Input:
   - Los datos iniciales de P y los que calcularemos con el modelo

   Output:
   - Chi cuadrado

"""
    
def likelihood(p, p_model):   
    chi_2 = 0.0
    contador = 0
    for i in range(len(p_model)):
        for j in range(len(p)):
            chi_2 += (p[j] - p_model[i])**2
            contador += 1
    return chi_2

"""

   Procedemos a realizar la exploracion con el algoritmo de metropolis-hastings

"""

#Inicializacion de vectores para la exploracion

alpha_walk = np.empty((0))
beta_walk = np.empty((0))
km1_walk = np.empty((0))
km2_walk = np.empty((0))
s_walk = np.empty((0))
los_sigmas = np.empty((0))

alpha_walk = np.append(alpha_walk, np.random.random())
beta_walk = np.append(beta_walk, np.random.random())
km1_walk = np.append(km1_walk, np.random.random())
km2_walk = np.append(km2_walk, np.random.random())
s_walk = np.append(s_walk, np.random.random())

#Exploracion

n_iterations = 2000

for i in range(n_iterations):
    
    alpha_prime = np.abs(np.random.normal(alpha_walk[i], 0.1)) 
    beta_prime = np.abs(np.random.normal(beta_walk[i], 0.1))
    km1_prime = np.abs(np.random.normal(km1_walk[i], 0.1))
    km2_prime = np.abs(np.random.normal(km2_walk[i], 0.1))
    s_prime = np.abs(np.random.normal(s_walk[i], 0.1))

    t_init, p_init = rk4(alpha_walk[i], beta_walk[i], km1_walk[i], km2_walk[i], s_walk[i])
    t_prime, p_prime = rk4(alpha_prime, beta_prime, km1_prime, km2_prime, s_prime)

    sigma_1 = likelihood(p, p_init)
    sigma_2 = likelihood(p, p_prime)
    
    sigma = sigma_1 - sigma_2
    
    if(sigma >= 0.0):
        
        alpha_walk  = np.append(alpha_walk, alpha_prime)
        beta_walk  = np.append(beta_walk, beta_prime)
        km1_walk  = np.append(km1_walk, km1_prime)
        km2_walk  = np.append(km2_walk, km2_prime)
        s_walk  = np.append(s_walk, s_prime)
        los_sigmas = np.append(los_sigmas, sigma_2)
        
    else:
        
        otro = np.random.random()
        
        if(otro <= exp(-(0.5*(sigma)))):
            
            alpha_walk = np.append(alpha_walk, alpha_prime)
            beta_walk = np.append(beta_walk, beta_prime)
            km1_walk = np.append(km1_walk, km1_prime)
            km2_walk = np.append(km2_walk, km2_prime)
            s_walk = np.append(s_walk, s_prime)
            los_sigmas = np.append(los_sigmas, sigma_2)
            
        else:
            
            alpha_walk = np.append(alpha_walk, alpha_walk[i])
            beta_walk = np.append(beta_walk, beta_walk[i])
            km1_walk = np.append(km1_walk, km1_walk[i])
            km2_walk = np.append(km2_walk, km2_walk[i])
            s_walk = np.append(s_walk, s_walk[i])
            los_sigmas = np.append(los_sigmas, sigma_1)

#Imprimimos los parametros (alpha_walk, beta_walk, km1_walk, km2_walk, s_walk) en un archivo para luego ser graficados

np.savetxt('datos_arreglados.dat', zip(alpha_walk, beta_walk, km1_walk, km2_walk, s_walk,los_sigmas))

#Ahora sacamos los valores minimos para alpha, beta, km1, km2 y s, y los imprimimos

h = np.argmin(los_sigmas)

alpha_min = alpha_walk[h]
beta_min = beta_walk[h]
km1_min = km1_walk[h]
km2_min = km2_walk[h]
s_min = s_walk[h]
los_sigmas_min = los_sigmas[h]

print "Valores minimos de los parametros:"
print "alpha = " + str(alpha_walk[h])
print "beta = " + str(beta_walk[h])
print "km1 = " + str(km1_walk[h])
print "km2 = " + str(km2_walk[h])
print "s = " + str(s_walk[h])
print "Chi_2 = " + str(los_sigmas[h])

f = open('Valores_parametros_min.dat', 'w')
f.write("alpha_minimo = " + str(alpha_min))
f.write(" ")
f.write("beta_minimo = " + str(beta_min))
f.write(" ")
f.write("km1_minimo = " + str(km1_min))
f.write(" ")
f.write("km2_minimo = " + str(km2_min))
f.write(" ")
f.write("s_minimo = " + str(s_min))
f.write(" ")
f.write("Chi_2_minimo = " + str(los_sigmas_min))
f.close()

#Graficamos el modelo con los valores minimos y los datos iniciales

best_t, best_p = rk4(alpha_walk[h], beta_walk[h], km1_walk[h], km2_walk[h], s_walk[h]) 

fig2 = plt.figure(figsize = (12,9))
plt.plot(t, p, 'ko')
plt.plot(t, p, label = "Datos iniciales")
plt.plot(best_t, best_p, label = "Modelo aproximado") 
fig2.suptitle("Datos iniciales contra el modelo", fontsize=20, horizontalalignment='center')
plt.xlabel("Tiempo (t)", size=16)
plt.ylabel("Sustrato (P)", size=16)
plt.axis([5, 28, 17, 35])

fig2.savefig("Datos_modelo_min.png")
            
#Ahora determinamos los valores promedio de los parametros, los imprimimos y los graficamos

average_alpha = np.average(alpha_walk)
average_beta = np.average(beta_walk)
average_km1 = np.average(km1_walk)
average_km2 = np.average(km2_walk)
average_s = np.average(s_walk)
average_los_sigmas = np.average(los_sigmas)

print "Valores promedio de los parametros:"
print "alpha = " + str(average_alpha)
print "beta = " + str(average_beta)
print "km1 = " + str(average_km1)
print "km2 = " + str(average_km2)
print "s = " + str(average_s)
print "Chi_2 = " + str(average_los_sigmas)

f = open('Valores_parametros_prom.dat', 'w')
f.write("alpha_promedio = " + str(average_alpha))
f.write(" ")
f.write("beta_promedio = " + str(average_beta))
f.write(" ")
f.write("km1_promedio = " + str(average_km1))
f.write(" ")
f.write("km2_promedio = " + str(average_km2))
f.write(" ")
f.write("s_promedio = " + str(average_s))
f.write(" ")
f.write("Chi_2_promedio = " + str(average_los_sigmas))
f.close()

tiempo, mejor_p = rk4(average_alpha, average_beta, average_km1, average_km2, average_s)

fig3 = plt.figure(figsize = (12,9))
plt.plot(t, p, 'ko')
plt.plot(t, p, label = "Datos iniciales")
plt.plot(tiempo, mejor_p, label = "Modelo aproximado") 
fig3.suptitle("Datos iniciales contra el modelo", fontsize=20, horizontalalignment='center')
plt.xlabel("Tiempo (t)", size=16)
plt.ylabel("Sustrato (P)", size=16)
plt.axis([5, 28, 17, 35])

fig3.savefig("Datos_modelo_prom.png")
