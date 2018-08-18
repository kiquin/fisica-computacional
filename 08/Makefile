all : datos_arreglados.dat

datos_arreglados.dat: param

param:Exploracion_MCMC.py
		python Exploracion_MCMC.py
		python contour.py datos_arreglados.dat
		
clean:
		rm -f datos_arreglados.dat *.png Valores_parametros_prom.dat Valores_parametros_min.dat
