#include <stdio.h>
#include <stdlib.h>
#include <math.h>

double cart_x(double r, double theta){
	return r * cos(theta);
}

double cart_y(double r, double theta){
	return r * sin(theta);
}

int main (int argc, char **argv){

//Primero, las constantes

double G = 4.51922;
double M = 1;

//Se inicializan los valores iniciales a partir de la entrada

double x0 = atof(argv[1]);
double y0 = atof(argv[2]);
double vx = atof(argv[3]) * (1.022);
double vy = atof(argv[4]) * (1.022);

//Abro un nuevo archivo para imprimir los valores iniciales

FILE * archivo;
archivo = fopen("IC.txt", "w");

//Imprimo las condiciones del objeto central

fprintf (archivo, "%d %f %f %f %f\n", -1, x0, y0, vx, vy);

//Realizo ciclos para calcular e imprimir inmediatamente las condiciones de cada masa

for(int i=0; i<5; i++){

int n_cuerpos = 12 + 6*i;
double radio = 10 + 10*i;
double velocidad = sqrt(G*M/radio);

	for(int j=0; j<n_cuerpos; j++){
		double angulo = 2 * M_PI * j/n_cuerpos;
		
		double pos_x = cart_x(radio,angulo) + x0;
		double pos_y = cart_y(radio,angulo) + y0;
		double vel_x = cart_y(velocidad,angulo) + vx;
		double vel_y = cart_x(velocidad,angulo) + vy;
		
		int id = j + (12 + 3*(i-1)) * i;
		
		fprintf (archivo, "%d %f %f %f %f\n", id, pos_x, pos_y, vel_x, vel_y);
	}
}

//Cierro mi archivo al final del programa

fclose(archivo);  
return 0;
}
