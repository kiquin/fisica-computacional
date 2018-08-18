#include <stdio.h>
#include <stdlib.h>
#include <math.h>

double f_prime1 (double velocidad){
	return velocidad;
}

double f_prime2 (int serie, double pos1, double pos2, double *mpos1, double *mpos2, double G){
	double respuesta;
	
	respuesta = -G * (pos1 - mpos1[0]) / pow( pow(pos1 - mpos1[0],2.0) + pow(pos2 - mpos2[0],2.0) ,1.5);
	
	if (serie == 2){
		respuesta += -G * (pos1 - mpos1[1]) / pow( pow(pos1 - mpos1[1],2.0) + pow(pos2 - mpos2[1],2.0) ,1.5);
	}
	return respuesta;
}

double f_masa (int serie, double mx1, double my1, double mx2, double my2, double G){
	if(serie == 2){
		return -G * (mx1 - mx2) / pow( pow(mx1 - mx2,2.0) + pow(my1 - my2,2.0) ,1.5);
	}
	else{
	 	return 0;
	}
}

void rk4step (double * x, double * vx, double * y, double * vy, double * xm, double * ym, int n, int serie, double h, double G){

	for (int i=0; i<n; i++){

	//Primer paso, para X y Y
	double xk1_1 = f_prime1(vx[i]);
	double xk1_2 = f_prime2(serie, x[i], y[i], xm, ym, G);
	double yk1_1 = f_prime1(vy[i]);
	double yk1_2 = f_prime2(serie, y[i], x[i], ym, xm, G);
	
	//Segundo paso, calculando primero posiciones y velocidades, luego los valores de K
	double x1_1 = x[i] + (h/2.0) * xk1_1;
	double x1_2 = vx[i] + (h/2.0) * xk1_2;
	double y1_1 = y[i] + (h/2.0) * yk1_1;
	double y1_2 = vy[i] + (h/2.0) * yk1_2;
	
	double xk2_1 = f_prime1(x1_2);	
	double xk2_2 = f_prime2(serie, x1_1, y1_1, xm, ym, G);	
	double yk2_1 = f_prime1(y1_2);	
	double yk2_2 = f_prime2(serie, y1_1, x1_1, ym, xm, G);	
	
	//Tercer Paso
	double x2_1 = x[i] + (h/2.0) * xk2_1;
	double x2_2 = vx[i] + (h/2.0) * xk2_2;
	double y2_1 = y[i] + (h/2.0) * yk2_1;
	double y2_2 = vy[i] + (h/2.0) * yk2_2;
		
	double xk3_1 = f_prime1(x2_2);	
	double xk3_2 = f_prime2(serie, x2_1, y2_1, xm, ym, G);	
	double yk3_1 = f_prime1(y2_2);	
	double yk3_2 = f_prime2(serie, y2_1, x2_1, ym, xm, G);	
    
    //Cuarto paso
	double x3_1 = x[i] + h * xk3_1;
	double x3_2 = vx[i] + h * xk3_2;
	double y3_1 = y[i] + h * yk3_1;
	double y3_2 = vy[i] + h * yk3_2;
	
	double xk4_1 = f_prime1(x3_2);	
	double xk4_2 = f_prime2(serie, x3_1, y3_1, xm, ym, G);	
	double yk4_1 = f_prime1(y3_2);	
	double yk4_2 = f_prime2(serie, y3_1, x3_1, ym, xm, G);	
    
    //Calculo de promedios
    double xk_promedio_1 = (1.0/6.0)*(xk1_1 + 2.0*xk2_1 + 2.0*xk3_1 + xk4_1);
    double xk_promedio_2 = (1.0/6.0)*(xk1_2 + 2.0*xk2_2 + 2.0*xk3_2 + xk4_2);
    
    double yk_promedio_1 = (1.0/6.0)*(yk1_1 + 2.0*yk2_1 + 2.0*yk3_1 + yk4_1);
    double yk_promedio_2 = (1.0/6.0)*(yk1_2 + 2.0*yk2_2 + 2.0*yk3_2 + yk4_2);

	//Actualizo mis vectores
    x[i] = x[i] + h * xk_promedio_1;
    vx[i] = vx[i] + h * xk_promedio_2;
    
    y[i] = y[i] + h * yk_promedio_1;
    vy[i] = vy[i] + h * yk_promedio_2;
	}
}

void rkmasa (double *xm, double*ym, double *vxm, double *vym, int serie, double h, double G){

	//Primer paso, para X y Y de la masa 1
	double xk1_1 = f_prime1(vxm[0]);
	double xk1_2 = f_masa(serie, xm[0],  ym[0],  xm[1], ym[1], G);
	double yk1_1 = f_prime1(vym[0]);
	double yk1_2 = f_masa(serie, ym[0],  xm[0],  ym[1], xm[1], G);
	
	//Primer paso, para X y Y de la masa 2
	double bxk1_1 = f_prime1(vxm[1]);
	double bxk1_2 = f_masa(serie, xm[1],  ym[1],  xm[0], ym[0], G);
	double byk1_1 = f_prime1(vym[1]);
	double byk1_2 = f_masa(serie, ym[1],  xm[1],  ym[0], xm[0], G);
	
	//Segundo paso, calculando primero posiciones y velocidades, luego los valores de K
	double x1_1 = xm[0] + (h/2.0) * xk1_1;
	double x1_2 = vxm[0] + (h/2.0) * xk1_2;
	double y1_1 = ym[0] + (h/2.0) * yk1_1;
	double y1_2 = vym[0] + (h/2.0) * yk1_2;
	
	double bx1_1 = xm[1] + (h/2.0) * xk1_1;
	double bx1_2 = vxm[1] + (h/2.0) * xk1_2;
	double by1_1 = ym[1] + (h/2.0) * yk1_1;
	double by1_2 = vym[1] + (h/2.0) * yk1_2;
	
	double xk2_1 = f_prime1(x1_2);	
	double xk2_2 = f_masa(serie, x1_1,  y1_1,  bx1_1, by1_1, G);
	double yk2_1 = f_prime1(y1_2);	
	double yk2_2 = f_masa(serie, y1_1,  x1_1,  by1_1, bx1_1, G);
	
	double bxk2_1 = f_prime1(bx1_2);	
	double bxk2_2 = f_masa(serie, bx1_1,  by1_1,  x1_1, y1_1, G);
	double byk2_1 = f_prime1(by1_2);	
	double byk2_2 = f_masa(serie, by1_1,  bx1_1,  y1_1, x1_1, G);
	
	//Tercer Paso
	double x2_1 = xm[0] + (h/2.0) * xk2_1;
	double x2_2 = vxm[0] + (h/2.0) * xk2_2;
	double y2_1 = ym[0] + (h/2.0) * yk2_1;
	double y2_2 = vym[0] + (h/2.0) * yk2_2;
	
	double bx2_1 = xm[1] + (h/2.0) * bxk2_1;
	double bx2_2 = vxm[1] + (h/2.0) * bxk2_2;
	double by2_1 = ym[1] + (h/2.0) * byk2_1;
	double by2_2 = vym[1] + (h/2.0) * byk2_2;
		
	double xk3_1 = f_prime1(x2_2);	
	double xk3_2 = f_masa(serie, x2_1,  y2_1,  bx2_1, by2_1, G);
	double yk3_1 = f_prime1(y2_2);	
	double yk3_2 = f_masa(serie, y2_1,  x2_1,  by2_1, bx2_1, G);
    
    double bxk3_1 = f_prime1(bx2_2);	
	double bxk3_2 = f_masa(serie, bx2_1,  by2_1,  x2_1, y2_1, G);
	double byk3_1 = f_prime1(by2_2);	
	double byk3_2 = f_masa(serie, by2_1,  bx2_1,  y2_1, x2_1, G);
	
    //Cuarto paso
	double x3_1 = xm[0] + h * xk3_1;
	double x3_2 = vxm[0] + h * xk3_2;
	double y3_1 = ym[0] + h * yk3_1;
	double y3_2 = vym[0] + h * yk3_2;
	
	double bx3_1 = xm[1] + h * bxk3_1;
	double bx3_2 = vxm[1] + h * bxk3_2;
	double by3_1 = ym[1] + h * byk3_1;
	double by3_2 = vym[1] + h * byk3_2;
	
	double xk4_1 = f_prime1(x3_2);	
	double xk4_2 = f_masa(serie, x3_1,  y3_1,  bx3_1, by3_1, G);
	double yk4_1 = f_prime1(y3_2);	
	double yk4_2 = f_masa(serie, y3_1,  x3_1,  by3_1, bx3_1, G);
	
	double bxk4_1 = f_prime1(bx3_2);	
	double bxk4_2 = f_masa(serie, bx3_1,  by3_1,  x3_1, y3_1, G);
	double byk4_1 = f_prime1(y3_2);	
	double byk4_2 = f_masa(serie, by3_1,  bx3_1,  y3_1, x3_1, G);
    
    //Calculo de promedios
    double xk_promedio_1 = (1.0/6.0)*(xk1_1 + 2.0*xk2_1 + 2.0*xk3_1 + xk4_1);
    double xk_promedio_2 = (1.0/6.0)*(xk1_2 + 2.0*xk2_2 + 2.0*xk3_2 + xk4_2);
    
    double yk_promedio_1 = (1.0/6.0)*(yk1_1 + 2.0*yk2_1 + 2.0*yk3_1 + yk4_1);
    double yk_promedio_2 = (1.0/6.0)*(yk1_2 + 2.0*yk2_2 + 2.0*yk3_2 + yk4_2);
    
    double bxk_promedio_1 = (1.0/6.0)*(bxk1_1 + 2.0*bxk2_1 + 2.0*bxk3_1 + bxk4_1);
    double bxk_promedio_2 = (1.0/6.0)*(bxk1_2 + 2.0*bxk2_2 + 2.0*bxk3_2 + bxk4_2);
    
    double byk_promedio_1 = (1.0/6.0)*(byk1_1 + 2.0*byk2_1 + 2.0*byk3_1 + byk4_1);
    double byk_promedio_2 = (1.0/6.0)*(byk1_2 + 2.0*byk2_2 + 2.0*byk3_2 + byk4_2);

	//Actualizo mis vectores
    xm[0] = xm[0] + h * xk_promedio_1;
    vxm[0] = vxm[0] + h * xk_promedio_2;
    
    ym[0] = ym[0] + h * yk_promedio_1;
    vym[0] = vym[0] + h * yk_promedio_2;
    
    xm[1] = xm[1] + h * bxk_promedio_1;
    vxm[1] = vxm[1] + h * bxk_promedio_2;
    
    ym[1] = ym[1] + h * byk_promedio_1;
    vym[1] = vym[1] + h * byk_promedio_2;
}

int main(int argc, char **argv) {

//El programa revisa que se utilice el numero correcto de argumento
	if( argc != 2){
  		printf("Este programa necesita un unico archivo de datos para ser ejecutado, intentelo de nuevo");
 		exit(1);
	}

//Abre el archivo de condiciones iniciales 
	FILE *in;

	if(!(in=fopen(argv[1], "r"))){
		printf("Ha ocurrido un problema al abrir el archivo %s\n", argv[1]);
		exit(1);
	}

//Cuento las lineas del archivo
	int n_lines = 0;
	int c;
 	do{
    	c = fgetc(in);
    	if(c=='\n'){
      	n_lines++;
    	}
  	}while(c!=EOF);

	rewind(in);
	
//Inicializo los vectores para guardar la información de cada cuerpo y las masas
	double G = 4.51922;
	
	double * x = calloc(n_lines, sizeof(double));
	double * y = calloc(n_lines, sizeof(double)); 
	double * vx = calloc(n_lines, sizeof(double));
	double * vy = calloc(n_lines, sizeof(double));

	double * xm = calloc(2, sizeof(double));
	double * ym = calloc(2, sizeof(double));
	double * vxm = calloc(2, sizeof(double));
	double * vym = calloc(2, sizeof(double));
	
	int i;
	
	int rid;
	double rx;
	double ry;
	double rvx;
	double rvy;
	
//Leo el archivo y guardo los datos en los vectores
	int serie = 0;
	int n_objetos = 0;
	
	for(i=0; i<n_lines; i++){
	
		fscanf (in, "%d %lf %lf %lf %lf\n", &rid, &rx, &ry, &rvx, &rvy);
		
		if (rid == -1){
			xm[serie] = rx;
			ym[serie] = ry;
			vxm[serie] = rvx;
			vym[serie] = rvy;
			
			serie++;
		}
		else{
			x[n_objetos] = rx;
			y[n_objetos] = ry;
			vx[n_objetos] = rvx;
			vy[n_objetos] = rvy;
			n_objetos++;
		}
	}

//Cierra el archivo al finalizar
	fclose(in);

//	printf ("%lf %lf %lf %lf %d\n", x2[26], y2[26], vx2[26], vy2[26], serie);
double h = 0.1;

FILE * archivo;

char myfile[256];

for(int j=0; j<50000; j++){
	rk4step (x, vx, y, vy, xm, ym, n_objetos, serie, h, G);

	rkmasa (xm, ym, vxm, vym, serie, h,  G);

	if((j+1)%10000 == 0){
		int z;
		z = sprintf(myfile, "output_%d.dat", (j+1)/10 );
		archivo = fopen(myfile, "w");
		
		if(serie == 2){
			fprintf(archivo,"%d %lf %lf %lf %lf\n", -2, xm[1],  ym[1], vxm[1], vym[1]);
		}
		
		fprintf(archivo,"%d %lf %lf %lf %lf\n", -1, xm[0],  ym[0], vxm[0], vym[0]);
		
		for(int k=0; k<n_objetos; k++){
			fprintf(archivo,"%d %lf %lf %lf %lf\n", k, x[k],  y[k], vx[k], vy[k]);
		}
	}
}
   
fclose(archivo);

	return 0;
}



