#include <iostream>
#include "math.h"

using namespace std;

void runge_kutta(double data[][100],double (*odes[])(double y[],double x),int series, int steps)
{
	
	int currSerie,currStep;
	int i,j;
	double k1,k2,k3,k4;
	double y[series],tempY[series];
	
	
	
	//runge_kutta method:
	for(currStep = 0; currStep < steps-1; currStep++){
		//initial values assigned to y[]
		for(currSerie = 0; currSerie < series; currSerie++){
			y[currSerie] = data[currSerie][currStep];
		}
		for(currSerie = 0; currSerie < series; currSerie++){
		
			k1 = (*odes[currSerie])(y,currStep);
			for(i = 0; i < series; i++){
				tempY[i] = y[i] + k1/2.;
			}
			k2 = (*odes[currSerie])(tempY,currStep+0.5);
			for(i = 0; i < series; i++){
				tempY[i] = y[i] + k2/2.;
			}
			k3 = (*odes[currSerie])(tempY,currStep+0.5);
			for(i = 0; i < series; i++){
				tempY[i] = y[i] + k3;
			}
			k4 = (*odes[currSerie])(tempY,currStep+1);
			data[currSerie][currStep+1] = y[currSerie] + 1/6.0*(k1+2*k2+2*k3+k4);
		}
	}
}


double ode1(double y[2],double x){
	return 2*x*y[0]+y[1];
}

double ode2(double y[2],double x){
	return -y[0]+2*x*y[1];
}


int main(int argc, char *argv[]) {
	
	double (*ode[2])(double y[2],double x);
	double data[2][100];
	int i,j;
	data[0][0] = 0;
	data[1][0] = 1;
	ode[0] = ode1;
	ode[1] = ode2;
	runge_kutta(data, ode, 2, 20);
	for (i = 0;i<50;i++) {
		printf("{%f,%f},",data[0][i],data[1][i]);
		
	}
	
}