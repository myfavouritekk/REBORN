#include <iostream>
#include "math.h"
#define MAXSTEPS 1000

using namespace std;
void runge_kutta(double data[][MAXSTEPS], double (*odes[])(double y[],double x), int series, int steps)
{//100 is steps
	
	int currSerie, currStep;
	int i;
	double k1, k2, k3, k4;
	double *y = new double[series], *tempY = new double[series];
	
	
	
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
    delete y;
    delete tempY;
}
