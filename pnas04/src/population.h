#ifndef POPULATION_H
#define POPULATION_H

#include "cell.h"

class Population {

    public:

		Population (const int& _ncell = 100);
		~Population ();

		/*
		 * ATT: It is the user's responsibility to ensure the
		 * format of the input dynamics is acceptable by the program.
		 * Thus, we adopt a simple format:
		 * *
		 * numr, numc
		 * tN(1), tN(2), ..., tN(numc)
		 * tT(1), tT(2), ..., tT(numc)
		 * x(1), y(1), y(2), ..., y(numc)
		 * x(2), y(1), y(2), ..., y(numc)
		 * ...
		 * x(numr), y(1), y(2), ..., y(numc)
		 * *
		 * in which numr is the number of x points and numc is the 
		 * number of y points. tN represents targetName and tT
		 * represents targetType, where 0 is for inducer and 1 is 
		 * for protein. Only dynamics of inducers and proteins are
		 * required as inputs to reverse engineer the underlying 
		 * biochemical networks. 
		 */
		void readDynamics (const string& fn);

        void growth ();
        void selection ();
   
    private:

		//	dynamics of network nodes
		int numr;	//	number of x points
		double*  xpoints;

		int numc;
		double** ypoints;	//	 number of y points per x

		//	target nodes
		string* targetName;	//	name of targets
		int* targetType;	//	0 is inducer, 1 is protein
        
		//	cells in the population
		int ncell;
        Cell** cells;
};

#endif
