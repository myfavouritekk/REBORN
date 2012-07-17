#ifndef POPULATION_H
#define POPULATION_H

#include "cell.h"
#include "scorefunction.h"
#include <iostream>
#include <stdexcept>
#include <fstream>
#include <libxml/tree.h>
#include <libxml/parser.h>
#pragma comment(lib,"libxml2.lib")

class Population {

    public:

		Population (const int& _ncell = 100, const int& _evolution = 1000);//_ncell used as a reference to 100 which is the number of total cells
		~Population ();
    
    

		/*
		 * ATT: It is the user's responsibility to ensure the
		 * format of the input dynamics is acceptable by the program.
		 * Thus, we adopt a simple format:
		 * *
		 * numr, numind, numprot
		 * tN(1), tN(2), ..., tN(numc)
		 * tT(1), tT(2), ..., tT(numc)
		 * x(1), y(1), y(2), ..., y(numc)
		 * x(2), y(1), y(2), ..., y(numc)
		 * ...
		 * x(numr), y(1), y(2), ..., y(numc)
		 * *
		 * in which numr is the number of x points and numc (=numind+numprot) is the 
		 * number of y points. tN represents targetName and tT
		 * represents targetType, where 0 is for inducer and 1 is 
		 * for protein. Only dynamics of inducers and proteins are
		 * required as inputs to reverse engineer the underlying 
		 * biochemical networks. 
		 */
		void readDynamics (const string& fn);
        
        void readDynamicsFromConsole();

        //initializer
        void init();
    
        //get evolution number
        int getEvolution();
    
        //for growth phase
        void growth ();
    
        //population mutation for topology
        void mutation();
    
        //pupulation mutation for kinects
        void mut_parameters();
    
        //sort all cells after mutate parameters for 100 generations
        void sort();
    
        //for selection phase
        void selection ();
    
        //classify the surviving cells
        void classification();
    
        //judging whether the evolution is terminated
        bool isTerminate();
    
        //generating certain output format
        void genTikzFormat();
    
        //generating output
        void output();
    
        //generating XML output format
        void genXMLFormat();
       
    private:

		//	dynamics of network nodes
		int numr;               //	number of x points
		double*  xpoints; //dynamic array containing the xpoints

		int numind, numprot;    // number of inducers and proteins
		double** ypoints;       //	 target dynamic data

		int ncell;              //number of cells in the population
        Cell** cells;           //set of pointers of cells
    
        //classification of cells
        std::vector<std::vector<Cell*> > classifiedCells;
    
        int evolution;
    
        //score function for sorting
        ScoreFunc sfunc;
};

#endif
