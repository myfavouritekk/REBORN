#include "population.h"

Population::Population (const int& _ncell)
	:ncell(_ncell), numx(0), numind(0), numprot(0), xpoints(NULL), ypoints(NULL), cells(NULL) 
{
	ncell = ncell <=0 ? 100 : ncell;
}

Population::~Population () {
	if(xpoints != NULL) {delete [] xpoints;}
	if(ypoints != NULL) {
		for(int ir = 0; ir < numr; ir++) {
			if(ypoints[ir] != NULL) {
				delete [] ypoints[ir];
			}
		}
		delete [] ypoints;
	}
}
  
void Population::readDynamics (const string& fn) {
    
	/*
	 * READ INDUCER
	 */
    ifstream infile;
    infile.open (fn1.c_str());
    if (!infile) {
        cerr << "Error: unable to open input file: " << infile << endl;
		exit(1);
    }

	//	number of xpoints, inducers and proteins
	infile >>  numx >> numind >> numprot;
    if (infile.bad ()) throw runtime_error ("IO stream corrupted");
	if (infile.fail ()) throw runtime_error ("bad data");
	if (!numx) {
		cerr << "Error: empty data" << endl;
		exit(1);
	}
	
	//	number of columns, inducer plus protein
	int numy = numind + numprot;

	//	read dynamic data
	xpoints = new double [numx];
	ypoints = new double*[numx];
	for(int ir = 0; ir < numx; ir++) {
		infile >> xpoints[ir];
		ypoints[ir] = new double*[numy];
		for(int ic = 0; ic < numy; ic++) {
			infile >> ypoints[ir][ic];
			if (infile.bad ()) throw runtime_error ("IO stream corrupted");
			if (infile.fail ()) throw runtime_error ("bad data");
		}
	}

	//	foreach cell, initialization
	for (int ic=0; ic < ncell; ic++) {
		cells[ic] = new Cell (numind, numprot);
	}

	return;
}

