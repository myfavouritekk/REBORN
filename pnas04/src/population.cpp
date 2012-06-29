#include "population.h"



//constructor, by default it will use score function 1 and evlute 100 generations
Population::Population (const int& _ncell)
:ncell(_ncell), numr(0), numind(0), numprot(0), xpoints(NULL), ypoints(NULL), cells(NULL) ,sfunc(ScoreFunc(1)),evolution(100)
{
	ncell = ncell <=0 ? 100 : ncell;
}

//destructor
Population::~Population () {
	if(xpoints != NULL) {delete [] xpoints;}
	if(ypoints != NULL) {
		for(int ir = 0; ir < numind + numprot; ir++) {
			if(ypoints[ir] != NULL) {
				delete [] (ypoints[ir]);
			}
		}
		delete [] ypoints;
	}
    
    
    //delete set of cells
    if (cells != NULL) {
        //ncell ~ ncell*2-1 have been deleted in selection()
        for (int i = 0; i < ncell; i++) {
            if (cells[i] != NULL) {
                delete cells[i];
            }
        }
        delete cells;
    }
}


//population initializer
void Population::init(){
    
    //read dynamic data
    std::string fn;
    std::cout << "Enter the input file name:" << std::endl;
    std::cin >> fn;
    
    //initialize the total set of cells that can contain 2*ncell cells
    cells = new Cell*[2 * ncell];
    
    readDynamics(fn);
    
    
   
}



//for growth phase
void Population::growth(){;
    Cell* currCell;
    for(int i = 0; i < ncell; i++){
        
        //a cell will duplicate itself
        Cell* aCell = new Cell(*(cells[i]));
        cells[ncell + i] = aCell;
        
        //mutation phase
        currCell = cells[ncell + i];
        currCell->mutation();
        
        //get its score
        currCell->getScore(sfunc, ypoints, numind + numprot, numr, false);
    }
    evolution--;//evolution once
}


void quickSort(Cell* cells[],int num){
    
    //end of recursion, no need to sort
    if(num <= 1)return;
    
    int numLess = 0, numGreater = 0;//to store lengths of two subgroups
    Cell** less = new Cell*[num];
    Cell** greater = new Cell*[num];
    
    
    for(int i = 1; i < num; i++){//starts from i = 1, the second member
        if(cells[i]->getCurrScore() > cells[0]->getCurrScore()){
            greater[numGreater] = cells[i];
            numGreater++;
        }else{
            less[numLess] = cells[i];
            numLess++;
        }
    }
    
    
    //recursion
    quickSort(less,numLess);
    quickSort(greater,numGreater);
    
    //final assignments
    cells[numLess] = cells[0];
    for(int i = 0; i < numLess; i++){
        cells[i] = less[i];
    }
    for(int i = 0; i + numLess + 1 < num; i++){
        cells[i+numLess+1]=greater[i];
    }
    
    delete [] less;
    delete [] greater;
}





//for seclction phase
void Population::selection(){
    quickSort(cells,ncell*2);
    
    //cells with lower score will extinct
    for (int i = ncell; i < 2 * ncell; i++) {
        if (cells[i] != NULL) {
            delete cells[i];
        }
    }
    
    std::cout << "Finished Evolution: " << 100 - evolution << std::endl;
    std::cout << "BestScore: " << cells[0]->getCurrScore() << std::endl;
}



bool Population::isTerminate(){
    return !evolution;//if evolution equals 0, then evolution terminates
}



//reading dynamic data set
using namespace std;
void Population::readDynamics (const string& fn) {
	/*
	 * READ INDUCER
	 */
    ifstream infile;
    infile.open (fn.c_str());
    if (!infile) {
        cerr << "Error: unable to open input file: " << infile << endl;
		exit(1);
    }

	//	number of xpoints, inducers and proteins
	infile >>  numr >> numind >> numprot;
    cout << numr << endl << numind << endl << numprot << endl;
    if (infile.bad ()) throw runtime_error ("IO stream corrupted");
	if (infile.fail ()) throw runtime_error ("bad data");
	if (!numr) {
		cerr << "Error: empty data" << endl;
		exit(1);
	}
	
	//	number of rows, inducer plus protein
	int numy = numind + numprot;

	//	read dynamic data
	xpoints = new double [numr];
	ypoints = new double*[numy];
    for (int i = 0; i < numy; i++) {
        ypoints[i] = new double[numr];
    }
	for(int ir = 0; ir < numr; ir++) {
		infile >> xpoints[ir];
        cout << xpoints[ir] << "\t";
		for(int ic = 0; ic < numy; ic++) {
			infile >> ypoints[ic][ir];
            cout << ypoints[ic][ir] << "\t";
			if (infile.bad ()) throw runtime_error ("IO stream corrupted");
			if (infile.fail ()) throw runtime_error ("bad data");
		}
        cout << endl;
	}

	//	for each cell, initialization
	for (int i = 0; i < ncell; i++) {
		cells[i] = new Cell(numind, numprot);
        cells[i]->getScore(sfunc, ypoints, numind + numprot, numr, false);//getScore in initialization
	}

	return;
}

void Population::readDynamicsFromConsole(){
    
}


void Population::genTikzFormat(){
    //to be implemented
    return;
}

void Population::output(){
    
    Cell* currCell;
    for (int i = 0; i < ncell; i++) {
        cout << "Cell " << i+1 << endl;
        currCell = cells[i];
        currCell->description();
    }
    
    //plot the best result
    currCell = cells[0];
    currCell->getScore(sfunc, ypoints, numind + numprot, numr, true);
    
    
}

