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
        delete [] cells;
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

#define BAD_CAST (xmlChar *)
void Population::genXMLFormat(){
    xmlDocPtr outputXML = xmlNewDoc(BAD_CAST"1.0");
    xmlNodePtr root_node = xmlNewNode(NULL, BAD_CAST"Survivals");
    xmlDocSetRootElement(outputXML, root_node);
    for (int i = 0; i < ncell; i++) {
        char cellIndex[10];
        sprintf(cellIndex, "%d",i + 1);
        xmlNodePtr aCell = xmlNewNode(NULL, BAD_CAST"Cell");
        xmlAddChild(root_node, aCell);
        xmlNewTextChild(aCell, NULL, BAD_CAST"Index", BAD_CAST cellIndex);
        xmlNodePtr nodes = xmlNewNode(NULL, BAD_CAST"Nodes");
        xmlAddChild(aCell, nodes);
        
        //adding every Node in nodes vector
        std::vector<Node*>::iterator iter_node = cells[i]->getNodesVector()->begin();
        std::vector<Node*>::iterator iter_node_end = cells[i]->getNodesVector()->end();
        while (iter_node != iter_node_end) {
            xmlNodePtr aNode = xmlNewNode(NULL, BAD_CAST"Node");
            xmlAddChild(nodes, aNode);
            char nodeIndex[10];
            sprintf(nodeIndex, "%d",(*iter_node)->getNindex());
            xmlNewTextChild(aNode, NULL, BAD_CAST"Index", BAD_CAST(nodeIndex));
            string nodeString = (*iter_node)->getNstring();
            xmlNewTextChild(aNode, NULL, BAD_CAST"NodeString", BAD_CAST((unsigned char*)nodeString.c_str()));
            iter_node++;
        }
        
        xmlNodePtr rlist = xmlNewNode(NULL, BAD_CAST"Reaction List");
        xmlAddChild(aCell, rlist);
        
        //adding every Reaction in rlist vector
        std::vector<Reaction*>::iterator iter_reaction = cells[i]->getRlistVector()->begin();
        std::vector<Reaction*>::iterator iter_reaction_end = cells[i]->getRlistVector()->end();
        while (iter_reaction != iter_reaction_end) {
            xmlNodePtr aReaction = xmlNewNode(NULL, BAD_CAST"Reaction");
            xmlAddChild(rlist, aReaction);
            char reactionType[10];
            sprintf(reactionType, "%d",(*iter_reaction)->getRtype());
            xmlNewTextChild(aReaction, NULL, BAD_CAST"Type", BAD_CAST reactionType);
            //reactants
            xmlNodePtr reactants = xmlNewNode(NULL, BAD_CAST"Reactants");
            xmlAddChild(aReaction, reactants);
            iter_node = (*iter_reaction)->getReactantsVector()->begin();
            iter_node_end = (*iter_reaction)->getReactantsVector()->end();
            while (iter_node != iter_node_end) {
                std::string nodeString = (*iter_node)->getNstring();
                xmlNewTextChild(reactants, NULL, BAD_CAST"Node", BAD_CAST((unsigned char*)nodeString.c_str()));
                iter_node++;
            }
            //modifiers
            xmlNodePtr modifiers = xmlNewNode(NULL, BAD_CAST"Modifiers");
            xmlAddChild(aReaction, modifiers);
            iter_node = (*iter_reaction)->getModifiersVector()->begin();
            iter_node_end = (*iter_reaction)->getModifiersVector()->end();
            while (iter_node != iter_node_end) {
                std::string nodeString = (*iter_node)->getNstring();
                xmlNewTextChild(modifiers, NULL, BAD_CAST"Node", BAD_CAST((unsigned char*)nodeString.c_str()));
                iter_node++;
            }
            //products
            xmlNodePtr products = xmlNewNode(NULL, BAD_CAST"Products");
            xmlAddChild(aReaction, products);
            iter_node = (*iter_reaction)->getProductsVector()->begin();
            iter_node_end = (*iter_reaction)->getProductsVector()->end();
            while (iter_node != iter_node_end) {
                std::string nodeString = (*iter_node)->getNstring();
                xmlNewTextChild(products, NULL, BAD_CAST"Node", BAD_CAST((unsigned char*)nodeString.c_str()));
                iter_node++;
            }
            
            char forwardRate[20];
            sprintf(forwardRate, "%.6f",(*iter_reaction)->getForwardRate());
            xmlNewTextChild(aReaction, NULL, BAD_CAST"Forward Rate", BAD_CAST forwardRate);
            char reverseRate[20];
            sprintf(reverseRate, "%.6f",(*iter_reaction)->getReverseRate());
            xmlNewTextChild(aReaction, NULL, BAD_CAST"Reverse Rate", BAD_CAST reverseRate);
            
            iter_reaction++;
        }
    }
    xmlKeepBlanksDefault(0);
    xmlSaveFormatFileEnc("XMLOutput.xml", outputXML, NULL, 1);
    xmlFreeDoc(outputXML);
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
        currCell->generateTimeCourses(ypoints, numind + numprot, numr);
        currCell->correlationMatrix(numr);
        currCell->getVariation(numr);
        currCell->fitnessVariation(numr);
        currCell->description(numr);
        
    }
    
    //plot the best result
    currCell = cells[0];
    currCell->getScore(sfunc, ypoints, numind + numprot, numr, true);
    
    
}

