#include "population.h"



//constructor, by default it will use score function 1 and evlute 100 generations
Population::Population (const int& _ncell, const int& _evolution)
:ncell(_ncell), numr(0), numind(0), numprot(0), xpoints(NULL), ypoints(NULL), cells(NULL) ,sfunc(ScoreFunc(1)),evolution(_evolution)
{
	ncell = ncell <=0 ? 100 : ncell;
    evolution = evolution <= 0? 1000 : evolution;
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

int Population::getEvolution(){
    return evolution;
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
        currCell->mutation();
        currCell->mutation();
        
        //get its score
        currCell->getScore(sfunc, ypoints, numind + numprot, numr, false);
    }
    evolution--;//evolution once
}

//mutation only for topology
void Population::mutation(){
    Cell* currCell;
    for(int i = 10; i < ncell; i++){
        
        currCell = cells[i];
        
        //mutation topology
        currCell->mut_topology();
        
    }
    evolution--;//evolution once
}

//mutation only for kinetics
void Population::mut_parameters(){
    Cell* curCell;
    for (int i = 0; i < ncell; i++) {
        Cell* aCell = new Cell(*(cells[i]));
        cells[ncell + i] = aCell;
        
        curCell = cells[ncell + i];
        curCell -> mut_parameters();
        
        curCell -> getScore(sfunc, ypoints, numind + numprot, numr, false);
        cells[i] -> getScore(sfunc, ypoints, numind + numprot, numr, false);
        
        if (curCell -> getCurrScore() < cells[i] -> getCurrScore()) {
            delete cells[i];
            cells[i] = curCell;
        }else {
            delete curCell;
        }
    }
    evolution--;
}

//mutation only for kinetics using simulated annealing algorithm
void Population::mut_parameters_simAnneal(){
    Cell* curCell;
    for (int i = 0; i < ncell; i++) {
        Cell* aCell = new Cell(*(cells[i]));
        cells[ncell + i] = aCell;
        
        curCell = cells[ncell + i];
        curCell -> mut_parameters_simAnneal();
        
        curCell -> getScore(sfunc, ypoints, numind + numprot, numr, false);
        cells[i] -> getScore(sfunc, ypoints, numind + numprot, numr, false);
        
        if (curCell -> getCurrScore() < cells[i] -> getCurrScore()) {
            delete cells[i];
            cells[i] = curCell;
        }else {
            delete curCell;
        }
    }
    evolution--;
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


//sort after changing kinects for several generations
void Population::sort(){
    quickSort(cells, ncell);
    for (int i = 0; i < ncell; i++) {
        cells[i]->addRanking(i + 1); //add ranking i to the ith cell
    }
}


//for seclction phase
void Population::selection(){
    quickSort(cells,ncell*2);
    
    for (int i = 0; i < ncell; i++) {
        cells[i]->addRanking(i + 1); //add ranking i to the ith cell
    }
    
    //cells with lower score will extinct
    for (int i = ncell; i < 2 * ncell; i++) {
        if (cells[i] != NULL) {
            delete cells[i];
        }
    }
    
    std::cout << "Finished Evolution: " << TOTAL_EVO - evolution + 1 << std::endl;
    std::cout << "BestScore: " << cells[0]->getCurrScore() << std::endl;
}




bool Population::isTerminate(){
    return (evolution <= 0) ? true : false ;//if evolution equals 0, then evolution terminates
}



//reading dynamic data set
void Population::readDynamics (const string& fn) {
	/*
	 * READ INDUCER
	 */
    std::ifstream infile;
    infile.open (fn.c_str());
    if (!infile) {
        std::cerr << "Error: unable to open input file: " << infile << std::endl;
		exit(1);
    }

	//	number of xpoints, inducers and proteins
	infile >>  numr >> numind >> numprot;
    std::cout << numr << std::endl << numind << std::endl << numprot << std::endl;
    if (infile.bad ()) throw std::runtime_error ("IO stream corrupted");
	if (infile.fail ()) throw std::runtime_error ("bad data");
	if (!numr) {
		std::cerr << "Error: empty data" << std::endl;
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
        std::cout << xpoints[ir] << "\t";
		for(int ic = 0; ic < numy; ic++) {
			infile >> ypoints[ic][ir];
            std::cout << ypoints[ic][ir] << "\t";
			if (infile.bad ()) throw std::runtime_error ("IO stream corrupted");
			if (infile.fail ()) throw std::runtime_error ("bad data");
		}
        std::cout << std::endl;
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
/*Population::genXMLFormat()
 *purpose: generate .xml file to store all the information of current Population
 *details: the xml output reads like this:
 *<Survivals>
    <Cell>
        <Index> cell index </Index>
        <Nodes>
            <Node>
                <Index> node index </Index>
                <NodeString> node string </NodeString>
            </Node>
        </Nodes>
        <Reaction List>
            <Reaction>
                <Type> reaction type </Type>
                <Reactants>
                    <Node> node string </Node>
                </Reactants>
                <Modifiers>
                    <Node> node string </Node>
                </Modifiers>
                <Products>
                    <Node> node string </Node>
                </Products>
                <Forward Rate> forward rate </Forward Rate>
                <Reverse Rate> reverse rate </Reverse Rate>
            </Reaction>
        </Reaction List>
    </Cell>
  </Survivals>
 *
 */
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
    xmlSaveFormatFileEnc("population.xml", outputXML, NULL, 1);
    xmlFreeDoc(outputXML);
}

//after evolution, classify those cells
void Population::classification(){
    for (int i = 0; i < ncell; i++) {
        Cell* currCell = cells[i];        
        if (classifiedCells.size() == 0) {
            std::vector<Cell*>* firstType = new std::vector<Cell*>;
            classifiedCells.push_back(*firstType);
        }
        std::vector<std::vector<Cell*> >::iterator iter_out = classifiedCells.begin();

        while (iter_out != classifiedCells.end()) {
            if ((*iter_out).size() == 0) {
                (*iter_out).push_back(currCell);
                break;
            }
            if (*(*iter_out)[0] == *currCell) {
                (*iter_out).push_back(currCell);
                break;
            }
            if (iter_out + 1 == classifiedCells.end()) {
                std::vector<Cell*>* newType = new std::vector<Cell*>;
                newType->push_back(currCell);
                classifiedCells.push_back(*newType);
                break;
            }
            iter_out++;
        }
    }
    
    //generate XML file containing classified cells
    xmlDocPtr classifiedCellsXMLOutput = xmlNewDoc(BAD_CAST"1.0");
    xmlNodePtr root_node = xmlNewNode(NULL, BAD_CAST"Classified Cells");
    xmlDocSetRootElement(classifiedCellsXMLOutput, root_node);
    int type = 1;
    std::vector<std::vector<Cell*> >::iterator iter_out = classifiedCells.begin();
    while (iter_out != classifiedCells.end()) {
        
        //constructing <Classified Cells>-<Type>
        char classificationType[10];
        sprintf(classificationType, "%d",type);
        xmlNodePtr aType = xmlNewNode(NULL, BAD_CAST"Type");
        xmlAddChild(root_node, aType);
        //constructing <Classified Cells>-<Type>-<Type Index>
        xmlNewTextChild(aType, NULL, BAD_CAST"Type Index", BAD_CAST classificationType);
        //constructing <Classified Cells>-<Type>-<Cells>
        xmlNodePtr cellsOfThisType = xmlNewNode(NULL, BAD_CAST"Cells");
        xmlAddChild(aType, cellsOfThisType);

        //constructing <Classified Cells>-<Type>-<Cells>-<Cell>
        std::vector<Cell*>::iterator iter_inner = (*iter_out).begin();
        while (iter_inner != iter_out->end()) {
            xmlNodePtr aCell = xmlNewNode(NULL, BAD_CAST"Cell");
            xmlAddChild(cellsOfThisType, aCell);
            xmlNodePtr nodes = xmlNewNode(NULL, BAD_CAST"Nodes");
            xmlAddChild(aCell, nodes);
            
            //adding every Node in nodes vector
            std::vector<Node*>::iterator iter_node = (*iter_inner)->getNodesVector()->begin();
            std::vector<Node*>::iterator iter_node_end = (*iter_inner)->getNodesVector()->end();
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
            std::vector<Reaction*>::iterator iter_reaction = (*iter_inner)->getRlistVector()->begin();
            std::vector<Reaction*>::iterator iter_reaction_end = (*iter_inner)->getRlistVector()->end();
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
            iter_inner++;
        }
        iter_out++;
        type++;
    }
    
    xmlKeepBlanksDefault(0);
    xmlSaveFormatFileEnc("ClassifiedCells.xml", classifiedCellsXMLOutput, NULL, 1);
    xmlFreeDoc(classifiedCellsXMLOutput);
}

void Population::genTikzFormat(){
    //to be implemented
    return;
}

void Population::output(){
    
    Cell* currCell;
    for (int i = 0; i < ncell; i++) {
        std::cout << "Cell " << i+1 << std::endl;
        currCell = cells[i];
        currCell->generateTimeCourses(ypoints, numind + numprot, numr);
        
        int generation = TOTAL_EVO - evolution + 1;// current generation of whole evolution process
        int ranking = i + 1; //ranking of this cell
        std::stringstream name;
        name << "generation_" << generation << "_cell_" << ranking << ".txt";
        currCell -> printCurrDataToAFile(name.str(), numr);
        
        currCell->description(numr);
        
    }
    
    //plot the best result
    currCell = cells[0];
    std::cout << "Best Time Courses:" << std::endl;
    currCell->getScore(sfunc, ypoints, numind + numprot, numr, true);
    std::cout << "Best Cell:" << std::endl;
    currCell->generateTimeCourses(ypoints, numind + numprot, numr);
    currCell->description(numr);
    
    
}

