#include "cell.h"
#include "reaction.h"
#include <time.h>
#include <vector>
#include <iostream>


//constructor
Cell::Cell(const int& _numind, const int& _numprot) {
	int currIndex = nodes.size();
	for(int im = 0; im < _numind; im++) {
		Node* inducer = new Node(currIndex, 1);
        nodes.push_back(inducer);
        currIndex++;
	}
	for(int ip = 0; ip < _numprot; ip++) {
		Node* gene = new Node(currIndex, 2);//Type 2 is gene
        nodes.push_back(gene);
        currIndex++;
		Node* prot = new Node(currIndex, 3);//type 3 is protein
        nodes.push_back(prot);
        currIndex++;
	}
}


//destructor
Cell::~Cell() {}



bool Cell::existsNode (const Node& node) {
	vector<Node*>::iterator iter = nodes.begin();
	while (iter != nodes.end()) {
		if(*(*iter) == node) {return true;}
		iter ++;
	}
	return false;
}

bool Cell::existsReaction (const Reaction& rxn) {
	vector<Reaction*>::iterator iter = rlist.begin();
	while (iter != rlist.end()) {
		if(*(*iter) == rxn) {return true;}
		iter++;
	}
	return false;
}

void Cell::mut_deg_prot () {

	if(!rlist.size()) return;

    //  the degradation rate of a protein is modified

    int num1 = 0;//contain number of degradation reactions
    int num2 = 0;//index of a certain reaction
    vector<int> indice;

    //  scan <rlist> and choose one degradation reaction (protein)
    std::vector<Reaction*>::iterator iter = rlist.begin();
    while (iter != rlist.end()) {
        if((*iter)->getRtype() == 1) { //   #1 is protein degradation
            num1++;
            indice.push_back(num2);
        }
        num2++;
        iter ++;
    }
	
	if(!num1) return;	//	no protein

    srand(time(NULL));
    int opIndex = indice[rand()%num1];//opIndex contains a certain degradation reaction

    //  modify its degradation rate
    Reaction* currR = rlist[opIndex];
    currR->modifyForwardRate ();

    return;
}

void Cell::mut_kin_const () {

	if(!rlist.size()) return;
	
	//	a kinetic constant of one reaction is modified
	
	srand(time(NULL));
	Reaction* currR = rlist[rand()%rlist.size()];

	if((double)rand()/RAND_MAX<=0.5 || !currR->isReversible()) {
		currR->modifyForwardRate();
	}
	else {
		currR->modifyReverseRate();
	}

	return;
}

void Cell::mut_add_gene () {//	add a gene
	
	//	create new nodes representing this gene and its protein
	int currNI = nodes.size();

	Node* gene = new Node(currNI,2);	//	gene
	nodes.push_back(gene);

	Node* prot = new Node(currNI+1,3);	//	prot
	nodes.push_back(prot);
	
	//	create reaction 0, transcription
	Reaction* r0 = new Reaction ();
	r0->setReversible(false);
	r0->initForwardRateRandomly();
	r0->addModifier(gene);
	r0->addProduct(prot);
	rlist.push_back(r0);

	//	create reaction 1, protein degradation
	Reaction* r1 = new Reaction ();
	r1->setReversible(false);
	r1->initForwardRateRandomly();
	r1->addReactant(prot);
	rlist.push_back(r1);

	return;
}

void Cell::mut_add_regu () {
	
	//A new interaction between a protein and a gene or a gene/protein complex is introduced
    int num1 = 0;//number of protein
    int num2 = 0;//number of gene/protein complex
	int num3 = 0;//index
    vector<int> protIndice;
	vector<int> cplxIndice;

    //  scan <nodes> and choose  gene and one gene/promoter complex
    std::vector<Node*>::iterator iter1 = nodes.begin();
    while (iter1 != nodes.end()) {
        if((*iter1)->getNtype() == 3) { //   #3 is a protein
            num1++;
            protIndice.push_back(num3);
        }
		
		if((*iter1)->getNtype() == 5) {
			num2++;
			cplxIndice.push_back(num3);
		}

        num3++;
        iter1++;
    }

	if(!num1 || !num2) return; // no gene or gene/protein complex

	//	sample random number
    srand(time(NULL));
    int opIndex1 = protIndice[rand()%num1];//contain index of a protein
    int opIndex2 = cplxIndice[rand()%num2];//contain index of a gene or a gene/protein complex
	Node* exGene = nodes[opIndex2]->extractFirstGene();
	if(exGene == NULL) return;
	
	//	find the index of the protein created by this extracted gene
	Node* exProt = NULL;
    std::vector<Reaction*>::iterator iter2 = rlist.begin();
    while (iter2 != rlist.end()) {
        if((*iter2)->getRtype() == 0) { //   #0 is gene transcription
			Node* sr = (*iter2)->getModifier(0);
			if(sr != NULL && (sr->getNindex() == exGene->getNindex())) {
				exProt = (*iter2)->getProduct(0);
				break;
			}
        }
        iter2 ++;
    }
	if(exProt == NULL) return;

	//	create a new node as the product
	Node* ncomplex = new Node(nodes.size(), 5, nodes[opIndex1], nodes[opIndex2]);
	nodes.push_back(ncomplex);

	//	create reaction 0, transcription
	Reaction* r0 = new Reaction ();
	r0->setReversible(false);
	r0->initForwardRateRandomly();
	r0->addModifier(ncomplex);
	r0->addProduct(exProt);
	rlist.push_back(r0);

	//	create reaction 1, binding/unbinding between protein and gene or gene/protein complex
	Reaction* r1 = new Reaction ();
	r1->setReversible(true);
	r1->initForwardRateRandomly();
	r1->initReverseRateRandomly();
	r1->addReactant(nodes[opIndex1]);
	r1->addReactant(nodes[opIndex2]);
	r1->addProduct (ncomplex);
	rlist.push_back(r1);

	return;
}

void Cell::mut_add_postmod () {
	//	a post modification is add
	
    srand(time(NULL));
	if((double)rand()/RAND_MAX <= 0.5) {//	a single protein case


	}
	else {//	two proteins (protein complex) case
		
		int num1 = 0;
		int num2 = 0;
		vector<int> protIndice;
			
		vector<Node*>::iterator iter = nodes.begin();
		while (iter != nodes.end()) {
			if((*iter)->getNode(0) == NULL) { 
				num1++;
				protIndice.push_back(num2);
			}
			num2++;
			iter ++;
		}

		if(!num1) return;	//	no protein or protein complex
		
		int opIndex1 = protIndice[rand()%num1];	//	protein 1
		int opIndex2 = protIndice[rand()%num1];	//	protein 2

		double possibility = (double)rand()/RAND_MAX;
		if(possibility < 1/3) {//	dimerization
			Node* dimer = new Node (nodes.size(), nodes[opIndex1], nodes[opIndex2]);
			Reaction* dimerization = new Reaction ();
			dimerization->setReversible(true);
			dimerization->initForwardRateRandomly();
			dimerization->initReverseRateRandomly();
			dimerization->addReactant(nodes[opIndex1]);
			dimerization->addReactant(nodes[opIndex2]);
			dimerization->addProduct(dimer);
			if(!existsReaction(*dimerization)) {
				nodes.push_back(dimer);
				rlist.push_back(dimerization);
			}
			else {
				delete dimer;
				delete dimerization;
			}
		}
		else if(possibility < 2/3) {
			Reaction* partialDeg = new Reaction ();
			partialDeg->setReversible(false);
			partialDeg->initForwardRateRandomly();
		}
		else {
            
            
            //to be continued
		}

	}
	
	return;
}

#define PROB1 0.5
#define PROB2 0.5
#define PROB3 0.5
#define PROB4 0.5
#define PROB5 0.0


//overall mutation method
void Cell::mutation(){
    srand(time(NULL));
    if (rand() < RAND_MAX*PROB1) {
        mut_deg_prot();
    }
    if (rand() < RAND_MAX*PROB2) {
        mut_kin_const();
    }
    if (rand() < RAND_MAX*PROB3) {
        mut_add_gene();
    }
    if (rand() < RAND_MAX*PROB4) {
        mut_add_regu();
    }
    if (rand() < RAND_MAX*PROB5) {
        mut_add_postmod();
    }
}

//runge-kutta method
#define MAXSTEPS 100
void runge_kutta(double data[][MAXSTEPS],double (*odes[])(double y[],double x),int series, int steps)
{//100 is steps
	
	int currSerie,currStep;
	int i;
    double k1,k2,k3,k4;
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
    delete y, tempY;
}




/*get score using the sfunc as score function and change its own currScore member
 *fist: using Runge-Kutta method to generate the time course and store all them in currData
 *second: using the sfunction as a score function, and passing currData and targetData as parameters to calculate score
 *third: assign the score to currScore
 */
void Cell::getScore(const ScoreFunc& sfunc, double** targetData){
    /*implement code here to calculate the score
    *the biggest problem here is how to turn reaction into methods
    *
    */ 
    return;
}

	
//return the currSore generated by getScore method
double Cell::getCurrScore(){
    std::cout << currScore << std::endl;
    return currScore;
}




//duplicate itself
Cell* Cell::aNewCopy(){
    Cell *newCell = new Cell(0,0);
    *newCell = *this;
    return newCell;
}














