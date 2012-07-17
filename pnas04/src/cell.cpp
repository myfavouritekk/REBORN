#include "cell.h"



/*constructor:
 input: _numind(number of inducers); _numprot(number of proteins)
 */
Cell::Cell(const int& _numind, const int& _numprot):numInducer(_numind) {
	int currIndex = nodes.size();
	int iop = 0;
	int* indexOfProt = new int[_numprot];
	for(int im = 0; im < _numind; im++) {
		Node* inducer = new Node(currIndex, 1);
        nodes.push_back(inducer);
		inputIndice.push_back(currIndex);           //push back inducer's index
       	currIndex++;		
	}
	for(int ip = 0; ip < _numprot; ip++) {
		Node* gene = new Node(currIndex, 2);//Type 2 is gene
        nodes.push_back(gene);
        currIndex++;
		Node* prot = new Node(currIndex, 3);//type 3 is protein
        nodes.push_back(prot);
		inputIndice.push_back(currIndex);             //push back protein's index
        indexOfProt[iop] = currIndex;
		iop++;
		currIndex++;

		Reaction *r0 = new Reaction(0);         //add transcription reaction
		r0->setReversible(false);
		r0->initForwardRateRandomly();
		r0->addModifier(gene);
		r0->addProduct(prot);

		Reaction* r1 = new Reaction(1);        //add degradation reaction
		r1->setReversible(false);
		r1->initForwardRateRandomly();
		r1->addReactant(prot);

		rlist.push_back(r0);
		rlist.push_back(r1);
	}
	for(int ioi =  0; ioi < _numind; ioi++){    //adding binding between an inducer and a protein
		int iopIndex = indexOfProt[rand()% _numprot];
		Node* selectedProt = nodes[iopIndex];
		Node* selectedInd = nodes[ioi];
		Node* complex = new Node(currIndex, 4, selectedProt, selectedInd);
		nodes.push_back(complex);
        
		Reaction* r8 = new Reaction(8); //reaction type 8 is binding between an inducer and a protein
		r8->setReversible(true);
		r8->setForwardRate (40);        //setting forward rate to a relatively large number to ensure the binding
		r8->initReverseRateRandomly();
		r8->addReactant(selectedProt);
		r8->addReactant(selectedInd);
		r8->addProduct(complex);
		rlist.push_back(r8);
	}
	delete [] indexOfProt;
}


//duplicate itself
Cell::Cell(Cell &cell){
    
    numInducer = cell.numInducer;
    
    //copy every node in the cell
    std::vector<Node*>::iterator iter_node = cell.nodes.begin();
    std::vector<Node*>::iterator iter_node_end = cell.nodes.end();
    while (iter_node != iter_node_end) {
        Node* nodeCopy = new Node(*(*iter_node));
        nodes.push_back(nodeCopy);
        iter_node++;
    }
    
    /*constructing the relationships of these cells, 
     *that is implementing components vector in each cell
     */
    int nodesSize = cell.nodes.size();
    for (int indexNinCell = 0; indexNinCell < nodesSize; indexNinCell++) {
        int componentsSize = cell.nodes[indexNinCell]->getNsize();
        for (int componentIndex = 0; componentIndex < componentsSize; componentIndex++) {
            if (componentIndex == 0 && !cell.nodes[indexNinCell]->getNode(componentIndex)) {//for no gene case, components[0]==NULL
                this->nodes[indexNinCell]->pushNode(NULL);
            }else {
                int aComponentIndex = cell.nodes[indexNinCell]->getNode(componentIndex)->getNindex();//store the index of the ith components of this node in the cell
                Node* correspondingComponent = this->nodes[aComponentIndex];
                this->nodes[indexNinCell]->pushNode(correspondingComponent);
            }
            
        }
    }
    
    
    
    //copy every reaction in the cell
    std::vector<Reaction*>::iterator iter_reaction = cell.rlist.begin();
    std::vector<Reaction*>::iterator iter_reaction_end = cell.rlist.end();
    while (iter_reaction != iter_reaction_end) {
        Reaction* reactionCopy = new Reaction(*(*iter_reaction));
        /*reconstructing the relationships of nodes and rlist,
         *based on the those relationships in the old cell called "cell"
         */
        int reactantsSize = (*iter_reaction)->getReactantsSize();
        for (int indexInReactants = 0; indexInReactants < reactantsSize; indexInReactants++) {
            int aNodeIndex = (*iter_reaction)->getReactant(indexInReactants)->getNindex();
            Node* coorespondingNode = this->nodes[aNodeIndex];
            reactionCopy->addReactant(coorespondingNode);
        }
        int modifiersSize = (*iter_reaction)->getModifiersSize();
        for (int indexInModifiers = 0; indexInModifiers < modifiersSize; indexInModifiers++) {
            int aNodeIndex = (*iter_reaction)->getModifier(indexInModifiers)->getNindex();
            Node* coorespondingNode = this->nodes[aNodeIndex];
            reactionCopy->addModifier(coorespondingNode);
        }
        int productsSize = (*iter_reaction)->getProductsSize();
        for (int indexInProducts = 0; indexInProducts < productsSize; indexInProducts++) {
            int aNodeIndex = (*iter_reaction)->getProduct(indexInProducts)->getNindex();
            Node* coorespondingNode = this->nodes[aNodeIndex];
            reactionCopy->addProduct(coorespondingNode);
        }
        
        rlist.push_back(reactionCopy);
        iter_reaction++;
    }
    
    inputIndice = cell.inputIndice;
    
}







//destructor
Cell::~Cell() {
    std::vector<Node*>::iterator iter_node = nodes.begin();
    std::vector<Node*>::iterator iter_node_end = nodes.end();
    while (iter_node != iter_node_end) {
        delete (*iter_node);
        iter_node++;
    }
    
    std::vector<Reaction*>::iterator iter_reaction = rlist.begin();
    std::vector<Reaction*>::iterator iter_reaction_end = rlist.end();
    while (iter_reaction != iter_reaction_end) {
        delete (*iter_reaction);
        iter_reaction++;
    }
    
}



bool Cell::existsNode (const Node& aNode) {
	vector<Node*>::iterator iter = nodes.begin();
    vector<Node*>::iterator iter_end = nodes.end();
	while (iter != iter_end) {
		if(*(*iter) == aNode) {return true;}
		iter ++;
	}
	return false;
}

bool Cell::existsReaction (const Reaction& aReaction) {
	vector<Reaction*>::iterator iter = rlist.begin();
    vector<Reaction*>::iterator iter_end = rlist.end();
	while (iter != iter_end) {
		if(*(*iter) == aReaction) {return true;}
		iter++;
	}
	return false;
}



bool Cell::operator==(Cell& aCell){
    if (nodes.size() != aCell.nodes.size() || rlist.size() != aCell.rlist.size()) {
        return false;
    }else {
        std::vector<Node*>::iterator iter_node = aCell.nodes.begin();
        std::vector<Node*>::iterator iter_node_end = aCell.nodes.end();
        while (iter_node != iter_node_end) {
            if (!existsNode(*(*iter_node))) {
                return false;
            }
            iter_node++;
        }
        std::vector<Reaction*>::iterator iter_reaction = aCell.rlist.begin();
        std::vector<Reaction*>::iterator iter_reaction_end = aCell.rlist.end();
        while (iter_reaction != iter_reaction_end) {
            if (!existsReaction(*(*iter_reaction))) {
                return false;
            }
            iter_reaction++;
        }
    }
    return true;
}

/*
 *five types of mutation:
 *1. deg_prot: change protein degradation 
 *2. kin_const: change kinetic constant
 *3. add_gene: create a new gene
 *4. add_regu: add new interaction between protein and a gene or a gene-protein complex
 *5. add_postmod: add post transcriptional modification
 */

//mut_deg_prot: change the protein degradation rate
void Cell::mut_deg_prot () {

	if(!rlist.size()) return;               //no reaction in the rlist, no modification

    //  the degradation rate of a protein to be modified
    int numOfDegReaction = 0;               //contain number of degradation reactions
    int indexOfDegReaction = 0;             //index of a certain reaction
    vector<int> indice;                     //contain the indice of degradation reactions

    //  scan <rlist> and choose one degradation reaction (protein)
    std::vector<Reaction*>::iterator iter = rlist.begin();
    std::vector<Reaction*>::iterator iter_end = rlist.end();
    while (iter != iter_end) {
        if((*iter)->getRtype() == 1) {      //#1 is protein degradation
            numOfDegReaction++;
            indice.push_back(indexOfDegReaction);
        }
        indexOfDegReaction++;
        iter ++;
    }
	
	if(!numOfDegReaction) return;           //no protein degradation reaction

    //srand((unsigned int)time(NULL));
    int opIndex = indice[rand()%numOfDegReaction];//opIndex contains a certain degradation reaction

    //  modify its degradation rate
    Reaction* currR = rlist[opIndex];
    currR->modifyForwardRate ();

    return;
}


//mut_kin_const: change the kinect constant of a reaction
void Cell::mut_kin_const () {

	if(!rlist.size()) return;
	
	//	a kinetic constant of one reaction to be modified
	
	//srand((unsigned int)time(NULL));
	Reaction* currR = rlist[rand()%rlist.size()];
	if((double)rand()/RAND_MAX <= 0.5 || !currR->isReversible()) {
		if(currR -> getRtype() != 8){
           currR->modifyForwardRate();
        }
	}
	else {
		currR->modifyReverseRate();
	}

	return;
}


//mut_add_gene: add new gene and its protein into nodes, and transcription and protein degradation reaction into rlist
void Cell::mut_add_gene () {            //	add a gene
	
	//	create new nodes representing this gene and its protein
	int currNI = nodes.size();

	Node* gene = new Node(currNI,2);	//	gene
	nodes.push_back(gene);
    
	Node* prot = new Node(currNI+1,3);	//	protein
	nodes.push_back(prot);
	
	//	create Reaction r0, transcription
	Reaction* r0 = new Reaction (0);
	r0->setReversible(false);// irreversible and no reactants
	r0->initForwardRateRandomly();
	r0->addModifier(gene);
	r0->addProduct(prot);

	//	create Reaction r1, protein degradation
	Reaction* r1 = new Reaction (1);
	r1->setReversible(false);
	r1->initForwardRateRandomly();
	r1->addReactant(prot);

	if(!existsReaction(*r0) && !existsReaction(*r1)){
		rlist.push_back(r0);
		rlist.push_back(r1);
	}
	else{
		delete r0;
        delete r1;
		return;
	}
		
	double random = (double)rand()/RAND_MAX;

	//add dimerization A+B->A:B
	if(random < 0.3){      

		int numProt=0;
		int index=0;
		std::vector<Node*>::iterator iter1 = nodes.begin();
		std::vector<Node*>::iterator iter1_end = nodes.end();
		vector<int> proIndice;
		while(iter1 != iter1_end){
			if((*iter1)->getNtype() == 3||(*iter1)->getNtype() == 6){
				numProt++;
				proIndice.push_back(index);
			}
			index ++;
			iter1 ++;
		}

		int opIndex=proIndice[rand()%numProt];
		Reaction* dimerization = new Reaction(5);
		Node* dimer = new Node(nodes.size(),prot, nodes[opIndex]);
		nodes.push_back(dimer);
		dimerization->setReversible(false);
		dimerization->initForwardRateRandomly();
		dimerization->addReactant(prot);
		dimerization->addReactant(nodes[opIndex]);
		dimerization->addProduct(dimer);

		Reaction* degrad = new Reaction(1);
		degrad->setReversible(false);
		degrad->initForwardRateRandomly();
		degrad->addReactant(dimer);

		rlist.push_back(degrad);
		rlist.push_back(dimerization);
	}

	//add regulation    A+b->A:b
	if(random >= 0.30 && random < 0.6){ 
		int numGene = 0;
		int index = 0;
		vector<int> geneIndice;

		std::vector<Node*>::iterator iter1 = nodes.begin();
		std::vector<Node*>::iterator iter1_end = nodes.end();
		while(iter1 != iter1_end){
			if((*iter1)->getNtype() == 2||(*iter1)->getNtype() == 5){
				numGene++;
				geneIndice.push_back(index);
			}
			index ++;
			iter1 ++;
		}

		if(!numGene) return;
		int opIndex = geneIndice[rand()%numGene];
		Node* exGene = nodes[opIndex]->extractFirstGene();
		if(exGene == NULL) return;

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

		Node* ncomplex = new Node(nodes.size(), 5, prot, nodes[opIndex]);
		nodes.push_back(ncomplex);

		//	create reaction 0, transcription
		Reaction* r0 = new Reaction (0);
		r0->setReversible(false);
		r0->initForwardRateRandomly();
		r0->addModifier(ncomplex);
		r0->addProduct(exProt);

		//	create reaction 1, binding/unbinding between protein and gene or gene/protein complex
		Reaction* r1 = new Reaction (2);
		r1->setReversible(true);
		r1->initForwardRateRandomly();
		r1->initReverseRateRandomly();
		r1->addReactant(prot);
		r1->addReactant(nodes[opIndex]);
		r1->addProduct (ncomplex);

		rlist.push_back(r0);
		rlist.push_back(r1);
	}


	//add catalytic degradation A+B->A
	if(random >= 0.6 && random <0.9){
		int numOfSingProt=0;
		int index = 0;
		vector<int> indiceOfSingProt;

		std::vector<Node*>::iterator iter1 = nodes.begin();
		std::vector<Node*>::iterator iter1_end = nodes.end();
		while(iter1 != iter1_end){
			if((*iter1)->getNtype() == 3){
				numOfSingProt++;
				indiceOfSingProt.push_back(index);
			}
			index ++;
			iter1 ++;
		}

		int opIndex = indiceOfSingProt[rand()%numOfSingProt];
		Reaction* catalyticDeg = new Reaction(6);
		catalyticDeg->setReversible(false);
		catalyticDeg->initForwardRateRandomly();
		catalyticDeg->addReactant(prot);
		catalyticDeg->addReactant(nodes[opIndex]);
		if(rand()%2 == 0)
			catalyticDeg->addProduct(prot);
		else
			catalyticDeg->addProduct(nodes[opIndex]);

		rlist.push_back(catalyticDeg);
	}


	//add partial catalytic degradation A+BC->B
	if(random >= 0.9){
		int index=0;
		int numOfCompProt=0;
		vector<int> indiceOfCompProt;
		std::vector<Node*>::iterator iter1 = nodes.begin();
		std::vector<Node*>::iterator iter1_end = nodes.end();
		while(iter1 != iter1_end){
			if((*iter1)->getNtype() == 3){
				numOfCompProt++;
				indiceOfCompProt.push_back(index);
			}
			index ++;
			iter1 ++;
		}

		int opIndex=indiceOfCompProt[rand()%numOfCompProt];
		Reaction* partialDeg=new Reaction(7);
		partialDeg->setReversible(false);
		partialDeg->initForwardRateRandomly();
		partialDeg->addReactant(prot);
		partialDeg->addReactant(nodes[opIndex]);
		int randComp=rand()%(nodes[opIndex]->getNsize()-1)+1;
		partialDeg->addProduct(nodes[opIndex]->getNode(randComp));   //choose the protein in the complex

		rlist.push_back(partialDeg);
	}

	return;
}

//mut_add_regu: add interaction
void Cell::mut_add_regu () {
	
	//A new interaction between a protein and a gene or a gene/protein complex is introduced
    int numPro = 0;                         //number of protein
    int numGenePro = 0;                     //number of gene/protein complex
	int index = 0;                          //index
    vector<int> protIndice;
	vector<int> cplxIndice;

    //  scan <nodes> and choose  gene and one gene/promoter complex
    std::vector<Node*>::iterator iter1 = nodes.begin();
    std::vector<Node*>::iterator iter1_end = nodes.end();
    while (iter1 != iter1_end) {
        if((*iter1)->getNtype() == 3) {     //   #3 is a protein
            numPro++;
            protIndice.push_back(index);
        }
		
		if((*iter1)->getNtype() == 5||(*iter1)->getNtype() == 2) {     //   #5 is gene/protein complex  and #2 is gene
			numGenePro++;
			cplxIndice.push_back(index);
		}

        index++;
        iter1++;
    }

	if(!numPro || !numGenePro) return; // no gene or gene/protein complex

	//	sample random number
    //srand((unsigned int)time(NULL));
    int opIndex1 = protIndice[rand()%numPro];//contain index of a protein
    int opIndex2 = cplxIndice[rand()%numGenePro];//contain index of a gene or a gene/protein complex
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
	if(existsNode(*ncomplex)){
		delete ncomplex;
		return;
	}
	nodes.push_back(ncomplex);

	//	create reaction 0, transcription
	Reaction* r0 = new Reaction (0);
	r0->setReversible(false);
	r0->initForwardRateRandomly();
	r0->addModifier(ncomplex);
	r0->addProduct(exProt);

	//	create reaction 1, binding/unbinding between protein and gene or gene/protein complex
	Reaction* r1 = new Reaction (2);
	r1->setReversible(true);
	r1->initForwardRateRandomly();
	r1->initReverseRateRandomly();
	r1->addReactant(nodes[opIndex1]);
	r1->addReactant(nodes[opIndex2]);
	r1->addProduct (ncomplex);

	if(!existsReaction(*r0) && !existsReaction(*r1)){
		rlist.push_back(r0);
		rlist.push_back(r1);
	}
	else {
		delete r0;
		delete r1;
	}

	return;
}



//mut_add_postmod: add a post-transcriptional regulation to the network
void Cell::mut_add_postmod () {
	//	a post modification is add
	
    //srand((unsigned int)time(NULL));
	 /* a protein is to be chosen from the existing ones, and a modified version of it is to be introduced */
	if((double)rand()/RAND_MAX <= 0.5) {                 //	single protein or single protein complex case
        int indexOfProt=0;
		int numOfProt=0;
		vector<int> protIndice;

		vector<Node*>::iterator iter = nodes.begin();
		vector<Node*>::iterator iter_end = nodes.end();
		while (iter !=iter_end){
			if((*iter)->getNtype() == 3 || (*iter)->getNtype() == 6){   //if this node is single protein or single protein complex
			 protIndice.push_back(indexOfProt);
			 numOfProt++;
			}
			indexOfProt++;
			iter++;
		}
		if(!numOfProt) return;   //no single protein
		int opIndex=protIndice[rand()%numOfProt];  //choose a random protein
		if(nodes[opIndex]->getNtype() == 3)       //if it is a single protein, the reaction type is A->A*
		{
			Node* modProt = new Node(nodes.size(),3);
			nodes.push_back(modProt);

			Reaction* r0=new Reaction(3);    //add modification reaction
			r0->setReversible(false);
			r0->initForwardRateRandomly();
			r0->addReactant(nodes[opIndex]);
			r0->addProduct(modProt);

			Reaction* r1=new Reaction(1);   //add degradation reaction of the modified protein
			r1->setReversible(false);
			r1->initForwardRateRandomly();
			r1->addReactant(modProt);

			if(!existsReaction(*r0) && !existsReaction(*r1)){
				rlist.push_back(r0);
				rlist.push_back(r1);
			}
			else{
				delete r0;
				delete r1;
			}
			return;
		}
		else{                                //if it is a protein complex, the reaction type is AB->A
			Reaction* r2=new Reaction(4);
			r2->setReversible(false);
			r2->initForwardRateRandomly();
			r2->addReactant(nodes[opIndex]);
			int randComp=rand()%(nodes[opIndex]->getNsize()-1)+1;
		    r2->addProduct(nodes[opIndex]->getNode(randComp));        //randomly choose a node in the complex as the product

			Reaction* r3=new Reaction(1);
			r3->setReversible(false);
			r3->initForwardRateRandomly();
			r3->addReactant(nodes[opIndex]->getNode(randComp));
			if(!existsReaction(*r2) && !existsReaction(*r3)){
				rlist.push_back(r2);
				rlist.push_back(r3);
			}
			else{
				delete r2;
				delete r3;
			}
			return;

		}


	}
	else {//	two proteins (protein complex) case
		
		int numOfProt = 0;
		int indexOfProt = 0;
		int numOfSingProt=0;
		int numOfCompProt=0;

		vector<int> protIndice;
		vector<int> singProtIndice;
		vector<int> compProtIndice;
			
		vector<Node*>::iterator iter = nodes.begin();
        vector<Node*>::iterator iter_end = nodes.end();
		while (iter != iter_end) {
			if((*iter)->getNode(0) == NULL) { 	
				if((*iter)->getNtype()==3){
					numOfSingProt++;
					singProtIndice.push_back(indexOfProt);
					numOfProt++;
					protIndice.push_back(indexOfProt);
				}
				if((*iter)->getNtype()==6){
					numOfCompProt++;
					compProtIndice.push_back(indexOfProt);
					numOfProt++;
					protIndice.push_back(indexOfProt);
				}
			}
			indexOfProt++;
			iter ++;
		}

		if(!numOfProt) return;	//	no protein or protein complex
		
		double possibility = (double)rand()/RAND_MAX;

		/*case 1, A+B->AB,dimerization*/

		if(possibility < 1./3) {
			int opIndex1 = protIndice[rand()%numOfProt];	//	protein 1
		    int opIndex2 = protIndice[rand()%numOfProt];	//	protein 2
			Node* dimer = new Node (nodes.size(), nodes[opIndex1], nodes[opIndex2]);
			if(existsNode(*dimer)){
				delete dimer;
				return;
			}
			Reaction* dimerization = new Reaction (5);
			dimerization->setReversible(false);
			dimerization->initForwardRateRandomly();
			dimerization->addReactant(nodes[opIndex1]);
			dimerization->addReactant(nodes[opIndex2]);
			dimerization->addProduct(dimer);

			Reaction* degrad = new Reaction(1);
			degrad->setReversible(false);
			degrad->initForwardRateRandomly();
			degrad->addReactant(dimer);

			if(!existsReaction(*dimerization) && !existsReaction(*degrad)) {
				nodes.push_back(dimer);
				rlist.push_back(dimerization);
				rlist.push_back(degrad);
			}
			else {
				delete dimer;
				delete degrad;
				delete dimerization;
			}
		}

		/*case 2, A+B->A, catalytic degradation*/

		else if(possibility < 2./3) {      
			if(!numOfSingProt) return;  //no single protein.
			int opIndex1=protIndice[rand()%numOfProt];   //protein 1
			int opIndex2=singProtIndice[rand()%numOfSingProt];   //single protein 2
			Reaction* catalyticDeg = new Reaction (6);
			catalyticDeg->setReversible(false);
			catalyticDeg->initForwardRateRandomly();
			catalyticDeg->addReactant(nodes[opIndex1]);
			catalyticDeg->addReactant(nodes[opIndex2]);
			catalyticDeg->addProduct(nodes[opIndex2]);
			if(!existsReaction(*catalyticDeg)){
			    rlist.push_back(catalyticDeg);}
			else
				delete catalyticDeg;
		}

		/*case 3, AB+C->A, partial degradation*/

		else {
			if(!numOfCompProt || !numOfSingProt)  return;  //no protein complex
            int opIndex1=compProtIndice[rand()%numOfCompProt];    //complex protein 1
			int opIndex2=singProtIndice[rand()%numOfSingProt];    //single protein 2
			Reaction* partialDeg=new Reaction(7);
			partialDeg->setReversible(false);
			partialDeg->initForwardRateRandomly();
			partialDeg->addReactant(nodes[opIndex1]);
			partialDeg->addReactant(nodes[opIndex2]);
			int randComp=rand()%(nodes[opIndex1]->getNsize()-1)+1;
		    partialDeg->addProduct(nodes[opIndex1]->getNode(randComp));   //choose the protein in the complex
			if(!existsReaction(*partialDeg)){
			   rlist.push_back(partialDeg);
			}
			else
				delete partialDeg; 
		}

	}
	
	return;
}

#define PROB1 0.5
#define PROB2 1.0
#define PROB3 0.4
#define PROB4 0.1
#define PROB5 0.1

void Cell:: mut_parameters(){
    if (rand() < RAND_MAX*PROB1) {
        mut_deg_prot();
    }
    if (rand() < RAND_MAX*PROB2) {
        mut_kin_const();
	}
}
void Cell:: mut_topology(){
    mut_add_regu();         // we want more regulation between genes and proteins, so topology-mutation always adds new regulation
	if(rand() < PROB3 * RAND_MAX){
		mut_add_gene();
    }
    if (rand() < PROB5 * RAND_MAX) {
        mut_add_postmod();
    }
}
	



//overall mutation method
void Cell::mutation(){
    //srand((unsigned int)time(NULL));
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
/*input requirement: first colomn of data[][] is the initial value
 *
 */
void runge_kutta(double **data,vector<Node*> nodes,vector<Reaction*> rlist ,int numInducers, int series, int steps)
{//100 is steps
	
	int currSerie,currStep;
    double *k1 = new double[series];
    double *k2 = new double[series];
    double *k3 = new double[series];
    double *k4 = new double[series];
    double *y = new double[series], *tempY1 = new double[series], *tempY2 = new double[series];
	
	
	
	//runge_kutta method:
	for(currStep = 0; currStep < steps-1; currStep++){
		//initial values assigned to y[]
		for(currSerie = 0; currSerie < series; currSerie++){
			y[currSerie] = data[currSerie][currStep];
		}
		for(currSerie = 0; currSerie < series; currSerie++){
			double delta = (nodes[currSerie]->ode)(rlist,y,currStep);
            k1[currSerie] = (currSerie < numInducers) ? 0. : delta;
			tempY1[currSerie] =  (y[currSerie] + k1[currSerie]/2.);//(y[currSerie] + k1[currSerie]/2. < 0.) ? 0. : (y[currSerie] + k1[currSerie]/2.);
        }
        for(currSerie = 0; currSerie < series; currSerie++){
            
			double delta = (nodes[currSerie]->ode)(rlist,tempY1,currStep + 0.5);
            k2[currSerie] = (currSerie < numInducers) ? 0. : delta;
			tempY2[currSerie] = (y[currSerie] + k2[currSerie]/2.);//(y[currSerie] + k2[currSerie]/2. < 0.) ? 0. : (y[currSerie] + k2[currSerie]/2.);
        }
        for(currSerie = 0; currSerie < series; currSerie++){
            
			double delta = (nodes[currSerie]->ode)(rlist,tempY2,currStep + 0.5);
            k3[currSerie] = (currSerie < numInducers) ? 0. : delta;
			tempY1[currSerie] = (y[currSerie] + k3[currSerie]);//(y[currSerie] + k3[currSerie] < 0.) ? 0. : (y[currSerie] + k3[currSerie]);
        }
        for(currSerie = 0; currSerie < series; currSerie++){
            
			double delta = (nodes[currSerie]->ode)(rlist,tempY1,currStep + 1.);
            k4[currSerie] = (currSerie < numInducers) ? 0. : delta;
        }
		for (currSerie = numInducers; currSerie < series; currSerie++) {
            data[currSerie][currStep + 1] = y[currSerie] + 1 / 6.0 * (k1[currSerie]+ 2 * k2[currSerie] + 2 * k3[currSerie] + k4[currSerie]);
            if (data[currSerie][currStep + 1] < 0.) {//in case of negative density
                data[currSerie][currStep + 1] = 0;
            }
            if (data[currSerie][currStep + 1] > 100000000.) {// in case of "infinite" density
                data[currSerie][currStep + 1] = 100000000.;
            }
        }	
			
    }
    
    delete [] k1;
    delete [] k2;
    delete [] k3;
    delete [] k4;
    delete [] y;
    delete [] tempY1;
    delete [] tempY2;
}



void Cell::generateTimeCourses(double** targetData,int numTargetNodes, int time){
    
    int size = nodes.size();//  how many nodes in this cell
    /* initialization: store the initial value in the first column of
     * currData, and for coloumn with index greater than number of target
     * nodes, initial value is 0
     */
    this->currData = new double*[size];
	for (int i = 0; i < size; i++) {
        this->currData[i] = new double[time];
        this->currData[i][0] = 1.;   // the initial value of gene is 1
    }
    for (int i = 0; i < numTargetNodes; i++) {
        this->currData[inputIndice[i]][0] = targetData[i][0];    //the initial value of inducers and proteins are the same as the input data.
    }
    for (int i = 0; i < numInducer; i++) {
        for (int j = 0; j < time; j++) {
            this->currData[i][j] = targetData[i][j];
        }
    }
    
    /* runge_kutta method, store the results in currData */
    runge_kutta(this->currData, nodes, rlist, numInducer, size, time);
    
    
}


/*get score using the sfunc as score function and change its own currScore member
 *fist: using Runge-Kutta method to generate the time course and store all them in currData
 *second: using the sfunction as a score function, and passing currData and targetData as parameters to calculate score
 *third: assign the score to currScore
 
 
 *prerequirements: nodes in the cell's "nodes" vector should be sorted by indice
 */
#define PARAMETER_NODE_SIZE 0.2
#define PARAMETER_COMPLEX_SIZE 0.05
#define PARAMETER_REACTION_SIZE 0.1
 
void Cell::getScore(ScoreFunc& sfunc, double** targetData, int numTargetNodes, int time, bool print){
   
    int size = nodes.size();//  how many nodes in this cell
    /* initialization: store the initial value in the first column of
     * currData, and for coloumn with index greater than number of target
     * nodes, initial value is 0
     */
    this->currData = new double*[size];
	for (int i = 0; i < size; i++) {
        this->currData[i] = new double[time];
        this->currData[i][0] = 1.;   // the initial value of gene is 1
    }
    for (int i = 0; i < numTargetNodes; i++) {
        this->currData[inputIndice[i]][0] = targetData[i][0];    //the initial value of inducers and proteins are the same as the input data.
    }
    for (int i = 0; i < numInducer; i++) {
        for (int j = 0; j < time; j++) {
            this->currData[i][j] = targetData[i][j];
        }
    }

    /* runge_kutta method, store the results in currData */
    runge_kutta(this->currData, nodes, rlist, numInducer, size, time);
    
    /* calculate total score uses the ScoreFunc sfunc: only numTargetNodes nodes
     * are calculated because there only those number nodes in targetData
     * therefore, numTargetNodes = numind + numprot
     */
    double totalScore = 0;
    for (int i = 0; i < numTargetNodes; i++) {
        totalScore += sfunc.getScore(this->currData[inputIndice[i]],targetData[i],time);  //compare the RK data and the input data, using score function.
    }
    
    currScore = totalScore;
	int complex_size = 0;
	int node_size = nodes.size();
	for(int i = 0; i < node_size; i++){
		if(nodes[i] -> getNtype() == 5)
			complex_size += ((nodes[i] -> getNsize()) - 2);
		else
			if(nodes[i] -> getNtype() == 6)
				complex_size += ((nodes[i] -> getNsize()) - 3);
	}
	currScore += PARAMETER_NODE_SIZE * nodes.size();
	currScore += PARAMETER_REACTION_SIZE * rlist.size();
	currScore += PARAMETER_COMPLEX_SIZE * complex_size;
	

    
    if (print) {
        for (int i = 0; i < numTargetNodes; i++) {
            for (int j = 0; j < time; j++) {
                std::cout << this->currData[inputIndice[i]][j] << "\t";
            }
            std::cout << std::endl;
        }
    }
    
    for (int i = 0; i < size; i++) {
        delete [] this->currData[i];
    }
    delete [] this->currData;
}

	
//return the currSore generated by getScore method
double Cell::getCurrScore(){
    //std::cout << currScore << std::endl;
    return currScore;
}



using namespace std;
//description method
void Cell::description(int time){
    
    //print nodes
    cout << " Nodes:" << endl;
    vector<Node*>::iterator iter = nodes.begin();
    while (iter != nodes.end()) {
        cout << "  Node " << (*iter)->getNindex() << ": " << (*iter)->getNstring() << endl;
        iter++;
    }
    
    cout << endl << endl << " Reactions:" << endl;
    vector<Reaction*>::iterator iter_reaction = rlist.begin();
    int index = 1;
    while (iter_reaction != rlist.end()) {
        (*iter_reaction)->description(index);
        index++;
        iter_reaction++;
    }
    
    
    cout << "Regulatory Matrix: " << endl;
    genRegulatoryRelationships();
    
    for (int i = 0; i < nodes.size(); i++) {
        delete [] currData[i];
    }
    delete [] currData;
}


std::vector<Node*>* Cell::getNodesVector(){
    return &nodes;
}

std::vector<Reaction*>* Cell::getRlistVector(){
    return &rlist;
}
void Cell::genRegulatoryRelationships(){
	int numOfProtLike = 0; 
	int numOfGene = 0;
	int sizeOfNode = nodes.size();
	// figure out the quantity of genes and prots
	for(int i = 0; i < sizeOfNode; i++){
		if(nodes[i] -> getNtype() == 3 || nodes[i] -> getNtype() == 6)
			numOfProtLike ++;
		else
			if(nodes[i] -> getNtype() == 2)
				numOfGene ++;
	}
	// initial function of regulatoryMatrix to zero
	regulatoryMatrix = new int* [numOfGene];
	for(int j = 0; j < numOfGene; j++){
		regulatoryMatrix[j] = new int [numOfProtLike];
	}
	for(int i = 0; i < numOfGene; i++)
		for(int j = 0; j < numOfProtLike; j++){
			regulatoryMatrix[i][j] = 0;
		}
	// to restore the index of genes and prots
	int* indexOfProtLike = new int [numOfProtLike];
	int* indexOfGene = new int[numOfGene];
	int h = 0;
	int k = 0;
	for(int i = 0; i < sizeOfNode; i++){
		if(nodes[i] -> getNtype() == 3 || nodes[i] -> getNtype() == 6){
			indexOfProtLike[h] = nodes[i] -> getNindex();
			h++;
		}
		else
			if(nodes[i] -> getNtype() == 2){
				indexOfGene[k] = nodes[i] -> getNindex();
				k++;
			}
	}
	int sizeOfReaction = rlist.size();
	int indexOfNewMod;
	int indexOfOldMod;
	int indexOfTargetProt;
	// find the reaction of type3
	for(int l = 0; l < sizeOfReaction; l++){
		if(rlist[l] -> getRtype() == 2){
			indexOfNewMod = rlist[l] -> getProduct(0) -> getNindex();
			if(rlist[l] -> getReactant(0) -> getNtype() == 1){
				indexOfOldMod = rlist[l] -> getReactant(0) -> getNindex();
		        indexOfTargetProt = rlist[l] -> getReactant(1) -> getNindex();
			}
			else{
				indexOfOldMod = rlist[l] -> getReactant(1) -> getNindex();
		        indexOfTargetProt = rlist[l] -> getReactant(0) -> getNindex();
			}
		double forwardrateOfNew;
		double forwardrateOfOld;
		// get the forwardrate of the old and the new
		for(int p = 0; p < sizeOfReaction; p++){
			if(rlist[p] -> getRtype() == 0){
				if(rlist[p] -> getModifier(0) -> getNindex() == indexOfNewMod)
					forwardrateOfNew = rlist[p] -> getForwardRate();
				else
					if(rlist[p] -> getModifier(0) -> getNindex() == indexOfOldMod)
						forwardrateOfOld = rlist[p] -> getForwardRate();
			}
		}
		int s = 0; // x 
		int t = 0; // y
		// find the possition of selected prot and gene
		for(int q = 0; q < numOfGene; q++){
			if(indexOfGene[q] == indexOfOldMod){
                s = q;
                break;
            }
		}
		for(int r = 0; r < numOfProtLike; r++){
			if(indexOfProtLike[r] == indexOfTargetProt){
				t = r;
                break;
            }
		}
		if(forwardrateOfNew > forwardrateOfOld)
			regulatoryMatrix[s][t] = 1;
		else
			regulatoryMatrix[s][t] = -1;
		}
	}
    for (int i = 0; i < numOfGene ; i++) {
        for (int j = 0; j < numOfProtLike; j++) {
            std::cout << regulatoryMatrix[i][j] << "\t";
        }
        std::cout << endl;
    }
	delete [] indexOfProtLike;
	delete [] indexOfGene;
}



				
		
		
		









