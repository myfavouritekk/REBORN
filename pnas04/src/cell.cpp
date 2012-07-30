#include "cell.h"

namespace ustc{

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
    rankings = cell.rankings;
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

//add ranking to the rankings vector
void Cell::addRanking(int ranking){
    rankings.push_back(ranking);
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
			vector<Node*>::iterator iter = nodes.begin();
			int indexOfModifier=0;
			int numOfModifier=0;
			vector<int> modifierIndex;
			while(iter != nodes.end()){     //choose nodes with type 1,3,4,6,then choose one as the modifier
				int _Ntype=(*iter)->getNtype();
				if(_Ntype ==1 || _Ntype ==3 || _Ntype==4 || _Ntype == 6){
					numOfModifier++;
					modifierIndex.push_back(indexOfModifier);
				}
				iter++;
				indexOfModifier++;
			}
			if(!numOfModifier) return;

			int opModifierIndex = modifierIndex[rand()%numOfModifier];  //choose one modifier

			Node* modProt = new Node(nodes.size(),3);
			nodes.push_back(modProt);

			Reaction* r0=new Reaction(3);    //add modification reaction
			r0->setReversible(false);
			r0->initForwardRateRandomly();
			r0->addReactant(nodes[opIndex]);
			r0->addModifier(nodes[opModifierIndex]);
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



void Cell:: mut_parameters(){
    if (rand() < RAND_MAX*PROB_MUT_DEG_PROT) {
        mut_deg_prot();
    }
    if (rand() < RAND_MAX*PROB_MUT_KIN_CONST) {
        mut_kin_const();
	}
}
void Cell:: mut_topology(){
    mut_add_regu();         // we want more regulation between genes and proteins, so topology-mutation always adds new regulation
	if(rand() < PROB_MUT_ADD_GENE * RAND_MAX){
		mut_add_gene();
    }
    if (rand() < PROB_MUT_ADD_POSTMOD * RAND_MAX) {
        mut_add_postmod();
    }
}
	

void Cell::mut_parameters_simAnneal(){
    
	if(!rlist.size()) return;
	
	//	all kinetic constantsto be modified
	
    std::vector<Reaction*>::iterator iter_reaction = rlist.begin();
    std::vector<Reaction*>::iterator iter_reaction_end = rlist.end();
    while (iter_reaction != iter_reaction_end) {
        if((*iter_reaction) -> getRtype() != 8){
            (*iter_reaction)->modifyForwardRate();
            (*iter_reaction)->modifyReverseRate();
        }
        iter_reaction++;
    }
    
	return;

}


//overall mutation method
void Cell::mutation(){
    //srand((unsigned int)time(NULL));
    if (rand() < RAND_MAX*PROB_MUT_DEG_PROT) {
        mut_deg_prot();
    }
    if (rand() < RAND_MAX*PROB_MUT_KIN_CONST) {
        mut_kin_const();
    }
    if (rand() < RAND_MAX*PROB_MUT_ADD_GENE) {
        mut_add_gene();
    }
    if (rand() < RAND_MAX*PROB4) {
        mut_add_regu();
    }
    if (rand() < RAND_MAX*PROB_MUT_ADD_POSTMOD) {
        mut_add_postmod();
    }
}

//runge-kutta method
/*input requirement: first colomn of data[][] is the initial value
 *
 */
void runge_kutta(double **data,vector<Node*> nodes,vector<Reaction*> rlist ,int numInducers, int series, int steps,int *timeOfAddInducers)
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
            k1[currSerie] = (currSerie < numInducers && timeOfAddInducers[currSerie] < currStep) ? 0. : delta; 
			tempY1[currSerie] =  (y[currSerie] + k1[currSerie]/2.);//(y[currSerie] + k1[currSerie]/2. < 0.) ? 0. : (y[currSerie] + k1[currSerie]/2.);
        }
        for(currSerie = 0; currSerie < series; currSerie++){
            
			double delta = (nodes[currSerie]->ode)(rlist,tempY1,currStep + 0.5);
            k2[currSerie] = (currSerie < numInducers && timeOfAddInducers[currSerie] < currStep) ? 0. : delta;
			tempY2[currSerie] = (y[currSerie] + k2[currSerie]/2.);//(y[currSerie] + k2[currSerie]/2. < 0.) ? 0. : (y[currSerie] + k2[currSerie]/2.);
        }
        for(currSerie = 0; currSerie < series; currSerie++){
            
			double delta = (nodes[currSerie]->ode)(rlist,tempY2,currStep + 0.5);
            k3[currSerie] = (currSerie < numInducers && timeOfAddInducers[currSerie] < currStep) ? 0. : delta;
			tempY1[currSerie] = (y[currSerie] + k3[currSerie]);//(y[currSerie] + k3[currSerie] < 0.) ? 0. : (y[currSerie] + k3[currSerie]);
        }
        for(currSerie = 0; currSerie < series; currSerie++){
            
			double delta = (nodes[currSerie]->ode)(rlist,tempY1,currStep + 1.);
            k4[currSerie] = (currSerie < numInducers && timeOfAddInducers[currSerie] < currStep) ? 0. : delta;
        }
		for (currSerie = 0; currSerie < series; currSerie++) {
            if(currSerie<numInducers && currStep <= timeOfAddInducers[currSerie]){//if currSerie>=numInducers, the right term will not be executed, notice"="!!. 
				continue;
			}

            data[currSerie][currStep + 1] = y[currSerie] + 1 / 6.0 * (k1[currSerie]+ 2 * k2[currSerie] + 2 * k3[currSerie] + k4[currSerie]);
            if (data[currSerie][currStep + 1] < 0.) {//in case of negative density
                data[currSerie][currStep + 1] = 0;
            }
            if (data[currSerie][currStep + 1] > 1000.) {// in case of "infinite" density
                data[currSerie][currStep + 1] = 1000.;
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
    int* timeOfAddInducers = new int[numInducer];
	for(int i = 0; i < numInducer; i++){
		for(int j = 0; j < time; j++){
			if(targetData[i][j]>0.0000001){
				timeOfAddInducers[i] = j;
				break;
			}
		}
	}

    int size = nodes.size();//  how many nodes in this cell
    /* initialization: store the initial value in the first column of
     * currData, and for coloumn with index greater than number of target
     * nodes, initial value is 0
     */
    this->currData = new double*[size];
	for (int i = 0; i < size; i++) {
        this->currData[i] = new double[time];
        this->currData[i][0] = 1.;   // the initial value of gene is 1
		if((*(this->getNodesVector()))[i]->getNtype() == 4)
			this->currData[i][0]=0;
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
    runge_kutta(this->currData, nodes, rlist, numInducer, size, time,timeOfAddInducers);
    delete[] timeOfAddInducers;
    
}

/*function: pring time course to a file with the name of a string parameter
 *prerequisite: currData must be generated, so use this method after generateTimeCourses method, and before description method
 *parameter: "name" contains the evolution generation and rankin of this cell, "time" is the total points of each time course
 */
 void Cell::printCurrDataToAFile(std::string name, int time){
    
    //print time courses to file
    std::ofstream timeCoursesFile;
    timeCoursesFile.open(name.c_str());
    int nodeSize = nodes.size();
    for (int i = 0; i < nodeSize; i++) {
        for (int j = 0; j < time; j++) {
            timeCoursesFile << this->currData[i][j] << "\t";
        }
        timeCoursesFile << std::endl;
    }
    timeCoursesFile.close();
}



/*get score using the sfunc as score function and change its own currScore member
 *fist: using Runge-Kutta method to generate the time course and store all them in currData
 *second: using the sfunction as a score function, and passing currData and targetData as parameters to calculate score
 *third: assign the score to currScore
 *prerequirements: nodes in the cell's "nodes" vector should be sorted by indice
 */
void Cell::getScore(ScoreFunc& sfunc, double** targetData, int numTargetNodes, int time, bool print){
	int* timeOfAddInducers = new int[numInducer];
	for(int i = 0; i < numInducer; i++){
		for(int j = 0; j < time; j++){
			if(targetData[i][j]>0.0000001){
				timeOfAddInducers[i] = j;
				break;
			}
		}
	}

    int size = nodes.size();//  how many nodes in this cell
    /* initialization: store the initial value in the first column of
     * currData, and for coloumn with index greater than number of target
     * nodes, initial value is 0
     */
    this->currData = new double*[size];
	for (int i = 0; i < size; i++) {
        this->currData[i] = new double[time];
        this->currData[i][0] = 1.;   // the initial value of gene is 1
		if((*(this->getNodesVector()))[i]->getNtype() == 4)
			this->currData[i][0]=0;
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
    runge_kutta(this->currData, nodes, rlist, numInducer, size, time,timeOfAddInducers);
    delete[] timeOfAddInducers;
    /* calculate total score uses the ScoreFunc sfunc: only numTargetNodes nodes
     * are calculated because there only those number nodes in targetData
     * therefore, numTargetNodes = numind + numprot
     */
    double totalScore = 0;
    for (int i = 0; i < numTargetNodes; i++) {
        totalScore += sfunc.getScore(targetData[i],this->currData[inputIndice[i]],time);  //compare the RK data and the input data, using score function.
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
	

    //print the time courses
    if (print) {
        
        //print time courses to file
        std::ofstream timeCoursesFile;
        timeCoursesFile.open("data.txt");
        for (int i = 0; i < numTargetNodes; i++) {
            for (int j = 0; j < time; j++) {
                std::cout << this->currData[inputIndice[i]][j] << "\t";
                timeCoursesFile << this->currData[inputIndice[i]][j] << "\t";
            }
            std::cout << std::endl;
            timeCoursesFile << std::endl;
        }
        timeCoursesFile.close();
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



//description method
void Cell::description(int time){
    
    //print score
    std::cout << " Score: " << currScore << std::endl;
    
    //print rankings
    std::cout << " Rankings: ";
    vector<int>::iterator iter_ranking = rankings.begin();
    while (iter_ranking != rankings.end()) {
        std::cout << *iter_ranking << "\t";
        iter_ranking++;
    }
    std::cout << std::endl;
    
    //print nodes
    std::cout << " Nodes:" << std::endl;
    vector<Node*>::iterator iter = nodes.begin();
    while (iter != nodes.end()) {
        std::cout << "  Node " << (*iter)->getNindex() << ": " << (*iter)->getNstring() << std::endl;
        iter++;
    }
    
    std::cout << std::endl << std::endl << " Reactions:" << std::endl;
    vector<Reaction*>::iterator iter_reaction = rlist.begin();
    int index = 1;
    while (iter_reaction != rlist.end()) {
        (*iter_reaction)->description(index);
        index++;
        iter_reaction++;
    }
    
    
    std::cout << "Regulatory Matrix: " << std::endl;
    genRegulatoryRelationships();
    
    std::vector<Motif*>::iterator iter_motif = motifs.begin();
    std::vector<Motif*>::iterator iter_motif_end = motifs.end();
    while (iter_motif != iter_motif_end) {
        (*iter_motif) -> description();
        iter_motif++;
    }
    
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

std::vector<int>* Cell::getInputIndiceVector(){
	return &inputIndice;
}

    
    
/*genRegulatoryRelationships() mehtod:
 *purpose: this method is to judge whether the genes in the cell have regulatory relationships,
 *  and put the regulatory relationships in the regulatoryMatrix array
 */
void Cell::genRegulatoryRelationships(){
    //int numOfProtLike = 0;
    int numOfGene = 0;
    int numberOfNodes = nodes.size();
    // figure out the quantity of genes
    for(int i = 0; i < numberOfNodes; i++){
        if(nodes[i] -> getNtype() == 2)
            numOfGene ++;
    }
    // initial function of regulatoryMatrix to 0
    regulatoryMatrix = new int* [numOfGene];
    for(int j = 0; j < numOfGene; j++){
        regulatoryMatrix[j] = new int [numOfGene];
    }
    for(int i = 0; i < numOfGene; i++)
        for(int j = 0; j < numOfGene; j++){
            regulatoryMatrix[i][j] = 0;
        }
    
    // to store the indice of genes and prots
    int* indexOfGene = new int[numOfGene];
    int k = 0;
    for(int i = 0; i < numberOfNodes; i++){
        if(nodes[i] -> getNtype() == 2){
            indexOfGene[k] = nodes[i] -> getNindex();
            k++;
        }
        
    }
    int numberOfReactions = rlist.size();
    int indexOfNewMod;
    int indexOfOldMod;
    int indexOfRegualtingProtLike;
    
    // find the reaction of Type 2: protein binding gene
    for(int l = 0; l < numberOfReactions; l++){
        if(rlist[l] -> getRtype() == 2){
            indexOfNewMod = rlist[l] -> getProduct(0) -> getNindex();//the new modifier is the a gene/protein complex, which is the first and only product fo a binding reaction
            if(rlist[l] -> getReactant(0) -> getNtype() == 2){  //node type 2 is gene
                indexOfOldMod = rlist[l] -> getReactant(0) -> getNindex();
                indexOfRegualtingProtLike = rlist[l] -> getReactant(1) -> getNindex();
            }
            else{
                indexOfOldMod = rlist[l] -> getReactant(1) -> getNindex();
                indexOfRegualtingProtLike = rlist[l] -> getReactant(0) -> getNindex();
            }
            
            double forwardrateOfNew;
            double forwardrateOfOld;
            // get the forwardrate of the old and the new transcription reactions
            for(int i = 0; i < numberOfReactions; i++){
                if(rlist[i] -> getRtype() == 0){//reaction type 0 is transcription
                    if(rlist[i] -> getModifier(0) -> getNindex() == indexOfNewMod)
                        forwardrateOfNew = rlist[i] -> getForwardRate();
                    else
                        if(rlist[i] -> getModifier(0) -> getNindex() == indexOfOldMod)
                            forwardrateOfOld = rlist[i] -> getForwardRate();
                }
            }
            int row = 0; // column of a gene
            int column = 0; // row of a gene
            // find the possition of selected prot and gene
            for(int i = 0; i < numOfGene; i++){
                if(indexOfGene[i] == indexOfOldMod){
                    row = i;
                    break;// get the row which the old modifier is in and stop looping
                }
            }
            
            //regualtion protein like nodes can be protein(type 3) or protein complex(type 6)
            int typeOfRegulatingProtLike = nodes[indexOfRegualtingProtLike] -> getNtype();
            switch (typeOfRegulatingProtLike) {
                case 3:{//single protain case
                    int indexOfRegulatingGene = 0;
                    for(int i = 0; i < numberOfReactions; i++){//find the gene that transcripts the regulating protein
                        if(rlist[i] -> getRtype() == 0){//transcription
                            if(rlist[i] -> getProduct(0) -> getNindex() == indexOfRegualtingProtLike){
                                indexOfRegulatingGene = rlist[i] -> getModifier(0) -> getNindex();
                                break;
                            }
                        }
                    }
                    for(int i = 0; i < numOfGene; i++){
                        if(indexOfGene[i] == indexOfRegulatingGene){
                            column = i;
                            if(forwardrateOfNew > forwardrateOfOld)
                                regulatoryMatrix[row][column] = 1;
                            else
                                regulatoryMatrix[row][column] = -1;
                            
                            break;
                        }
                    }
                    break;
                }
                    
                case 6:{//protein complex case
                    int numberOfRegulatingGenes = nodes[indexOfRegualtingProtLike] -> getNsize() - 1;//WRONG: the size of complex is not the number of regulating genes plus 1, because of cases like P1:P1:P3, but this will still work because it will assign to g0 twice
                    int* indiceOfRegualtingGenes = new int[numberOfRegulatingGenes];
                    for(int i = 1; i <= numberOfRegulatingGenes; i ++){//starts from 1 because components[0] is for gene or NULL, and "=" in "<=" is important!
                        indiceOfRegualtingGenes[i] = nodes[indexOfRegualtingProtLike] -> getNode(i) -> getNindex();
                        for(int i = 0; i < numberOfReactions; i++){
                            if(rlist[i] -> getRtype() == 0){
                                if(rlist[i] -> getProduct(0) -> getNindex() == indiceOfRegualtingGenes[i]){
                                    indiceOfRegualtingGenes[i] = rlist[i] -> getModifier(0) -> getNindex();
                                    for(int i = 0; i < numOfGene; i++){
                                        if(indexOfGene[i] == indiceOfRegualtingGenes[i]){
                                            column = i;
                                            if(forwardrateOfNew > forwardrateOfOld)
                                                regulatoryMatrix[row][column] = 1;
                                            else
                                                regulatoryMatrix[row][column] = -1;
                                            break;
                                        }
                                    }
                                }
                            }
                        }
                    }
                    delete[] indiceOfRegualtingGenes;
                    break;
                }

                default:
                    break;
            }
        }
    }
    
    
    
    
    //output the regulatoryMatix
    /*the structure fo output looks like
                g0  g2  g4
             g0 0   1   0
             g2 -1  0   0
             g4 -1  1   0
     */
    std::cout<< "\t";
    for (int i = 0; i < numOfGene ; i++){
        std::cout<< nodes[indexOfGene[i]] -> getNstring()<< "\t";
    }
    std::cout<< std::endl;
    
    for (int i = 0; i < numOfGene ; i++) {
        std::cout<< nodes[indexOfGene[i]] -> getNstring()<< "\t";
        for (int j = 0; j < numOfGene; j++) {
            std::cout << regulatoryMatrix[i][j] << "\t";
        }
        std::cout << std::endl;
    }
    delete [] indexOfGene;

    
    
    //find and print motifs
    findMotifs();
    
}

    

//find motifs consisting 1, 2 or 3 genes
void Cell::findMotifs(){
    int numOfGenes = 0;
    std::vector<Node*>::iterator iter_node = nodes.begin();
    std::vector<Node*>::iterator iter_node_end = nodes.end();
    while (iter_node != iter_node_end) {
        int nodeType = (*iter_node)->getNtype();
        if (nodeType == 2) {// Type 2 is gene
            numOfGenes++;
        }
        iter_node++;
    }
    
    //store the indice of genes
    int* indiceOfGenes = new int[numOfGenes];
    for (int i = 0 ; i < numOfGenes; i++) {
        int nodeType = nodes[i] -> getNtype();
        if (nodeType == 2) {//Type 2 is gene
            indiceOfGenes[i] = nodes[i] -> getNindex();
        }
    }
    
    //find and print motifs consisting of 1, 2 and 3 genes
    findSingleMotifs(numOfGenes, indiceOfGenes);
    findDoubleMotifs(numOfGenes, indiceOfGenes);
    findTripleMotifs(numOfGenes, indiceOfGenes);
    
    delete [] indiceOfGenes;

}

    
//find, store and print single-gene motifs
void Cell::findSingleMotifs(int numberOfGenes, int* indiceOfGenes){

    for(int i = 0; i < numberOfGenes ; i++){
		bool isSingle=1;

		if(regulatoryMatrix[i][i] == 0){
			continue;
		}

		for(int j = 0; j < numberOfGenes; j++){
			if(j==i) continue;
			if(regulatoryMatrix[i][j]!=0 || regulatoryMatrix[j][i]!=0 ){
				isSingle=0;
				break;
			}
		}

		if(isSingle){
			int** motifMatrix = new int*[1];
            motifMatrix[0] = new int[1];
            
            std::cout<<"Single Mofif:"<<std::endl;
			std::cout<<"\t";
			std::cout<<nodes[indiceOfGenes[i]] -> getNstring()<<"\t";
			std::cout<<std::endl;
			std::cout<<nodes[indiceOfGenes[i]] -> getNstring()<<"\t";
			std::cout<<regulatoryMatrix[i][i];
            motifMatrix[0][0] = regulatoryMatrix[i][i];
			std::cout<<std::endl;
            
            std::vector<Node*> motifNodes;
            motifNodes.push_back(nodes[indiceOfGenes[i]]);
            Motif* singleMotif = new Motif(&motifNodes, motifMatrix);
            motifs.push_back(singleMotif);
            
            delete [] motifMatrix[0];
            delete [] motifMatrix;
		}
	}
}

    
//find, store and print two-gene motifs
void Cell::findDoubleMotifs(int numberOfGenes, int* indiceOfGenes){
	for(int i = 0; i < numberOfGenes; i++){
		for(int j=i+1;j < numberOfGenes; j++){
			bool isDouble=1;

			if(regulatoryMatrix[i][j]==0 && regulatoryMatrix[j][i]==0){
				continue;
			}

			for(int k=0;k < numberOfGenes; k++){
				if(k==i||k==j) continue;
				if(regulatoryMatrix[i][k]!=0||regulatoryMatrix[k][i]!=0||regulatoryMatrix[j][k]!=0||regulatoryMatrix[k][j]!=0){
					isDouble=0;
					break;
				}
			}

			if(isDouble){
                int** motifMatrix = new int*[2];
                motifMatrix[0] = new int[2];
                motifMatrix[1] = new int[2];
                for (int i = 0; i < 2; i++) {
                    for (int j = 0;  j < 2; j++) {
                        motifMatrix[i][j] = 0;
                    }
                }
                
				std::cout<<"Double Motif:"<<std::endl;
				std::cout<<"\t";
				std::cout<<nodes[indiceOfGenes[i]] -> getNstring()<<"\t";
				std::cout<<nodes[indiceOfGenes[j]] -> getNstring()<<"\t"<<std::endl;
				std::cout<<nodes[indiceOfGenes[i]] -> getNstring()<<"\t"<<regulatoryMatrix[i][i]<<"\t"<<regulatoryMatrix[i][j];
				std::cout<<std::endl;
				std::cout<<nodes[indiceOfGenes[j]] -> getNstring()<<"\t"<<regulatoryMatrix[j][i]<<"\t"<<regulatoryMatrix[j][j];
				std::cout<<std::endl;
                
                //store motif matrix
                motifMatrix[0][0] = regulatoryMatrix[i][i];
                motifMatrix[0][1] = regulatoryMatrix[i][j];
                motifMatrix[1][0] = regulatoryMatrix[j][i];
                motifMatrix[1][1] = regulatoryMatrix[j][j];
                
                //constructing a double motif
                std::vector<Node*> motifNodes;
                motifNodes.push_back(nodes[indiceOfGenes[i]]);
                motifNodes.push_back(nodes[indiceOfGenes[j]]);
                Motif* doubleMotif = new Motif(&motifNodes, motifMatrix);
                motifs.push_back(doubleMotif);
                
                delete [] motifMatrix[0];
                delete [] motifMatrix[1];
                delete [] motifMatrix;

			}
		}
	}  
}

    
//find, store and print three-gene motifs
void Cell::findTripleMotifs(int numberOfGenes, int* indiceOfGenes){
    	
    for(int i = 0;i < numberOfGenes; i++){
		for(int j = i + 1; j < numberOfGenes; j ++){
			if(regulatoryMatrix[i][j] == 0 && regulatoryMatrix[j][i] == 0){
				continue;
			}
			for(int k = j + 1; k < numberOfGenes; k ++){
				if(regulatoryMatrix[i][k] == 0 && regulatoryMatrix[k][i] == 0 && regulatoryMatrix[j][k] == 0 && regulatoryMatrix[k][j] == 0){
					continue;
				}
				else{
                    
                    int** motifMatrix = new int*[3];
                    for (int i = 0; i < 3; i++) {
                        motifMatrix[i] = new int[3];
                        for (int j = 0; j < 3; j++) {
                            motifMatrix[i][j] = 0;
                        }
                    }
                    
					std::cout << "\t" << nodes[indiceOfGenes[i]] -> getNstring() << "\t" << nodes[indiceOfGenes[j]] -> getNstring() << "\t" << nodes[indiceOfGenes[k]] -> getNstring() << std::endl;
					std::cout << nodes[indiceOfGenes[i]] -> getNstring() << "\t" <<regulatoryMatrix[i][i] << "\t" << regulatoryMatrix[i][j] << "\t" << regulatoryMatrix[i][k] << std::endl;
					std::cout << nodes[indiceOfGenes[j]] -> getNstring() << "\t" <<regulatoryMatrix[j][i] << "\t" << regulatoryMatrix[j][j] << "\t" << regulatoryMatrix[j][k] << std::endl;
					std::cout << nodes[indiceOfGenes[k]] -> getNstring() << "\t" <<regulatoryMatrix[k][i] << "\t" << regulatoryMatrix[k][j] << "\t" << regulatoryMatrix[k][k] << std::endl;
                    
                    //store motif matrix
                    int index[3];
                    index[0] = i;
                    index[1] = j;
                    index[2] = k;
                    for (int m = 0; m < 3; m++) {
                        for (int n = 0; n < 3; n++) {
                            motifMatrix[m][n] = regulatoryMatrix[index[m]][index[n]];
                        }
                    }
                    
                    //constructing a triple motif
                    std::vector<Node*> motifNodes;
                    motifNodes.push_back(nodes[indiceOfGenes[i]]);
                    motifNodes.push_back(nodes[indiceOfGenes[j]]);
                    motifNodes.push_back(nodes[indiceOfGenes[k]]);
                    Motif* tripleMotif = new Motif(&motifNodes, motifMatrix);
                    motifs.push_back(tripleMotif);
                    
                    delete [] motifMatrix[0];
                    delete [] motifMatrix[1];
                    delete [] motifMatrix[2];
                    delete [] motifMatrix;

				}
			}
		}
	}
    

}

    
    
    

//void  addReaction(int _rtype,int index)(A -> A*)
void Cell::addReaction(int _rtype,int index){
	if(_rtype !=3){
		std::cout<<"Wrong Reaction Type!"<<std::endl;
		return;
	}
	if(!existsNode(*nodes[index])){
		std::cout<<"No Such Node!"<<std::endl;
		return;
	}

	vector<Node*>::iterator iter = nodes.begin();
	int indexOfModifier=0;
	int numOfModifier=0;
	vector<int> modifierIndex;
	while(iter != nodes.end()){     //choose nodes with type 1,3,4,6,then choose one as the modifier
		int _Ntype=(*iter)->getNtype();
		if(_Ntype ==1 || _Ntype ==3 || _Ntype==4 || _Ntype == 6){
			numOfModifier++;
			modifierIndex.push_back(indexOfModifier);
		}
		iter++;
		indexOfModifier++;
	}
	if(!numOfModifier) return;

	int opModifierIndex = modifierIndex[rand()%numOfModifier];  //choose one modifier

	Node* modifiedProt = new Node(nodes.size(),3);
	Reaction* modification = new Reaction(_rtype);

	modification -> setReversible(false);
	modification -> initForwardRateRandomly();
	modification -> addReactant(nodes[index]);
	modification -> addModifier(nodes[opModifierIndex]);
	modification -> addProduct(modifiedProt);

	Reaction* degradation = new Reaction(1);
	degradation -> setReversible(false);
	degradation -> initForwardRateRandomly();
	degradation -> addReactant(modifiedProt);

	nodes.push_back(modifiedProt);
	rlist.push_back(modification);
	rlist.push_back(degradation);
}

//void addReaction(int _rtype,int firstIndex,int secondIndex)
void Cell::addReaction(int _rtype,int firstIndex,int secondIndex){

	if(!existsNode(*nodes[firstIndex]) || !existsNode(*nodes[secondIndex])){
		std::cout<<"No Such Nodes!"<<std::endl;
		return;
	}
	
	switch(_rtype){

	case 8:             //indu + prot -> indu:prot
		{
		Reaction* r1=new Reaction(8);
		Node* induProt=new Node(nodes.size(),4,nodes[firstIndex],nodes[secondIndex]);
		nodes.push_back(induProt);

		r1->setReversible(true);
		r1->initForwardRateRandomly();
		r1->initReverseRateRandomly();
		r1->addReactant(nodes[firstIndex]);
		r1->addReactant(nodes[secondIndex]);
		r1->addProduct(induProt);

		rlist.push_back(r1);
		break;
		}
	case 5:             //prot1 + prot2 -> prot1:prot2
		{
		Reaction *r1 = new Reaction(5);
		Reaction *r2 = new Reaction(1); 
		Node* dimer = new Node(nodes.size(),6,nodes[firstIndex],nodes[secondIndex]);
		nodes.push_back(dimer);

		r1->setReversible(false);
		r1->initForwardRateRandomly();
		r1->addReactant(nodes[firstIndex]);
		r1->addReactant(nodes[secondIndex]);
		r1->addProduct(dimer);

		r2->setReversible(false);
		r2->initForwardRateRandomly();
		r2->addReactant(dimer);

		rlist.push_back(r1);
		rlist.push_back(r2);
		break;
		}
	case 2:             //prot + gene -> gene::prot
		{
		Reaction *r1 = new Reaction(2);
		Reaction *r2 = new Reaction(0);
		Node * binding = new Node(nodes.size(),5,nodes[firstIndex],nodes[secondIndex]);
		nodes.push_back(binding);

		r1->setReversible(true);
		r1->initForwardRateRandomly();
		r1->initReverseRateRandomly();
		r1->addReactant(nodes[firstIndex]);
		r1->addReactant(nodes[secondIndex]);
		r1->addProduct(binding);

		Node* exGene = binding->extractFirstGene();
		if(exGene == NULL) return;

		Node* exProt = NULL;
		std::vector<Reaction*>::iterator iter = rlist.begin();
		while (iter != rlist.end()) {
        if((*iter)->getRtype() == 0) { //   #0 is gene transcription
			Node* sr = (*iter)->getModifier(0);
			if(sr != NULL && (sr->getNindex() == exGene->getNindex())) {
				exProt = (*iter)->getProduct(0);
				break;
				}
		 }
			iter ++;
		}

		if(exProt == NULL) return;

		r2->setReversible(false);
		r2->initForwardRateRandomly();
		r2->addModifier(binding);
		r2->addProduct(exProt);

		rlist.push_back(r1);
		rlist.push_back(r2);
		break;
		}
	default:
		{
		std::cout<<"Wrong Reaction Type!"<<std::endl;
		break;
		}
	}

}

 
}   //namespace ustc
				
		
		
		









