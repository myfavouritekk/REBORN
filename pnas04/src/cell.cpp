#include "cell.h"

void Cell::mut_deg_prot () {

	if(!rlist.size()) return;

    //  the degradation rate of a protein is modified

    int num1 = 0;
    int num2 = 0;
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
    int opIndex = indice[rand()%num1];

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

	Node* gene = new Node(currNI,0);	//	gene
	nodes.push_back(gene);

	Node* prot = new Node(currNI+1,1);	//	prot
	nodes.push_back(prot);
	
	//	create reaction 0, transcription
	Reaction* r0 = new Reaction ();
	r0->setReversible(false);
	r0->initForwardRateRandomly();
	r0->addModifier(gene);
	r0->addProduct(prot);
	rlist->push_back(r0);

	//	create reaction 1, protein degradation
	Reaction* r1 = new Reaction ();
	r1->setReversible(false);
	r1->initForwardRateRandomly();
	r1->addReactant(prot);
	rlist->push_back(r1);

	return;
}

void Cell::mut_add_regu () {
	
	//A new interaction between a protein and a gene or a gene/protein complex is introduced
    int num1 = 0;
    int num2 = 0;
	int num3 = 0;
    vector<int> protIndice;
	vector<int> cplxIndice;

    //  scan <rlist> and choose  gene and one gene/promoter complex
    std::vector<Node*>::iterator iter1 = nodes.begin();
    while (iter1 != nodes.end()) {
        if((*iter1)->getNtype() == 1) { //   #1 is a protein
            num1++;
            protIndice.push_back(num3);
        }
		
		if((*iter1)->getNtype() == 2 && (*iter1)->isGeneProteinComplex()) {
			num2++;
			cplxIndice.push_back(num3);
		}

        num3++;
        iter1++;
    }

	if(!num1 || !num2) return; // no gene or gene/protein complex

	//	sample random number
    srand(time(NULL));
    int opIndex1 = protIndice[rand()%num1];
    int opIndex2 = cplxIndice[rand()%num2];
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
	Node* ncomplex = new Node(nodes.size(), 2, nodes[opIndex1], nodes[opIndex2]);
	nodes.push_back(ncomplex);

	//	create reaction 0, transcription
	Reaction* r0 = new Reaction ();
	r0->setReversible(false);
	r0->initForwardRateRandomly();
	r0->addModifier(ncomplex);
	r0->addProduct(exProt);
	rlist->push_back(r0);

	//	create reaction 1, binding/unbinding between protein and gene or gene/protein complex
	Reaction* r1 = new Reaction ();
	r1->setReversible(true);
	r1->initForwardRateRandomly();
	r1->initReverseRateRandomly();
	r1->addReactant(nodes[opIndex1]);
	r1->addReactant(nodex[opIndex2]);
	r1->addProduct (ncomplex);
	rlist->push_back(r1);

	return;
}

void Cell::add_mut_postmod () {
	//	a post modification is add
	
    srand(time(NULL));
	if((double)rand()/RAND_MAX <= 0.5) {//	a single protein case


	}
	else {//	two proteins (protein complex) case
		
		int num1 = 0;
		int num2 = 0;
		vector<int> protIndice;
			
		std::vector<Nodes*>::iterator iter = nodes.begin();
		while (iter != nodes.end()) {
			if((*iter)->getNtype() == 1 || (*iter)->isProteinComplex()) { 
				num1++;
				protIndice.push_back(num2);
			}
			num2++;
			iter ++;
		}

		if(!num1) return;	//	no protein or protein complex
		
		int opIndex1 = protIndice[rand()%num1];	//	protein 1
		int opIndex2 = protIndice[rand()%num1];	//	protein 2
	


		

	}
	
	return;
}
	
