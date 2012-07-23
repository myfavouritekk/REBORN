#include "reaction.h"
namespace ustc{
Reaction::Reaction ()
{}

Reaction::Reaction(int _rtype):rtype(_rtype)
{}

//copy constructor
Reaction::Reaction(Reaction &aReaction){
    rtype = aReaction.rtype;
    reversible = aReaction.reversible;
    forwardRate = aReaction.forwardRate;
    reverseRate = aReaction.reverseRate;
    
    /*here do not implement vectors in the reaction containing
     *reactants, modifiers and products,
     *leave them for cell to implement, because the relationships between
     *nodes in vector nodes and reactions in rlist need to be
     *implemented in the new cell
     */
}

Reaction::~Reaction ()
{}

bool Reaction::operator==(const Reaction& r1) const {
	if(rtype != r1.rtype) return false;
	//assert(reversible == r1.reversible);
    
    std::vector<Node*>::const_iterator iter, iter1;
    std::vector<int> indice, indice1;

	/*********************
	 * compare reactants
	 * ******************/
	iter = reactants.begin();
	iter1 = r1.reactants.begin();
	while(iter != reactants.end()) {
		indice.push_back((*iter)->getNindex());
		iter++;
	}
	std::sort(indice.begin(), indice.end());

	while(iter1 != r1.reactants.end()) {
		indice1.push_back((*iter1)->getNindex());
		iter1++;
	}
	std::sort(indice1.begin(), indice1.end()); 

	if(indice != indice1) return false;
	indice.clear ();      //clear members in indice
	indice1.clear ();     //clear members in indice1

	/*********************
	 * compare modifiers
	 * ******************/
	iter = modifiers.begin();
	iter1 = r1.modifiers.begin();
	while(iter != modifiers.end()) {
		indice.push_back((*iter)->getNindex());
		iter++;
	}
	std::sort(indice.begin(), indice.end());

	while(iter1 != r1.modifiers.end()) {
		indice1.push_back((*iter1)->getNindex());
		iter1++;
	}
	std::sort(indice1.begin(), indice1.end()); 

	if(indice != indice1) return false;
	indice.clear ();   //clear members in indice
	indice1.clear ();  //clear members in indice1

	/*********************
	 * compare products
	 * ******************/
	iter = products.begin();
	iter1 = r1.products.begin();
	while(iter != products.end()) {
		indice.push_back((*iter)->getNindex());
		iter++;
	}
	std::sort(indice.begin(), indice.end());

	while(iter1 != r1.products.end()) {
		indice1.push_back((*iter1)->getNindex());
		iter1++;
	}
	std::sort(indice1.begin(), indice1.end()); 

	if(indice != indice1) return false;
	indice.clear ();         //clear members in indice
	indice1.clear ();       //clear members in indice1
	//	else
	return true;
}

int Reaction::getRtype () {
	return rtype;
}

bool Reaction::isReversible () {
	return reversible;
}

void Reaction::setReversible(const bool& rev) {
	reversible = rev;
    if (!rev) {
        reverseRate = 0.;
    }
	return;
}

void Reaction::modifyForwardRate () {
    //srand((unsigned int)time(NULL));
    double rn = (double)rand() / RAND_MAX;
    rn = rn * ((1 / RATE_MULTI) - RATE_MULTI) + RATE_MULTI; // rn is uniform in [0.8, 1.25]
    
    forwardRate = forwardRate * rn;
    /*
    if ((rtype == 1 || rtype == 0) && forwardRate < 0.01) {//dagaration type
        forwardRate = 0.01;
    }*/
    return;
}

void Reaction::modifyReverseRate () {
    //srand((unsigned int)time(NULL));
    double rn = (double)rand()/RAND_MAX;
    
    rn = rn * ((1 / RATE_MULTI) - RATE_MULTI) + RATE_MULTI;// rn is uniform in [0.8, 1.25]
    
    reverseRate = reverseRate * rn;
    return;
}

void Reaction::initForwardRateRandomly () {//0-1
    //srand((unsigned int)time(NULL));
    forwardRate = (double)rand()/RAND_MAX;
    /*
    if ((rtype == 1 || rtype == 0) && forwardRate < 0.01) {//dagaration type
        forwardRate = 0.01;
    }*/
	return;
}

void Reaction::initReverseRateRandomly () {//0-1
    //srand((unsigned int)time(NULL));
    reverseRate = (double)rand()/RAND_MAX;
	return;
}


void Reaction::setForwardRate(double rate){
	forwardRate = rate;
}


double Reaction::getForwardRate(){
    return forwardRate;
}

double Reaction::getReverseRate(){
    return reverseRate;
}

void Reaction::addReactant( Node* sr) {
	reactants.push_back(sr);
	return;
}

void Reaction::addModifier( Node* sm) {
	modifiers.push_back(sm);
	return;
}

void Reaction::addProduct ( Node* sp) {
	products.push_back(sp);
	return;
}

Node* Reaction::getReactant (const int& ir) {
	if(ir < 0 || ir >= reactants.size()) return NULL;
	else return reactants[ir];
}

Node* Reaction::getModifier (const int& im) {
	if(im < 0 || im >= modifiers.size()) return NULL;
	else return modifiers[im];
}

Node* Reaction::getProduct  (const int& ip) {
	if(ip < 0 || ip >= products.size()) return NULL;
	else return products[ip];
}

bool Reaction::containNode(Node* aNode){
    if (containNodeAsModifier(aNode) || containNodeAsReactant(aNode) || containNodeAsProduct(aNode)) {
        return true;
    }
    return false;
}

int Reaction::containNodeAsReactant(Node* aNode){
    std::vector<Node*>::iterator iter = reactants.begin();
    std::vector<Node*>::iterator iter_end = reactants.end();
    int counts = 0;
    while (iter != iter_end) {
        if ((*iter) == aNode) {
            counts++;
        }
        iter++;
    }
    return counts;
}



int Reaction::containNodeAsModifier(Node* aNode){
    std::vector<Node*>::iterator iter = modifiers.begin();
    std::vector<Node*>::iterator iter_end = modifiers.end();
    int counts = 0;
    while (iter != iter_end) {
        if ((*iter) == aNode) {
            counts++;
        }
        iter++;
    }
    return counts;
}

int Reaction::containNodeAsProduct(Node* aNode){
    std::vector<Node*>::iterator iter = products.begin();
    std::vector<Node*>::iterator iter_end = products.end();
    int counts = 0;
    while (iter != iter_end) {
        if ((*iter) == aNode) {
            counts++;
        }
        iter++;
    }
    return counts;
}


int Reaction::getReactantsSize(){
    return reactants.size();
}


int Reaction::getModifiersSize(){
    return modifiers.size();
}


int Reaction::getProductsSize(){
    return products.size();
}

std::vector<Node*>* Reaction::getReactantsVector(){
    return &reactants;
}

std::vector<Node*>* Reaction::getModifiersVector(){
    return &modifiers;
}

std::vector<Node*>* Reaction::getProductsVector(){
    return &products;
}


//output method
void Reaction::description(int reactionIndex){
    
    std::cout << "  Reaction " << reactionIndex << ":" << std::endl;
    
    //Print Type:
    std::cout << "   Type: " << "Type " << rtype << "  ";
    switch (rtype) {
        case 0:
            std::cout << "Transcription";
            break;
        case 1:
            std::cout << "Protein Degradation";
            break;
        case 2:
            std::cout << "Binding";
            break;
        case 3:
            std::cout << "Modification";
            break;
        case 4:
            std::cout << "Partial Degradation";
            break;
        case 5:
            std::cout << "Dimerization";
            break;
        case 6:
            std::cout << "Catalytic Degradation";
            break;
        case 7:
            std::cout << "Partial Catalytic Degradation";
            break;
        default:
            break;
    }
    std::cout << std::endl;
    
    //Print reactant, modifiers and products
    std::cout << "   Reactants: ";
    vector<Node*>::iterator iter = reactants.begin();
    while (iter != reactants.end()) {
        std::cout << (*iter)->getNstring() << "\t";
        iter++;
    }
    std::cout << std::endl;
    std::cout << "   Modifiers: ";
    iter = modifiers.begin();
    while (iter != modifiers.end()) {
        std::cout << (*iter)->getNstring() << "\t";
        iter++;
    }
    std::cout << std::endl;
    std::cout << "   Products: ";
    iter = products.begin();
    while (iter != products.end()) {
        std::cout << (*iter)->getNstring() << "\t";
        iter++;
    }
    std::cout << std::endl;
    
    //Print kinect rates
    std::cout << "   Forward rate: " << forwardRate << std::endl;
    std::cout << "   Reverse rate: " << reverseRate << std::endl << std::endl;
    
    
}

}//namespace ustc


