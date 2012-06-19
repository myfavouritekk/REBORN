#include "reaction.h"
#include <time.h>
#include <algorithm>
#include <vector>

Reaction::Reaction () 
{}

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
		indice1.push_back((*iter)->getNindex());
		iter1++;
	}
	std::sort(indice1.begin(), indice1.end()); 

	if(indice != indice1) return false;
	indice.empty ();
	indice1.empty ();

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
		indice1.push_back((*iter)->getNindex());
		iter1++;
	}
	std::sort(indice1.begin(), indice1.end()); 

	if(indice != indice1) return false;
	indice.empty ();
	indice1.empty ();

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
		indice1.push_back((*iter)->getNindex());
		iter1++;
	}
	std::sort(indice1.begin(), indice1.end()); 

	if(indice != indice1) return false;

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
	return;
}

void Reaction::modifyForwardRate () {
    srand(time(NULL));
    double rn = (double)rand() / RAND_MAX * 2.0;
    forwardRate = forwardRate*rn;
    return;
}

void Reaction::modifyReverseRate () {
    srand(time(NULL));
    double rn = (double)rand()/RAND_MAX*2.0;
    reverseRate = reverseRate*rn;
    return;
}

void Reaction::initForwardRateRandomly () {//0-1
    srand(time(NULL));
    forwardRate = (double)rand()/RAND_MAX;
	return;
}

void Reaction::initReverseRateRandomly () {//0-1
    srand(time(NULL));
    reverseRate = (double)rand()/RAND_MAX;
	return;
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
