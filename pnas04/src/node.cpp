#include "node.h"



void quickSort(vector<Node*> _components, int num){
    
    int numLess = 0, numGreater = 0;    //to store lengths of two subgroups
    vector<Node*> less;                 //Node** less = new Node*[num];
    vector<Node*> greater;              //Node** greater = new Node*[num];
    
    //end of recursion, no need to sort
    if(num <= 1) return;
    
    vector<Node*>::iterator iter = _components.begin();
    iter++;
    if (_components.size() !=1 ) {
        if ((*iter)->getNindex() > (*_components.begin())->getNindex()) {
            greater.push_back(*iter);
            iter = _components.erase(iter);
            numGreater++;
        }
        else {
            less.push_back(*iter);
            iter = _components.erase(iter);
            numLess++;
        }
    }
    
    
    //recursion
    quickSort(less, numLess);
    quickSort(greater, numGreater);

    
    //final assignments
    //first, for less members
    std::vector<Node*>::iterator iter1 = less.end()-1;
    if (iter1 != less.begin()) {
        _components.insert(_components.begin(),*iter1);
        iter1--;
    }
    _components.insert(_components.begin(),*iter1);
    //second, for greater members
    std::vector<Node*>::iterator iter2 = greater.begin();
    if (iter2 != greater.end()) {
        _components.push_back(*iter2);
        iter2++;
    }
}



void sort(std::vector<Node*> _components, int num_of_members[3]){
    
    //count each components and put the numbers in members
    for (int i = 0; i < 3; i++) {
        num_of_members[i] = 0;
    }
    vector<Node*>::iterator iter = _components.begin();
    if (*iter != NULL) {
        num_of_members[2-1] = 1;                //have one gene
        iter++;
    }
    vector<Node*>::iterator iter_end = _components.end();
    do {
        num_of_members[(*iter)->getNtype()-1]++;
    }
    while(iter++ != iter_end);
    
    //sort the vector by components' indice
    int num = _components.size()-1;

    //delete the first one
    //either gene or NULL, should not be sorted
    Node* temp = *_components.begin();
    _components.erase(_components.begin());
    //sort the rest
    quickSort(_components, num);
    //put the first one back
    _components.insert(_components.begin(), temp);

}


//basic node constructor
Node::Node(int _nindex, int _ntype):nindex(_nindex),ntype(_ntype),components(NULL){
    
    
    if (ntype == 1) {               //inducer
        std::stringstream ss;
        ss << nindex;
        nstring = "indu" + ss.str();
    }   
    if (ntype == 2) {               //gene
        std::stringstream ss;
        ss << nindex;
        nstring = "g" + ss.str();
    }
    
    if (ntype == 3) {               //protein
        std::stringstream ss;
        ss << nindex;
        nstring = "P" + ss.str();
    }    
}


//complex node constructor
Node::Node(int _nindex, Node* _nleft, Node* _nright)
:nindex(_nindex)
{
	//	merge N_left and N_right
	
	//	gene
	if(_nleft->getNode(0) == NULL) {
		if(_nright->getNode(0) == NULL) {
			components[0] = NULL;
		}
		else {
			components[0] = _nright->getNode(0);
		}
	}
	else {
		if(_nright->getNode(0) == NULL) {
			components[0] = _nleft->getNode(0);
		}
		else{
            std::cerr << "Error: more than one genes appear in a complex!" << std::endl;
			std::exit (1);
		}
	}
    //the first component is set whether its a gene or NULL
    
	//	proteins
	int lsize = _nleft->getNsize ();
	int rsize = _nright->getNsize ();
    
	for(int ileft = 1; ileft < lsize; ileft++) {
		components.push_back(_nleft->getNode(ileft));
	}
	for(int iright = 1; iright < rsize; iright++) {
		components.push_back(_nright->getNode(iright));
	}
    
	//	sort
	//components.sort ();
    int members[3]={0};
    sort(components,members);
    
    ntype = 0;
    //assign ntype
    if (members[1]) {               //gene
        if (members[2]) {           //protein
            ntype = 5;              //gene/protein complex
        }else {
            ntype = 2;              //only gene
        }
    }
    else if (members[2]) {          //protein
            if (members[0]) {       //inducer
                ntype = 4;          //inducer/protein complex
            }else if(members[2]>=2) {
                    ntype = 6;      //pure protein complex    
                }else {
                    ntype = 3;      //only protein
                }
         }else {
            ntype = 1;              //only inducer
    }
    
    
	//	assign representations
	nstring = write();		
}


//complex node constructon of certain type
Node::Node(int _nindex, int _ntype, Node* _nleft, Node* _nright):ntype(_ntype){
    Node(_nindex, _nleft, _nright);
}


Node::~Node()
{}


int Node::getNindex () {
	return nindex;
}

int Node::getNtype () {
    return ntype;
}
int Node::getNsize () {
	return components.size ();
}

Node* Node::getNode (const int& index) {
	if(index >= 0 && index < components.size()) {
		return components[index];
	}
	else return NULL;
}

string Node::getNstring () {
	return nstring;
}

bool Node::operator==(const Node& n1) const {
	return nstring == n1.nstring;
}

string Node::write () {

	string nodestr = NULL;
	if(components[0] != NULL) {                     //gene
		nodestr += components[0]->getNstring()+":";
	}

	std::vector<Node*>::iterator iter = components.begin ();
	iter++;                                         //components[0] is for gene
	while (iter != components.end ()) {             //proteins
		nodestr = nodestr + (*iter)->getNstring() + ":";
		iter++;
	}

	return nodestr.substr(0,nodestr.length()-1);    //delete last ":"
}


Node* Node::extractFirstGene(){
    return components[0];
}
