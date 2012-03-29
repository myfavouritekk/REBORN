#include "node.h"


Node::Node(int _nindex, Node* _nleft, Node* _nright)
    :nindex(_nindex)
{
	//	merge N left and N right
	
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
			components[0] = _left->getNode(0);
		}
		else{
			stderr << "Error: more than one genes appear in a complex!" << endl;
			std::exit (1);
		}
	}

	//	proteins
	int lsize = _nleft->getNsize ();
	int rsize = _nright->getNsize ();

	for(int ileft = 0; ileft < lsize; ileft++) {
		components.push_back(_nleft->getNode(ileft));
	}
	for(int iright = 0; iright < rsize; iright++) {
		components.push_back(_nright->getNode(iright));
	}

	//	sort
	components.sort ();
	
	//	assign ntype
	if(components[0] == NULL) {
		if(components.size() > 2) {//	protein complex
			ntype = 3;
		}
		else ntype = 1;	//	protein
	}
	else {
		if(components.size() > 1) {//	gene/protein complex
			ntype = 2;
		}
		else ntype = 0;	//	gene
	}

	//	assign representations
	nstring = write();		
}

Node::~Node()
{}

bool Node::operator==(const Node& n1) const {
	return nstring == n1.nstring;
}

int Node::getNindex () {
	return nindex;
}

int Node::getNtype () {
    return ntype;
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

int Node::getNsize () {
	return components.size ();
}

string Node::write () {

	string nodestr('');
	if(components[0] != NULL) {//gene
		nodestr .= components[0].getNodeString().":";
	}

	std::vector<Node*>::iterator iter = components.begin ();
	iter++;
	while (iter != components.end ()) {//proteins
		nodestr .= (*iter)->getNodeString().":";
		iter++;
	}

	return nodestr.substr(0,nodestr.length-1);
}
