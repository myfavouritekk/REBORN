#include "node.h"


Node::Node(int _nindex, int _ntype, Node* _nleft, Node* _nright)
    :nindex(_nindex), ntype(_ntype), nleft(_nleft), nright(_nright)
{}

Node::~Node()
{}

bool Node::operator==(const Node& n1) const {
	return ntype == n1.ntype;
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
