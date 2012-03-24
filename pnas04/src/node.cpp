#include "node.h"


Node::Node(int _ntype, Node* _nleft, Node* _nright)
    :ntype(_ntype), nleft(_nleft), nright(_nright)
{}

Node::~Node()
{}

Node* Node::getNright () {
    if(nright != NULL) {
        return nright;
    }
    else {
        return NULL;
    }
}


Node* Node::getNleft(){
    if(nleft != NULL){
        return nleft;
    }
    else{
        return NULL;
    }
}


int Node::getNtype(){
    return ntype;
}

