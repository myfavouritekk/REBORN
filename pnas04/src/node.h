#ifndef NODE_H
#define NODE_H

class Node {

    public:

        Node (int _ntype=0, Node* _nleft=NULL, Node* _nright=NULL); //constructor
        ~Node ();   //deconstructor

        int getNtype (); //return node type
        Node* getNleft (); //return component 1
        Node* getNright (); //return component 2

    private:

        int ntype; //0 is gene, 1 is protein, 2 is complex
        Node* nleft, nright;
};

#endif
