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
    
    std::stringstream ss;
    ss << nindex;

    switch (ntype) {
    case 1:
        nstring = "indu" + ss.str();
        break;
    case 2:
        nstring = "g" + ss.str();
        break;
    case 3:
        nstring = "P" + ss.str();
        break;
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



//copy constructor
Node::Node(Node &aNode){
    nstring = aNode.nstring;
    nindex = aNode.nindex;
    ntype = aNode.ntype;
    
    /* do not copy components, leave it for cell to finish,
     * because it needs the relationships with other node,
     * with cannot judge right now.
     */
    
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
        nodestr += components[0]->getNstring() + ":";
    }

    std::vector<Node*>::iterator iter = components.begin ();
    std::vector<Node*>::iterator iter_end = components.end();
    //components[0] is for gene
    while (++iter != iter_end) {                    //proteins
        nodestr += (*iter)->getNstring() + ":";
    }

    return nodestr.substr(0, nodestr.length()-1);   //delete last ":"
}


Node* Node::extractFirstGene(){
    return components[0];
}


//for cell to access components vector
void Node::pushNode(Node* aNode){
    components.push_back(aNode);
}




#include "reaction.h"
/*
 ode:
 input: reactionList is the rlist variable in the cell; y contains the current value of each node, x is current time
 basic concepts: ode has the form of dN/dx = ode(Y,x), y is a vector contains current nodes' values, and node index is the subscript i in y[i]
 */
double Node::ode(std::vector<Reaction*> reactionList, double *y, double t){
    double result = 0;
    std::vector<Reaction*>::iterator iter = reactionList.begin();
    std::vector<Reaction*>::iterator iter_end = reactionList.end();
    
    /* iterate over all the reactions in the reactionList, and every reaction that contains this node will have influence on result */
    while (iter != iter_end) {
        bool containThisNode;
        containThisNode = std::find((*iter)->reactants.begin(), (*iter)->reactants.end(), this) != (*iter)->reactants.end() || std::find((*iter)->modifiers.begin(), (*iter)->modifiers.end(), this) != (*iter)->modifiers.end() || std::find((*iter)->products.begin(), (*iter)->products.end(), this) != (*iter)->products.end();
        
        
        /* current reaction contains this node in reactants, modifiers or products */
        if (containThisNode) {
            int reactionType = (*iter) -> getRtype();
            switch (reactionType) {
                case 0:                 //Transcription of a gene from a free promoter or a bound promoter/protein complex
                    if (!(std::find((*iter)->products.begin(), (*iter)->products.end(), this)==(*iter)->products.end())) {
                        //this node is in prodcuts of type 0 reaction , that is , this node is a protein, and is to be transcripted
                        int DNAIndex;
                        double forwardRate = (*iter)->forwardRate;
                        DNAIndex = (*iter)->modifiers[0]->nindex;//for protein transcription, a DNA is regarded as a modifier rather than a reactant
                        result += forwardRate * y[DNAIndex];
                    }
                    break;
                    
                case 1:                 //degradation of protein
                    if (!(std::find((*iter)->reactants.begin(), (*iter)->reactants.end(), this)==(*iter)->reactants.end())) {
                        //this node is in reactants of type 1 reaction , that is , this node is a protein, and is to degrade
                        int proteinIndex;
                        double forwardRate = (*iter)->forwardRate;
                        proteinIndex = this->nindex;
                        result -= forwardRate * y[proteinIndex];
                    }
                    break;
                    
                    
                default:                //all the rest reaction types could have the same routine
                    
                    /* if this node is in reactants */
                    if (!(std::find((*iter)->reactants.begin(), (*iter)->reactants.end(), this)==(*iter)->reactants.end())) {
                        //this node is a reactant
                        double reactantsMultiply = 1., productsMultiply = 1.;
                        
                        /* calculate reactantsMultiply */
                        std::vector<Node*>::iterator reactantsIter = (*iter)->reactants.begin();
                        std::vector<Node*>::iterator reactantsIter_end = (*iter)->reactants.end();
                        while (reactantsIter != reactantsIter_end) {
                            int index = (*reactantsIter)->nindex;
                            reactantsMultiply *= y[index];
                            reactantsIter++;
                        }
                        
                        /* calculate productsMultiply */
                        std::vector<Node*>::iterator productsIter = (*iter)->products.begin();
                        std::vector<Node*>::iterator productsIter_end = (*iter)->products.end();
                        while (productsIter != productsIter_end) {
                            int index = (*productsIter)->nindex;
                            productsMultiply *= y[index];
                            productsIter++;
                        }
                        
                        double forwardRate = (*iter)->forwardRate;
                        double reverseRate = (*iter)->reverseRate;
                        
                        result = result -  forwardRate * reactantsMultiply + reverseRate * productsMultiply;
                    }
                    
                    
                    /* if this node is in products, 
                     *can't include the following block in else of the previous if statement, 
                     *for this node may be in both reactants and products 
                     */
                    if (!(std::find((*iter)->products.begin(), (*iter)->products.end(), this)==(*iter)->products.end())) {
                            //this node is a reactant
                        double reactantsMultiply = 1., productsMultiply = 1.;
                        
                        /* calculate reactantsMultiply */
                        std::vector<Node*>::iterator reactantsIter = (*iter)->reactants.begin();
                        std::vector<Node*>::iterator reactantsIter_end = (*iter)->reactants.end();
                        while (reactantsIter != reactantsIter_end) {
                            int index = (*reactantsIter)->nindex;
                            reactantsMultiply *= y[index];
                            reactantsIter++;
                        }
                        
                        /* calculate productsMultiply */
                        std::vector<Node*>::iterator productsIter = (*iter)->products.begin();
                        std::vector<Node*>::iterator productsIter_end = (*iter)->products.end();
                        while (productsIter != productsIter_end) {
                            int index = (*productsIter)->nindex;
                            productsMultiply *= y[index];
                            productsIter++;
                        }
                        
                        double forwardRate = (*iter)->forwardRate;
                        double reverseRate = (*iter)->reverseRate;
                        
                        result = result +  forwardRate * reactantsMultiply - reverseRate * productsMultiply;
                    }
                break;
            }
        }
    }
    
    return result;
}



