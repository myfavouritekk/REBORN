#ifndef NODE_H
#define NODE_H
#include "cstring"
#include <string>
#include <vector>
using std::string;
using std::vector;

/*
 * Each node represents a species in a bio-reaction network. A species 
 * could be a gene, a protein or a complex composed of gene and proteins or proteins only.
 * We assign each node a specific string having the pattern a:A:B:C:... with
 * the first lower case representing the gene and following upper cases the proteins.
 * Each one of the six types, Inducer 1, Gene 2, Protein 3, Inducer/Protein Complex 4, Gene/Protein Complex 5, Pure Protein Complex 6, 
 * has one associated number, 1, 2, 3, 4, 5, 6, respectively.
 */


class Node {

    public:

		Node (int _nindex, int _ntype); //	basic node constructor
		Node (int _nindex, Node* _nleft, Node* _nright); //complex node constructor
        Node (int _nindex, int _ntype, Node* _nleft, Node* _nright);//complex node constructor of certain type
        ~Node ();   //deconstructor
		
		int getNtype ();	//  return node type 1-6
		int getNindex ();	//	return node index
		int getNsize ();	//	return size of components
		Node* getNode (const int& index);	//	return indexed Node 		
		string getNstring ();	//	return node string

		//	check equality with given node
		bool operator==(const Node& n1) const;

		//	assign representations
		string write ();
    
    
        //for gene/protein complex
        Node* extractFirstGene();

    private:

		int nindex;	//index in the node vector in the cell
        int ntype;	//node type
		string nstring;	//	node string
		vector<Node*> components;	//	nodes constituting the complex, the first component, components[0], is a gene or nil
};

#endif
