#ifndef NODE_H
#define NODE_H

/*
 * Each node represents a species in a bio-reaction network. A species 
 * could be a gene, a protein or a complex composed of gene and proteins or proteins only.
 * We assign each node a specific string having the pattern a:A:B:C:... with
 * the first lower case representing the gene and following upper cases the proteins.
 * Each one of the four types, Gene, Protein, Gene/Protein Complex, Pure Protein Complex,
 * has one associated number, 0, 1, 2, 3, respectively.
 */

class Node {

    public:

		Node (int _nindex, Node* _nleft=NULL, Node* _nright=NULL); //constructor
        ~Node ();   //deconstructor
		
		int getNtype ();	//  return node type
		int getNindex ();	//	return node index
		Node* getNode (const int& index);	//	return indexed Node 		
		Node* getNstring ();	//	return node string

		//	check equality with given node
		bool operator==(const Node& n1) const;

		//	output
		string write ();

    private:

		int nindex;	//index in the node vector
        int ntype;	//node type
		string nstring;	//	node string
		vector<Node*> components;	//	nodes constituting the complex
};

#endif
