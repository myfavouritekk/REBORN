#ifndef CELL_H
#define CELL_H

#include "node.h"
#include "reaction.h"

class Cell {

    public:

		//	no default constructor
		Cell (const int& _numind, const int& _numprot);

		//	deconstructor
		~Cell ();

		//	call for six submutations successively
        Cell* mutation ();

    private:

        //  The degradation rate of a protein is modified
        void mut_deg_prot();

		//	A kinetic constant of one reaction is modified
		void mut_kin_const();

		//	A new gene is created
		void mut_add_gene();

		//	A new interaction between a protein and a gene or gene/protein complex is introduced
		void mut_add_regu();

		//	A post transcriptional modification is added
		void mut_add_postmod();

		//	check if a node has been added previously
		bool existsNode (const Node& node);
		
		//	check if a reaction has been added previously
		bool existsReaction (const Reaction& rxn);

        //	Global Storage
        vector<Node*> nodes;
        vector<Reaction*> rlist;	//	document operations made to develop the network
};

#endif
