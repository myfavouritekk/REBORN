#ifndef CELL_H
#define CELL_H

#include "node.h"
#include "reaction.h"

class Cell {

    public:

        Cell* mutation ();
        double getScore (const scoreFunc&, const bool&);

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
		
        //  Global Array of Coefficients
        std::vector<Node*> nodes;
        std::vector<Reaction*> rlist;
};

#endif
