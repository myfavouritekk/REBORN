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

        double currScore;

        //  Global Array of Coefficients
        std::vector<Node*> nodes;
        std::vector<Reaction*> rlist;
};

#endif
