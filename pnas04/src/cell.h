#ifndef CELL_H
#define CELL_H

#include "node.h"
#include "reaction.h"

class Cell {

    public:

        Cell* mutation ();
        double getScore (const scoreFunc&, const bool&);

    private:

        double currScore;

        //  Global Array of Coefficients
        double** coefs;
        Node** nodes;
        Reaction* ReactionList;
};

#endif
