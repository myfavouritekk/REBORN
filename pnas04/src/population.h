#ifndef POPULATION_H
#define POPULATION_H

#include "cell.h"

class Population {

    public:

        void growth ();
        void selection ();
   
    private:
        
        Cell* cells;
        double* targetData;
        scoreFunction sfunc;
};

#endif
