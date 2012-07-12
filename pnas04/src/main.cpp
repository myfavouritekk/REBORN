#include "population.h"
#include <ctime>


int main (void) {

    srand((unsigned int)time(NULL));
    //  initialization
    Population myPop(20);
    myPop.init ();

    //  evolution
    while (!myPop.isTerminate ()) {
        if (myPop.getEvolution() % 100 == 0) {
            //sort
            myPop.sort();
            //  output
            myPop.output();
            myPop.mutation();
        }
        myPop.mut_parameters();
    }

    //myPop.classification();
    
    //  output
    myPop.sort();
    myPop.output();
    //myPop.genXMLFormat();

    return 0;
}
