#include "population.h"
#include <ctime>


int main (void) {

    srand((unsigned int)time(NULL));
    //  initialization
    Population myPop(100);
    myPop.init ();

    //  evolution
    while (!myPop.isTerminate ()) {
        if (myPop.getEvolution() % 200 == 0) {
            //sort
            myPop.sort();
            //  output
            myPop.output();
            myPop.mutation();
        }
        myPop.mut_parameters();
        std::cout << "Finished Evolution: " << 999 - myPop.getEvolution() << std::endl;
    }

    //myPop.classification();
    
    //  output
    myPop.sort();
    myPop.output();
    //myPop.genXMLFormat();

    return 0;
}
