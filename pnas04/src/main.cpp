#include "population.h"
#include <ctime>
#define TOTAL_EVO 1000

int main (void) {

    srand((unsigned int)time(NULL));
    //  initialization
    Population myPop(100, TOTAL_EVO);
    myPop.init ();

    //  evolution
    while (!myPop.isTerminate ()) {
        if (myPop.getEvolution() % 200 == 0) {
            //sort
            myPop.sort();
            myPop.output();
            myPop.mutation();
        }
        myPop.mut_parameters();
        std::cout << "Finished Evolution: " << TOTAL_EVO - myPop.getEvolution() << std::endl;
    }

    
    //  output
    myPop.sort();
    myPop.output();

    return 0;
}
