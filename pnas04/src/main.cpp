#include "population.h"
#include <ctime>
#define TOTAL_EVO 1000

int main (void) {

    srand((unsigned int)time(NULL));
    
    int evolutionForTopology = 10;
    
    //  initialization
    Population myPop(20, TOTAL_EVO);
    myPop.init ();

    //  evolution
    while (!myPop.isTerminate ()) {
        
        if (myPop.getEvolution() <= evolutionForTopology) {
            myPop.mutation();
        }else{
            myPop.mut_parameters();
        }
        std::cout << "Finished Evolution: " << TOTAL_EVO - myPop.getEvolution() << std::endl;
    }

    
    //  output
    myPop.sort();
    myPop.output();

    return 0;
}
