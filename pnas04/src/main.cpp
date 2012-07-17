#include "population.h"
#include <ctime>
#define TOTAL_EVO 750

int main (void) {

    srand((unsigned int)time(NULL));
    //  initialization
    Population myPop(200, TOTAL_EVO);
    myPop.init ();

    int i = 1;
    int sum = 1;
    //  evolution
    while (!myPop.isTerminate ()) {
        if ((TOTAL_EVO + 1 - myPop.getEvolution()) % (50 * sum) == 0) {
            //sort
            myPop.sort();
            myPop.output();
            myPop.mutation();
            i++;
            sum += i;
        }
        myPop.mut_parameters();
        std::cout << "Finished Evolution: " << TOTAL_EVO - myPop.getEvolution() << std::endl;
    }

    
    //  output
    myPop.sort();
    myPop.output();

    return 0;
}
