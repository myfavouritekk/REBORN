#include "population.h"
#include <ctime>
#include "consts.h"

int main (void) {

    srand(1);
    //  initialization
    ustc::Population myPop(POPULATION, TOTAL_EVO);
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
        myPop.mut_parameters_simAnneal();
        std::cout << "Finished Evolution: " << TOTAL_EVO - myPop.getEvolution() << std::endl;
    }

    
    //  output
    myPop.sort();
    myPop.output();

    return 0;
}
