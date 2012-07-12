#include "population.h"
#include <ctime>


int main (void) {

    srand((unsigned int)time(NULL));
    //  initialization
    // initialize several Populations and run simultaneously
    
    Population** myPops = new Population*[5];
    for (int i = 0; i < 5; i++) {
        myPops[i] = new Population(20);
        myPops[i] -> init();
    }
    

    //  evolution
    while (!myPops[0] -> isTerminate ()) {
        for (int i = 0 ; i < 5; i++) {
            myPops[i] -> growth ();
            myPops[i] -> selection ();

        }
    }

    //myPops.classification();
    
    //  output
    
    for (int i = 0; i < 5; i++) {
        myPops[i] -> output();

    }
    //myPops.genXMLFormat();

    return 0;
}
