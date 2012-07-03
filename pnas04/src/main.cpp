#include "population.h"
#include <ctime>


int main (void) {

    srand((unsigned int)time(NULL));
    //  initialization
    Population myPop;
    myPop.init ();

    //  evolution
    while (!myPop.isTerminate ()) {
        myPop.growth ();
        myPop.selection ();
    }

    //  output
    myPop.output();
    myPop.genXMLFormat();

    return 0;
}
