#include "population.h"
#include <ctime>


int main (void) {

    srand((unsigned int)time(NULL));
//  initialization
    Population mypop;
    mypop.init ();

//  evolution
    while (!mypop.isTerminate ()) {
        mypop.growth ();
        mypop.selection ();
    }

//  output
    mypop.output();

    return 0;
}
