#include "population.h"

int main (void) {

//  initialization
    Population mypop;
    mypop.init ();

//  evolution
    while (mypop.isTerminate ()) {
        mypop.growth ();
        mypop.selection ();
    }

//  output
    mypop.genTikzFormat ();

    return 0;
}
