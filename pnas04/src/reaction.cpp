#include "reaction.h"

void Reaction::modifyForwardRate () {
    srand(time(NULL));
    double rn = (double)rand()/RAND_MAX*2.0;
    forwardRate = forwardRate*rn;
    return;
}

void Reaction::modifyReverseRate () {
    srand(time(NULL));
    double rn = (double)rand()/RAND_MAX*2.0;
    reverseRate = reverseRate*rn;
    return;
}
