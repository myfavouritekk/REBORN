#include "cell.h"

void Cell::mut_deg_prot () {
    //  the degradation rate of a protein is modified

    int num1 = 0;
    int num2 = 0;
    vector<int> indice;

    //  scan <rlist> and choose one degradation reaction (protein)
    std::vector<Reaction*>::iterator iter = rlist.begin();
    while (iter != rlist.end()) {
        if(iter->getType() == 1) { //   #1 is protein degradation
            num1++;
            indice.push_back(num2);
        }
        num2++;
        iter ++;
    }

    srand(time(NULL));
    int opIndex = indice[rand()%num1];

    //  modify its degradation rate
    Reaction* currR = rlist[opIndex];
    currR->modifyForwardRate ();

    return;
}

