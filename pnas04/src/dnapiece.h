//
//  dnapiece.h
//  PNAS04Team
//
//  Created by Jack Kang on 12-8-1.
//
//

#ifndef __PNAS04Team__dnapiece__
#define __PNAS04Team__dnapiece__

#include <iostream>
#include "part.h"
#include "consts.h"

namespace ustc {

class DNAPiece{
    
private:
    
    Part* promoter;
    Part* rbs;
    Part* gene;
    Part* terminator;
    
public:
    
    Part* getPromoter();
    Part* getRbs();
    Part* getGene();
    Part* getTerminator();
   
    
};  //  class DNAPiece
    
    
    
    
}// namespace ustc



#endif /* defined(__PNAS04Team__dnapiece__) */
