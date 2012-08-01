//
//  part.h
//  PNAS04Team
//
//  Created by Jack Kang on 12-8-1.
//
//

#ifndef __PNAS04Team__part__
#define __PNAS04Team__part__

#include <iostream>
#include "consts.h"

namespace ustc {

enum part_Type{PROMOTER, RBS, GENE, TERMINATOR};
    
class Part{
  
private:
    
    part_Type pType;
    
public:
    
    part_Type getPtype();
    
};  // class Part



}   //  namespace ustc

#endif /* defined(__PNAS04Team__part__) */
