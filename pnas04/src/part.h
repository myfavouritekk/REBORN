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

enum part_Type{PROMOTER, RBS, GENE_PART, TERMINATOR};
    
class Part{
  
private:
    
    part_Type pType;
    std::string id;
    std::string name;
    
public:
    
    //  constructor
    Part(
            const part_Type& _pType,
            const std::string& _id,
            const std::string& _name
            );
    
    
    part_Type getPtype();
    std::string getId();
    std::string getName();
    
};  // class Part



}   //  namespace ustc

#endif /* defined(__PNAS04Team__part__) */
