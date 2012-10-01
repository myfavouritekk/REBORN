//
//  networkinference.h
//  PNAS04Team
//
//  Created by Jack Kang on 12-8-3.
//
//

#ifndef __PNAS04Team__networkinference__
#define __PNAS04Team__networkinference__

#include <iostream>
#include "population.h"
#include "consts.h"

namespace ustc {
    
    
class NetworkInference{
    
public:
    
    NetworkInference();
    ~NetworkInference();
    
    
    //
    //      the how process of generating genetic regulatory networks
    //      based on the input time courses
    //
    void reverseEngineering(std::string);
    
private:
    
    //
    //      methods to help network inference process
    //
    //
    void askInformation(ustc::Population* targetPop);
};

    
}   //      namespace ustc




#endif /* defined(__PNAS04Team__networkinference__) */
