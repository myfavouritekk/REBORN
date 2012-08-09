//
//  operon.h
//  PNAS04Team
//
//  Created by Jack Kang on 12-8-7.
//
//

#ifndef __PNAS04Team__operon__
#define __PNAS04Team__operon__

#include <iostream>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <vector>



namespace ustc {
    
    
class Operon{
    
public:
    //      constructor and destructor
    Operon(const std::string _name);
    ~Operon();
    
    //      description method
    std::string description();
    
    //      whether this Operon has data in the database
    bool isAvailableInDatabase();
    
    
    //      get components numbers
    int getNumGenes();
    int getNumPromoters();
    int getNumTerminators();
    
private:
    std::string name;
    std::vector<std::string> genes;
    std::vector<std::string> promoters;
    std::vector<std::string> terminators;
    
    bool hasData;
    
};
    
    
    
    
}   //      namespace ustc


#endif /* defined(__PNAS04Team__operon__) */
