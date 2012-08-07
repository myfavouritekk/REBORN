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
    Operon(const std::string _identifier);
    ~Operon();
    
    //      description method
    void description();
    
private:
    std::string name;
    std::vector<std::string> genes;
    
    
};
    
    
    
    
}   //      namespace ustc


#endif /* defined(__PNAS04Team__operon__) */
