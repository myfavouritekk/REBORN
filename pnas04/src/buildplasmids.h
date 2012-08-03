//
//  buildplasmids.h
//  PNAS04Team
//
//  Created by Jack Kang on 12-8-3.
//
//

#ifndef __PNAS04Team__buildplasmids__
#define __PNAS04Team__buildplasmids__

#include <iostream>
#include "plasmid.h"

namespace ustc {
    
    
class BuildPlasmids{

    friend class plasmid;
    
public:
    BuildPlasmids();
    ~BuildPlasmids();
    
    
    //
    //      the whole process of building the plasmids with
    //      genetic networks generated from network inference
    //      process and real parts form partsregistry
    //
    void buildProcess();
    
private:
    
    void loadDatabase();
    
    //
    //      data loaded in from the external database
    //
    int numOfRegulators;
    int numOfRegulatees;
    std::string* regulatorNames;
    std::string* regulateeNames;
    int** regulatoryMatix;
    
    
};
    
    
}   //      namespace ustc




#endif /* defined(__PNAS04Team__buildplasmids__) */
