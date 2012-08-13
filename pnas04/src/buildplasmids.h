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
    
    
    //      build process based on finding operon-operon relationships
    void buildUsingOperons();
    
    //      build process based on finding promoter-gene relationships
    void buildUsingBioBricks();
    
    //
    //      data loaded in from the external database
    //
    int numOfRegulators;
    int numOfRegulatees;
    std::string* regulatorNames;
    std::string* regulateeNames;
    int** wholeRegulatoryMatrixInDataBase;
    
    
};
    
    
}   //      namespace ustc




#endif /* defined(__PNAS04Team__buildplasmids__) */
