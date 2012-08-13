//
//  biobrick.cpp
//  PNAS04Team
//
//  Created by Jack Kang on 12-8-11.
//
//

#include "biobrick.h"
#include <fstream>
#include <sstream>

namespace ustc {
    
    
BioBrick::BioBrick(const std::string _geneName, const std::string _promoterName):hasData(true){
    itsGeneName = _geneName;
    itsPromoterName = _promoterName;
    std::ifstream geneFile, promoterFile;
    std::stringstream geneFileName, promoterFileName;
    geneFileName << GENES_PATH << _geneName;
    promoterFileName << PROMOTERS_PATH << _promoterName;
    
    geneFile.open(geneFileName.str().c_str());
    promoterFile.open(promoterFileName.str().c_str());

    if (!geneFile.is_open() || !promoterFile.is_open()) {
        hasData = false;
    }
    geneFile.close();
    promoterFile.close();
    
    
}

std::string BioBrick::getGeneName(){
    return itsGeneName;
}

std::string BioBrick::getPromoterName(){
    return itsPromoterName;
}


bool BioBrick::isAvailableInDatabase(){
    return hasData;
}
   

std::string BioBrick::description(){
    
    std::stringstream result;
    
    if (isAvailableInDatabase()) {
        result << "BioBrick" << std::endl
            << "    Promoter: " << itsPromoterName << std::endl
            << "    Gene: " << itsGeneName << std::endl;
    }
    
    return result.str();
}


    
    
    
    
}   //      namespace ustc