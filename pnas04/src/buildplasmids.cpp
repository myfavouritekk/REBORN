//
//  buildplasmids.cpp
//  PNAS04Team
//
//  Created by Jack Kang on 12-8-3.
//
//

#include "buildplasmids.h"
#include <fstream>

namespace ustc {
 

BuildPlasmids::BuildPlasmids()
{}


BuildPlasmids::~BuildPlasmids()
{}


void BuildPlasmids::buildProcess(){
    
    //      loading data from database
    loadDatabase();
    
    ustc::Plasmid** plasmids = new Plasmid*[NUM_SBMLMODEL];
    for (int i = 0; i < NUM_SBMLMODEL; i++) {
        plasmids[i] = new Plasmid();
        plasmids[i] -> readMotifs(i);
    }
    
}
    

void BuildPlasmids::loadDatabase(){
    std::ifstream database;
    database.open("USTC_SOFTWARE_PARTS_DATA.txt");
    database >> numOfGeneParts >> numOfGeneParts;
    
    //      allocate arrays to store database
    geneNames = new std::string[numOfGeneParts];
    promoterNames = new std::string[numOfPromoterParts];
    for (int i = 0; i < numOfGeneParts; i++) {
        database >> geneNames[i];
    }
    for (int i = 0; i < numOfPromoterParts; i++) {
        database >> promoterNames[i];
    }
    
    regulatoryMatix = new int*[numOfPromoterParts];
    for (int i = 0; i < numOfPromoterParts; i++) {
        regulatoryMatix[i] = new int[numOfGeneParts];
        for (int j = 0; j < numOfGeneParts; j++) {
            database >> regulatoryMatix[i][j];
        }
    }
    
    database.close();
}
    
    
}   //      namespace ustc