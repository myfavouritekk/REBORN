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
    database >> numOfRegulators >> numOfRegulatees;
    
    //      allocate arrays to store database
    regulatorNames = new std::string[numOfRegulators];
    regulateeNames = new std::string[numOfRegulatees];
    for (int i = 0; i < numOfRegulators; i++) {
        database >> regulatorNames[i];
    }
    for (int i = 0; i < numOfRegulatees; i++) {
        database >> regulateeNames[i];
    }
    
    regulatoryMatix = new int*[numOfRegulatees];
    for (int i = 0; i < numOfRegulatees; i++) {
        regulatoryMatix[i] = new int[numOfRegulators];
        for (int j = 0; j < numOfRegulators; j++) {
            database >> regulatoryMatix[i][j];
        }
    }
    
    database.close();
}
    
    
}   //      namespace ustc
