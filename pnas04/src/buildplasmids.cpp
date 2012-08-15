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
    
    
    std::cout << "Which build method do you what?" << std::endl
        << "1. Using Operons as basic elements." << std::endl
        << "2. Using promoters and genes as basic elements." << std::endl;
    int buildChoice;
    std::cin >> buildChoice;
    
    switch (buildChoice) {

        case 2:{
            buildUsingBioBricks();
            break;
        }
        default:{
            
            //==============================================================//
            //      build plasmid based on operon-operon relationships      //
            //==============================================================//
            buildUsingOperons();
            
            break;
        }
    }
}
    

void BuildPlasmids::loadDatabase(){
   
	std::ifstream databaseOfOperons;
	std::ifstream databaseOfGenesAndPromoters;
    std::stringstream databaseNameOfOperons;
	std::stringstream databaseNameOfGenesAndPromoters;

    databaseNameOfOperons << DATABASE_PATH << "USTC_SOFTWARE_PARTS_DATA.txt";
	databaseNameOfGenesAndPromoters << DATABASE_PATH << "USTC_SOFTWARE_BIOBRICKS_DATA.txt";

    databaseOfOperons.open(databaseNameOfOperons.str().c_str());
	databaseOfGenesAndPromoters.open(databaseNameOfGenesAndPromoters.str().c_str());

    if (!databaseOfOperons || !databaseOfGenesAndPromoters) {
        std::cerr << "Error, unable to load database!" << std::endl;
		return;
    }
	
	//load database with operons
    databaseOfOperons >> numOfRegulators >> numOfRegulatees;
    
    //      allocate arrays to store database
    regulatorNames = new std::string[numOfRegulators];
    regulateeNames = new std::string[numOfRegulatees];
    for (int i = 0; i < numOfRegulators; i++) {
        databaseOfOperons >> regulatorNames[i];
    }
    for (int i = 0; i < numOfRegulatees; i++) {
        databaseOfOperons >> regulateeNames[i];
    }
    
    wholeRegulatoryMatrixInDataBase = new int*[numOfRegulatees];
    for (int i = 0; i < numOfRegulatees; i++) {
        wholeRegulatoryMatrixInDataBase[i] = new int[numOfRegulators];
        for (int j = 0; j < numOfRegulators; j++) {
            databaseOfOperons >> wholeRegulatoryMatrixInDataBase[i][j];
        }
    }
    
    databaseOfOperons.close();

	//load database with promoters and genes
	databaseOfGenesAndPromoters >> numOfGenes >> numOfPromoters;

	geneNames = new std::string[numOfGenes];
	promoterNames = new std::string[numOfPromoters];
	for(int i = 0;i < numOfGenes; i++){
		databaseOfGenesAndPromoters >> geneNames[i];
        std::cout << geneNames[i] << "\t";
	}
    
    std::cout << std::endl << std::endl;
    
	for(int i = 0;i < numOfPromoters;i++){
		databaseOfGenesAndPromoters >> promoterNames[i];
        std::cout << promoterNames[i] << "\t";
	}

	wholeRegulatoryMatrixOfGenesAndPromoters = new int*[numOfPromoters];
	for(int i = 0; i < numOfPromoters; i++){
		wholeRegulatoryMatrixOfGenesAndPromoters[i] = new int[numOfGenes];
		for(int j = 0;j < numOfGenes; j++){
			databaseOfGenesAndPromoters >> wholeRegulatoryMatrixOfGenesAndPromoters[i][j];
            std::cout << wholeRegulatoryMatrixOfGenesAndPromoters[i][j] << "\t";
		}
        std::cout << std::endl << std::endl;
	}

	databaseOfGenesAndPromoters.close();

}
    

    
//      build plasmid based on operon-operon relationships
void BuildPlasmids::buildUsingOperons(){
    ustc::Plasmid** plasmids = new Plasmid*[NUM_SBMLMODEL];
    for (int plasmidIndex = 0; plasmidIndex < NUM_SBMLMODEL; plasmidIndex++) {
        plasmids[plasmidIndex] = new Plasmid();
        
        //  read complete regulatory matrix
        plasmids[plasmidIndex] -> readCompleteMatrix(plasmidIndex);
        
        //  find complete regulatory matrix in the database
        plasmids[plasmidIndex] -> findCompleteCandidates(
                                                         numOfRegulatees,
                                                         numOfRegulators,
                                                         regulatorNames,
                                                         regulateeNames,
                                                         (const int**)wholeRegulatoryMatrixInDataBase
                                                         );
        
        //  generate plans to build the plasmid
        plasmids[plasmidIndex] -> generatePlans();
        
        //  output those plans into files
        plasmids[plasmidIndex] -> generatePlanOutputs(plasmidIndex);
        
    }
}
    
    
//      build plasimid based on gene-promoter relationships
void BuildPlasmids::buildUsingBioBricks(){
    
    ustc::Plasmid** plasmids = new Plasmid*[NUM_SBMLMODEL];
    
    for (int plasmidIndex = 0; plasmidIndex < NUM_SBMLMODEL; plasmidIndex++) {
        plasmids[plasmidIndex] = new Plasmid();
        
        //  read complete regulatory matrix
        plasmids[plasmidIndex] -> readCompleteMatrix(plasmidIndex);
        
        //  find complete regulatory matrix in the database
        plasmids[plasmidIndex] -> findCompleteCandidatesUsingBiobricks(
                                                         numOfPromoters,
                                                         numOfGenes,
                                                         regulatorNames,
                                                         regulateeNames,
                                                         (const int**)wholeRegulatoryMatrixOfGenesAndPromoters
                                                         );
        
        //  generate plans to build the plasmid
        //plasmids[plasmidIndex] -> generatePlans();
        
        //  output those plans into files
        //plasmids[plasmidIndex] -> generatePlanOutputs(plasmidIndex);
        
    }
   
}
    
    
}   //      namespace ustc
