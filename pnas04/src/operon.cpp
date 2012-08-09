//
//  operon.cpp
//  PNAS04Team
//
//  Created by Jack Kang on 12-8-7.
//
//

#include "operon.h"
#include "consts.h"


namespace ustc {
    
    
Operon::Operon(const std::string _name):hasData(true){
    
    
    //      using the identifier as the name of a file and query other information
    
    std::ifstream tuDataFile;
    std::stringstream tuDataFileName;
    tuDataFileName << OPERONS_PATH << _name << ".txt";
    tuDataFile.open(tuDataFileName.str().c_str());
    
    //      in case of no file for this operon
    if (!tuDataFile) {
        
        name = _name;
        
        hasData = false;    //      to show not available in database
        
        return;
        
    }
    
    
    //      initialize the members
    std::string geneNames;
    int numGene;
    int numPromoter;
    int numTerminator;
    tuDataFile >> name >> numGene >> numPromoter >> numTerminator;

    //      adding genes
    for (int i = 0; i < numGene; i++) {
        std::string aGeneName;
        tuDataFile >> aGeneName;
        genes.push_back(aGeneName);
    }
    
    //      addding promoters
    for (int i = 0; i < numPromoter; i++) {
        std::string aPromoterName;
        tuDataFile >> aPromoterName;
        promoters.push_back(aPromoterName);
    }
    
    //      addding terminators
    for (int i = 0; i < numTerminator; i++) {
        std::string aTerminatorName;
        tuDataFile >> aTerminatorName;
        promoters.push_back(aTerminatorName);
    }

    
    tuDataFile.close();
}



std::string Operon::description(){

    std::stringstream result;
    
    if (isAvailableInDatabase()) {// has data
        result << "Operon: " << name << std::endl;
        
        //      print genes' information
        result << "    Genes: ";
        int numberOfGenes = (int)genes.size();
        for (int i = 0; i < numberOfGenes; i++) {
            result << genes[i] << "\t";
        }
        result << std::endl;
        
        //      print promoters' information
        result << "    Promoters: ";
        int numberOfPromoters = (int)promoters.size();
        for (int i = 0; i < numberOfPromoters; i++) {
            result << promoters[i] << "\t";
        }
        result << std::endl;
        
        //      print terminators' information
        result << "    Terminators: ";
        int numberOfTerminators = (int)terminators.size();
        for (int i = 0; i < numberOfTerminators; i++) {
            result << terminators[i] << "\t";
        }
        result << std::endl;
        
        
        return result.str();
    }
    
    result << "Operon: " << name << std::endl;
    result << "    Not available in the database,\n    Please update the database and try again.\n";
    return result.str();
}
    
    
int Operon::getNumGenes(){
    return (int)genes.size();
}

int Operon::getNumPromoters(){
    return (int)promoters.size();
}

int Operon::getNumTerminators(){
    return (int)terminators.size();
}
   
bool Operon::isAvailableInDatabase(){
    return hasData;
}
    
    
}   //      namespace ustc