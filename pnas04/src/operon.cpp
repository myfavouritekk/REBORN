//
//  operon.cpp
//  PNAS04Team
//
//  Created by Jack Kang on 12-8-7.
//
//

#include "operon.h"


namespace ustc {
    
    
Operon::Operon(const std::string _name){
    
    
    //      using the identifier as the name of a file and query other information
    
    std::ifstream tuDataFile;
    std::stringstream tuDataFileName;
    tuDataFileName << _name << ".txt";
    tuDataFile.open(tuDataFileName.str().c_str());
    
    //      initialize the members
    std::string geneNames;
    int numGene;
    tuDataFile >> name >> numGene >> geneNames;
    int* nameIndex = new int[numGene + 1];
    nameIndex[0] = 0;
    for (int i = 0; i < numGene; i++) {
        nameIndex[i + 1] = (int)geneNames.find(",", nameIndex[i] + 1, (int)geneNames.size());
        int start, end;
        if (i == 0) 
            start = 0;
        else
            start = nameIndex[i] + 1;
        if (i == numGene - 1)
            end = (int)geneNames.size();
        else
            end = nameIndex[i + 1] - 1;
            
        std::string newGeneString = geneNames.substr(start, end);
        genes.push_back(newGeneString);
    }
    
    tuDataFile.close();
}



void Operon::description(){
    
    std::cout << "Transcritional Unit: " << name << std::endl;
    std::cout << "\tGenes: ";
    int numberOfGenes = (int)genes.size();
    for (int i = 0; i < numberOfGenes; i++) {
        std::cout << genes[i] << "\t";
    }
    
}
    
    
    
}   //      namespace ustc