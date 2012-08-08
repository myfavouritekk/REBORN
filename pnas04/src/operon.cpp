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
    
    
Operon::Operon(const std::string _name){
    
    
    //      using the identifier as the name of a file and query other information
    
    std::ifstream tuDataFile;
    std::stringstream tuDataFileName;
    tuDataFileName << OPERONS_PATH << _name << ".txt";
    tuDataFile.open(tuDataFileName.str().c_str());
    
    //      initialize the members
    std::string geneNames;
    int numGene;
    tuDataFile >> name >> numGene >> geneNames;
    int* nameIndex = new int[numGene + 1];
    nameIndex[0] = 0;
    for (int i = 0; i < numGene; i++) {
        nameIndex[i + 1] = (int)geneNames.find(",", nameIndex[i] + 1);
        if (i == numGene - 1) {//   there is no last ",", nameIndex[i+1] is -1
            nameIndex[i + 1] = (int)geneNames.size();   //  assign it to the end of the string
        }
        int start, length;
        if (i == 0){
            start = 0;
            length = nameIndex[i + 1] - nameIndex[i];
        }
        else{
            start = nameIndex[i] + 1;
            length = nameIndex[i + 1] - nameIndex[i] - 1;
        }
        if (i == numGene - 1)
            length = (int)geneNames.size(); //  certainly out of range and will read to the end of the string
            
        std::string newGeneString = geneNames.substr(start, length);
        genes.push_back(newGeneString);
    }
    
    tuDataFile.close();
}



std::string Operon::description(){
    
    std::stringstream result;
    result << "Transcritional Unit: " << name << std::endl;
    result << "\tGenes: ";
    int numberOfGenes = (int)genes.size();
    for (int i = 0; i < numberOfGenes; i++) {
        result << genes[i] << "\t";
    }
    result << std::endl;
    return result.str();
}
    
    
    
}   //      namespace ustc