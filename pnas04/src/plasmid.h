//
//  plasmid.h
//  PNAS04Team
//
//  Created by Jack Kang on 12-8-1.
//
//

#ifndef __PNAS04Team__plasmid__
#define __PNAS04Team__plasmid__

#include <iostream>
#include <vector>
#include <sstream>
#include <fstream>
#include <stdexcept>
#include "operon.h"
#include "consts.h"
#include "biobrick.h"

namespace ustc{
    
    
    //==========================================================//
    //      structs used to find whole structure                //
    //==========================================================//
    
    
struct Regulation{
    std::string itsRegulateeNames;
    int regulateDirection;
};
    
struct RegulatorCandidates{
    std::string name;
    std::vector<Regulation* > itsRegulatees;
};


    
struct CompleteCandidate{
    int numOfRegulators;
    int numOfRegulatees;
    std::vector<RegulatorCandidates*> regulators;
    std::vector<std::string> regulateeStrings;
};
 




//******************************************************************//
//*     class Plasmid: constructed by real parts                   *//
//******************************************************************//
    
    
class Plasmid{
    
    
public:
    
    //      constructor
    Plasmid();
    

    
    //==========================================================//
    //      PART 1: methods related to models of the sofrware   //
    //==========================================================//
    
    
    //      read complete regulatory matrix from a text file
    void readCompleteMatrix(const int cellIndex);
    
    //      find suitable strctures to rebuild itself based on the whole matrix of operon-operon relationships
    void findCompleteCondidates(
                                const int numRow,
                                const int numColumn,
                                const std::string* namesOfRegualters,
                                const std::string* namesOfRegulatees,
                                const int** database
                                );

    //      find suitable strctures to rebuild itself based on the whole matrix of promoter-gene relationships
    void findCompleteCondidatesUsingBiobricks(
                                const int numRow,
                                const int numColumn,
                                const std::string* namesOfRegualters,
                                const std::string* namesOfRegulatees,
                                const int** database
                                );

    
    //==========================================================//
    //      PART 2: generate build plans based on the cadidates //
    //==========================================================//
    
    void generatePlans();
    void generatePlanOutputs(const int plasmidIndex);
    
private:
    
    //==========================================================//
    //      members related to real parts                       //
    //==========================================================//
    std::vector<std::vector<Operon*> > plans;
    std::vector<std::vector<BioBrick*> > biobrickPlans;
    
    
    //==========================================================//
    //      members used to store whole matrix                  //
    //==========================================================//
    int** wholeRegulatoryMatrix;
    int numOfGenes;
    std::vector<CompleteCandidate*> completeCandidates;
    
    
};



}// namespace ustc



#endif /* defined(__PNAS04Team__plasmid__) */
