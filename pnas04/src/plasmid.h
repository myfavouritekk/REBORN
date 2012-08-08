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
#include "dnapiece.h"
#include "operon.h"
#include "consts.h"

namespace ustc{
    
    //==========================================================//
    //      enum and structs used to find motifs                //
    //==========================================================//
    
enum motifSize{ _1X1 = 1, _2X2, _3X3};

struct GeneCandidate{
    std::string name;
    std::vector<std::string> itsPromoterStrings;
};


struct MotifCandidate{
    motifSize size;
    std::vector<std::string> geneString;
    std::vector<GeneCandidate*> genes;
    std::vector<std::string> promoterStrings;
};
    
    
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
    //      PART 1: methods related to real parts               //
    //==========================================================//
    
    DNAPiece* getDnaPiece(const int dnaPieceIndex);
    
    void addDnaPiece(DNAPiece* newDnaPiece);
    
    
    //==========================================================//
    //      PART 2: methods related to models of the sofrware   //
    //==========================================================//
    
    //      read motifs from a text file and store them
    void readMotifs(const int cellIndex);
    
    //      find suitable structures to rebuild itself based on regulatory matrice of motifs
    void findMotifsCandidates(
                                const int numRow,
                                const int numColumn,
                                const std::string* namesOfGenes,
                                const std::string* namesOfPromoters,
                                const int** database
                                );
    
    //      read complete regulatory matrix from a text file
    void readCompleteMatrix(const int cellIndex);
    
    //      find suitable strctures to rebuild itself based on the whole matrix
    void findCompleteCondidates(
                                const int numRow,
                                const int numColumn,
                                const std::string* namesOfRegualters,
                                const std::string* namesOfRegulatees,
                                const int** database
                                );
    
    //==========================================================//
    //      PART 3: generate build plans based on the cadidates //
    //==========================================================//
    
    void generatePlans();
    void generatePlanOutputs(const int plasmidIndex);
    
private:
    
    //==========================================================//
    //      members related to real parts                       //
    //==========================================================//
    std::vector<DNAPiece*> dnaPieces;
    std::vector<std::vector<Operon*> > plans;
    
    //==========================================================//
    //      members used to store motifs generated by software  //
    //==========================================================//
    int numSingleMotifs;
    int numDoubleMotifs;
    int numTripleMotifs;
    
    int*** singleMotifsMatrice;
    int*** doubleMotifsMatrice;
    int*** tripleMotifsMatrice;
    
    std::vector<MotifCandidate*> motifCandidates;
    
    //==========================================================//
    //      members used to store whole matrix                  //
    //==========================================================//
    int** wholeRegulatoryMatrix;
    int numOfGenes;
    std::vector<CompleteCandidate*> completeCandidates;
    
    
};



}// namespace ustc



#endif /* defined(__PNAS04Team__plasmid__) */
