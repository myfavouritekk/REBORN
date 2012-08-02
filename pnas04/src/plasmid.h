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
#include "consts.h"

namespace ustc{
    
    
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
    
    
class Plasmid{
    
private:
    
    std::vector<DNAPiece*> dnaPieces;
    
    int numSingleMotifs;
    int numDoubleMotifs;
    int numTripleMotifs;
    int numOfGenes;
    
    int*** singleMotifsMatrice;
    int*** doubleMotifsMatrice;
    int*** tripleMotifsMatrice;
    int** wholeMotifMatrix;
    
    std::vector<MotifCandidate*> candidates;
    
public:
    
    //  constructor
    Plasmid();
    
    DNAPiece* getDnaPiece(const int& dnaPieceIndex);
    
    void addDnaPiece(DNAPiece* newDnaPiece);
    
    void readMotifs(const int& cellIndex);
    
    void findCandidates(
                            const int& numRow,
                            const int& numColumn,
                            const std::string* namesOfGenes,
                            const std::string* namesOfPromoters,
                            const int** database
                            );
};



}// namespace ustc



#endif /* defined(__PNAS04Team__plasmid__) */
