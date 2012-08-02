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
#include "dnapiece.h"
#include "consts.h"

namespace ustc{

    
class Plasmid{
    
private:
    
    std::vector<DNAPiece*> dnaPieces;
    
    int numSingleMotifs;
    int numDoubleMotifs;
    int numTripleMotifs;
    
    int*** singleMotifsMatrice;
    int*** doubleMotifsMatrice;
    int*** tripleMotifsMatrice;
    
    
public:
    
    //  constructor
    Plasmid();
    
    DNAPiece* getDnaPiece(const int& dnaPieceIndex);
    
    void addDnaPiece(DNAPiece* newDnaPiece);
    
    void readMotifs(const int& cellIndex);
};



}// namespace ustc



#endif /* defined(__PNAS04Team__plasmid__) */
