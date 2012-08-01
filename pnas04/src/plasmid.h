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
    
public:
    
    //  constructor
    Plasmid();
    
    DNAPiece* getDnaPiece(const int& dnaPieceIndex);
    
    void addDnaPiece(DNAPiece* newDnaPiece);
    
};



}// namespace ustc



#endif /* defined(__PNAS04Team__plasmid__) */
