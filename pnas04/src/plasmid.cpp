//
//  plasmid.cpp
//  PNAS04Team
//
//  Created by Jack Kang on 12-8-1.
//
//

#include "plasmid.h"


namespace ustc {


Plasmid::Plasmid(){
    
}


DNAPiece* Plasmid::getDnaPiece(const int& dnaPieceIndex){
    int size = (int)dnaPieces.size();
    if (dnaPieceIndex >= 0 && dnaPieceIndex < size) {
        return dnaPieces[dnaPieceIndex];
    }
    return NULL;
}

void Plasmid::addDnaPiece(DNAPiece* newDNAPiece){
    dnaPieces.push_back(newDNAPiece);
}
    
    
}// namespace ustc