//
//  dnapiece.cpp
//  PNAS04Team
//
//  Created by Jack Kang on 12-8-1.
//
//

#include "dnapiece.h"


namespace ustc {

//  constructor
DNAPiece::DNAPiece(
                   Part* _promoter,
                   Part* _rbs,
                   Part* _gene,
                   Part* _terminator
                   ){
    
    
    promoter = _promoter;
    rbs = _rbs;
    gene = _gene;
    terminator = _terminator;
    
    direction = CLOCKWISE;  // default direction
}
    
Part* DNAPiece::getPromoter(){
    return promoter;
}

Part* DNAPiece::getRbs(){
    return rbs;
}
    
Part* DNAPiece::getGene(){
    return gene;
}
    
Part* DNAPiece::getTerminator(){
    return terminator;
}
    
void DNAPiece::setDirection(const piece_direction& _direction){
    direction = _direction;
}
    
piece_direction DNAPiece::getDirection(){
    return direction;
}
    
}// namespace ustc