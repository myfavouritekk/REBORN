//
//  dnapiece.cpp
//  PNAS04Team
//
//  Created by Jack Kang on 12-8-1.
//
//

#include "dnapiece.h"


namespace ustc {

//  constructors
    
DNAPiece::DNAPiece(){
    promoter = NULL;
    rbs = NULL;
    gene = NULL;
    terminator = NULL;
    
    direction = CLOCKWISE;
}
    
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
    
DNAPiece::DNAPiece(
                   Part* _promoter,
                   Part* _rbs,
                   Part* _gene,
                   Part* _terminator,
                   piece_direction _direction
                   ){
    
    
    promoter = _promoter;
    rbs = _rbs;
    gene = _gene;
    terminator = _terminator;
    
    direction = _direction;
}
    
    
DNAPiece::DNAPiece(
                   std::string promoterName,
                   std::string rbsName,
                   std::string geneName,
                   std::string terminatorName
                   ){
    promoter = new Part(PROMOTER, promoterName);
    rbs = new Part(RBS, rbsName);
    gene = new Part(GENE_PART, geneName);
    terminator = new Part(TERMINATOR, terminatorName);
    
    direction = CLOCKWISE;
}
DNAPiece::DNAPiece(
                   std::string promoterName,
                   std::string rbsName,
                   std::string geneName,
                   std::string terminatorName,
                   piece_direction _direction
                   ){
    promoter = new Part(PROMOTER, promoterName);
    rbs = new Part(RBS, rbsName);
    gene = new Part(GENE_PART, geneName);
    terminator = new Part(TERMINATOR, terminatorName);
    
    direction = _direction;

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