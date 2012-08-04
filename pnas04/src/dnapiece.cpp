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
                   std::string promoterId,
                   std::string rbsId,
                   std::string geneId,
                   std::string terminatorId
                   ){
    promoter = new Part(PROMOTER, promoterId);
    rbs = new Part(RBS, rbsId);
    gene = new Part(GENE_PART, geneId);
    terminator = new Part(TERMINATOR, terminatorId);
    
    direction = CLOCKWISE;
}
DNAPiece::DNAPiece(
                   std::string promoterId,
                   std::string rbsId,
                   std::string geneId,
                   std::string terminatorId,
                   piece_direction _direction
                   ){
    promoter = new Part(PROMOTER, promoterId);
    rbs = new Part(RBS, rbsId);
    gene = new Part(GENE_PART, geneId);
    terminator = new Part(TERMINATOR, terminatorId);
    
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