//
//  dnapiece.h
//  PNAS04Team
//
//  Created by Jack Kang on 12-8-1.
//
//

#ifndef __PNAS04Team__dnapiece__
#define __PNAS04Team__dnapiece__

#include <iostream>
#include "part.h"
#include "consts.h"

namespace ustc {

enum piece_direction{CLOCKWISE, ANTI_CLOCKWISE};
    
class DNAPiece{
    
private:
    
    Part* promoter;
    Part* rbs;
    Part* gene;
    Part* terminator;
    
    piece_direction direction;
    
public:
    
    //  constructors
    //
    //  constructors for top-down design
    //
    DNAPiece();     //      default constructor
    DNAPiece(
             Part* _promoter,
             Part* _rbs,
             Part* _gene,
             Part* _terminator
             );
    DNAPiece(
             Part* _promoter,
             Part* _rbs,
             Part* _gene,
             Part* _terminator,
             piece_direction _direction
             );
    
    //  constructors for bottom-up design
    DNAPiece(
             std::string promoterId,
             std::string rbsId,
             std::string geneId,
             std::string terminatorId
             );
    DNAPiece(
             std::string promoterId,
             std::string rbsId,
             std::string geneId,
             std::string terminatorId,
             piece_direction _direction
             );
    
    //  destructor
    ~DNAPiece();
    
    
    Part* getPromoter();
    Part* getRbs();
    Part* getGene();
    Part* getTerminator();
   
    void setDirection(
                const piece_direction& _direction
                );
    
    piece_direction getDirection();
    
};  //  class DNAPiece
    
    
    
    
}// namespace ustc



#endif /* defined(__PNAS04Team__dnapiece__) */
