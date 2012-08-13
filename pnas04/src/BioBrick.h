//
//  biobrick.h
//  PNAS04Team
//
//  Created by Jack Kang on 12-8-11.
//
//

#ifndef __PNAS04Team__biobrick__
#define __PNAS04Team__biobrick__

#include <iostream>
#include "consts.h"


namespace ustc {


class BioBrick{
  
public:
    
    //      constructor and destructor
    BioBrick(const std::string _geneName, const std::string _promoterName);
    ~BioBrick();
    
    //      description method
    std::string description();
    
    //      whether promoter and gene is available in the database
    bool isAvailableInDatabase();
    
    
    //      getters
    std::string getPromoterName();
    std::string getGeneName();
    
    
    
private:
    
    std::string itsPromoterName;
    std::string itsGeneName;
    
    bool hasData;
    
    
};

    
    
    
    
}   //      namespace ustc



#endif /* defined(__PNAS04Team__biobrick__) */
