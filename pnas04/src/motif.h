//
//  motif.h
//  PNAS04Team
//
//  Created by Jack Kang on 12-7-30.
//
//

#ifndef __PNAS04Team__motif__
#define __PNAS04Team__motif__

#include <iostream>
#include "consts.h"
#include "node.h"

namespace ustc {
    
class Motif{
    public:
        //constructor: no default constructor
        Motif(std::vector<Node*>* constructingNodes, int** relationshipMatrix);
        ~Motif();
        
        //return the motif type
        int getMtype();
        
        //return motif size
        int motifSize();
        
    private:
        //private method to generate motif type based on the motifMatrix
        void generateMtype();
    
        //containing the pointers of the constructing nodes of the motif
        std::vector<Node*> motifNodes;
    
        //storing the regulatory relationships in the motif
        int** motifMatrix;
    
        //motif type
        int mtype;
    
};  //class Motif definition
    
    
}   // namespace ustc


#endif /* defined(__PNAS04Team__motif__) */
