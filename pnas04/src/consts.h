//
//  consts.h
//  PNAS04Team
//
//  Created by Jack Kang on 12-7-19.
//
//

#ifndef PNAS04Team_consts_h
#define PNAS04Team_consts_h

//Total evolution number of the program
#define TOTAL_EVO 750

//Total number of cells
#define POPULATION 1000


//parameters for mutation in Cell
#define PROB_MUT_DEG_PROT 0.5
#define PROB_MUT_KIN_CONST 1.0
#define PROB_MUT_ADD_GENE 0.4
#define PROB4 0.1
#define PROB_MUT_ADD_POSTMOD 0.1


//parameter multiplied to rate when modify forward or reverse rates in Reaction
#define RATE_MULTI 0.7

//parameters that add nodes size, components complexity and rlist size in the score function
#define PARAMETER_NODE_SIZE 0.2
#define PARAMETER_COMPLEX_SIZE 0.05
#define PARAMETER_REACTION_SIZE 0.1



#endif
