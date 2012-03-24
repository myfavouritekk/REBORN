#ifndef REACTION_H
#define REACTION_H

/*
 * Representations
 * (0) Transcription of a gene from a free 
 *     promoter or a bound promoter/protein complex
 * (1) degradation of protein
 * (2) binding reaction between a protein and 
 *     a gene or an existing gene/protein complex
 * (3) posttranscriptional modification
 * (4) partial degradation, with one protein of the 
 *     complex chosen as the degradation product
 * (5) dimerization
 * (6) catalytic degradation
 * (7) partial catalytic degradation
 */

#include "node.h"

class Reaction{

    public:
        
        Reaction();
        ~Reaction();

        void modifyForwardRate();
        void modifyReverseRate();
        :A

:A

    private:

        //  reversibility
        bool isReversible;

        //  reaction rates
        double forwardRate;
        double reverseRate;

        //  species involved
        std::vector<Node*> reactants;
        std::vector<Node*> modifiers;
        std::vector<Node*> products;
};




#endif
