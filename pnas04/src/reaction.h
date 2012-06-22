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

#include <ctime>
#include <algorithm>
#include <vector>

class Reaction{

    public:
        
        Reaction();
        Reaction(int _rtype);
        ~Reaction();

		void setReversible(const bool& rev);
		bool isReversible ();

        void modifyForwardRate();
        void modifyReverseRate();
		void initForwardRateRandomly();	//	0~1
		void initReverseRateRandomly();	//	0~1

		void addReactant( Node* sr);//no const
		void addModifier( Node* sm);//no const
		void addProduct ( Node* sp);//no const

		int  getRtype ();
		Node* getReactant(const int& ir);
		Node* getModifier(const int& im);
		Node* getProduct (const int& ip);

		//	check equality with given reaction
		bool operator==(const Reaction& r1) const;
    
        friend double Node::ode(std::vector<Reaction*> reactionList, double *y, double t);
    
    
        //get a copy of this reaction for growth phase
        Reaction* aNewCopy();

    private:

		//	reaction type
		int rtype;

        //  reversibility
        bool reversible;

        //  reaction rates
        double forwardRate;
        double reverseRate;

        //  species involved
        std::vector<Node*> reactants;
        std::vector<Node*> modifiers;
        std::vector<Node*> products;
};


#endif
