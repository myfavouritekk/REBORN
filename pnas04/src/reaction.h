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
 * (8) binding reaction between an inducer and a protein
 */

#include "node.h"
#include "consts.h"

#include <algorithm>
#include <vector>

namespace ustc{

enum reaction_type{
            TRANSCRIPTION,
            DEGRADATION,
            BINDING,
            MODIFICATION,
            PARTIAL_DEGRADATION,
            DIMERIZATION,
            CATALYTIC_DEGRADATION,
            PARTIAL_CATALYTIC_DEGRADATION,
            INDU_PROT_BINDING
            };
    
class Reaction{

    public:
        
        Reaction();
        Reaction(int _rtype);
        Reaction(Reaction &reaction);//copy constructor
        ~Reaction();

		void setReversible(const bool& rev);
		bool isReversible ();

        void modifyForwardRate();
        void modifyReverseRate();
		void initForwardRateRandomly();	//	0~1
		void initReverseRateRandomly();	//	0~1
		void setForwardRate(double rate);

		void addReactant( Node* sr);//no const
		void addModifier( Node* sm);//no const
		void addProduct ( Node* sp);//no const

		int  getRtype ();
		Node* getReactant(const int& ir);
		Node* getModifier(const int& im);
		Node* getProduct (const int& ip);
        
        std::vector<Node*>* getReactantsVector();
        std::vector<Node*>* getModifiersVector();
        std::vector<Node*>* getProductsVector();

        double getForwardRate();
        double getReverseRate();
    
        int getReactantsSize();
        int getModifiersSize();
        int getProductsSize();
    
        //adding method to judge whether a node is in this reaction
        bool containNode( Node* aNode);
        int containNodeAsReactant( Node* aNode);
        int containNodeAsModifier( Node* aNode);
        int containNodeAsProduct( Node* aNode);
        
        //	check equality with given reaction
		bool operator==(const Reaction& r1) const;
    
        friend double Node::ode(std::vector<Reaction*> reactionList, double *y, double t);

        //output method
        void description(int reactionIndex);
    

    private:

		//	reaction type
		reaction_type rtype;

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
} //namespace ustc

#endif
