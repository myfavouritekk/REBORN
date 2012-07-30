//
//  motif.cpp
//  PNAS04Team
//
//  Created by Jack Kang on 12-7-30.
//
//

#include "motif.h"

namespace ustc {
    
Motif::Motif(std::vector<Node*>* constructingNodes, int** relationshipMatrix){
	vector<Node*>::iterator iter = (*constructingNodes).begin();
	vector<Node*>::iterator iter_end = (*constructingNodes).end();
	int matrixSize =(*constructingNodes).size();
	while (iter!=iter_end){
		motifNodes.push_back(*iter);
		iter++;
	}
	motifMatrix = new int*[matrixSize];
	for(int i = 0; i<matrixSize;i++){
		motifMatrix[i] = new int[matrixSize];
	}
	for(int i = 0; i<matrixSize; i++){
		for(int j=0; j<matrixSize; j++){
			motifMatrix[i][j]=relationshipMatrix[i][j];
		}
	}
	generateMtype();
}

Motif::~Motif(){
    
}


int Motif::getMtype(){
    return mtype;
}

//return motif size
int Motif::motifSize(){
    return (int)motifNodes.size();
}


//private method to generate motif type based on the motifMatrix
void Motif::generateMtype(){
    
}

    
    
}   //namespace ustc