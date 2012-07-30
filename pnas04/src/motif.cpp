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
    int size = motifNodes.size();
	switch (size){ //single gene regulation
		case 1:{
			if (motifMatrix[0][0] == 1){
				mtype = 1; // gene promotes itself
			}
			else{
				mtype = 2; // gene represses itself
			}
			break;
		}
		case 2:{ // two genes regulation
			int testValue = (motifMatrix[0][1] + motifMatrix[1][0]);
			switch (testValue){
				case 1:{
					/* 0 1  or  0 0
					   0 0      1 0
					*/
					mtype = 3;
					break;
				}
				case -1:{
					/* 0 -1  or  0 0
					   0  0     -1 0
					*/
					mtype = 4;
					break;
				}
				case 2:{
					/* 0 1   
					   1 0      
					*/
					mtype = 5;
					break;
				}
				case -2:{
					/* 0 -1  
					  -1  0      
					*/
					mtype = 6;
					break;
				}
				case 0:{
					/* 0 1  or  0 -1
					  -1 0      1  0
					*/
					mtype = 7;
					break;
				}
				default:{
					std::cout << "wrong input" << std::endl;
					break;
				}
			}
			break;
		}
		case 3:{
			mtype = 8;
			break;
		}
		default:{
			std:: cout << "wrong input" << std::endl;
			break;
		}
		
	}
	    
}

//print the basic information of this motif
void Motif::description(){
    int size = (int)motifNodes.size();
    std::cout << "Motif:" << std::endl << "\tType: " << mtype << std::endl;
    std::cout << "\tNodes: ";
    for (int i = 0; i < size; i++) {
        std::cout << motifNodes[i] -> getNstring() << "\t";
    }
    std::cout << std::endl << "\tMatrix:" << std::endl << "\t\t";
    for (int i = 0; i < size; i++) {
        std::cout << motifNodes[i] -> getNstring() << "\t";
    }
    for (int i = 0; i < size; i++) {
        std::cout << std::endl << "\t" << motifNodes[i] -> getNstring() << "\t";
        for (int j = 0; j < size; j++) {
            std::cout << motifMatrix[i][j] << "\t";
        }
    }
    std::cout << std::endl;
}

    
    
}   //namespace ustc