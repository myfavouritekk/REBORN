//
//  plasmid.cpp
//  PNAS04Team
//
//  Created by Jack Kang on 12-8-1.
//
//

#include "plasmid.h"


namespace ustc {


Plasmid::Plasmid(){
    
}


DNAPiece* Plasmid::getDnaPiece(const int& dnaPieceIndex){
    int size = (int)dnaPieces.size();
    if (dnaPieceIndex >= 0 && dnaPieceIndex < size) {
        return dnaPieces[dnaPieceIndex];
    }
    return NULL;
}

void Plasmid::addDnaPiece(DNAPiece* newDNAPiece){
    dnaPieces.push_back(newDNAPiece);
}
    
    
void Plasmid::readMotifs(const int& cellIndex){

	std::stringstream ss;
	ss << "Cell_"<<cellIndex <<"_Motifs.txt";
	std::ifstream infile;
	infile.open(ss.str().c_str());
	if (!infile) {
        std::cerr << "Error: unable to open input file: " << infile << std::endl;
		exit(1);
    }

	infile >> numSingleMotifs >> numDoubleMotifs >> numTripleMotifs;
	if (infile.bad ()) throw std::runtime_error ("IO stream corrupted");
	if (infile.fail ()) throw std::runtime_error ("bad data");
	if (numSingleMotifs + numDoubleMotifs + numTripleMotifs == 0) {
		std::cerr << "Error: empty data" << std::endl;
		exit(1);
	}

	// creat singleMotifMatrice
	singleMotifsMatrice = new int**[numSingleMotifs];
	for(int i = 0; i < numSingleMotifs; i++){
		singleMotifsMatrice[i] = new int*[1];
		singleMotifsMatrice[i][0] = new int[1];
	}
	
	// creat doubleMotifsMatrice
	doubleMotifsMatrice = new int**[numDoubleMotifs];
	for(int i = 0; i < numDoubleMotifs; i++){
		doubleMotifsMatrice[i] = new int*[2];
		for(int j = 0; j < 2; j++){
			doubleMotifsMatrice[i][j]=new int[2];
		}
	}

	// creat tripleMotifsMatrice
	tripleMotifsMatrice = new int**[numTripleMotifs];
	for(int i = 0; i < numTripleMotifs; i++){
		tripleMotifsMatrice[i] = new int*[3];
		for(int j = 0; j < 3; j++){
			tripleMotifsMatrice[i][j] = new int[3];
		}
	}

	// read singleMotifMatrice
	for(int i = 0; i < numSingleMotifs; i++){
		infile >> singleMotifsMatrice[i][0][0];
	}

	// read doubleMotifMatrice
	for(int i = 0; i < numDoubleMotifs; i++){
		infile >> doubleMotifsMatrice[i][0][0] >> doubleMotifsMatrice[i][0][1];
		infile >> doubleMotifsMatrice[i][1][0] >> doubleMotifsMatrice[i][1][1];
	}
	
	// read tripleMotifMatrice
	for(int i = 0;i < numTripleMotifs; i++){
		infile >> tripleMotifsMatrice[i][0][0] >> tripleMotifsMatrice[i][0][1] >> tripleMotifsMatrice[i][0][2];
		infile >> tripleMotifsMatrice[i][1][0] >> tripleMotifsMatrice[i][1][1] >> tripleMotifsMatrice[i][1][2];
		infile >> tripleMotifsMatrice[i][2][0] >> tripleMotifsMatrice[i][2][1] >> tripleMotifsMatrice[i][2][2];
	}
}

}// namespace ustc