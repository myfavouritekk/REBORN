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


DNAPiece* Plasmid::getDnaPiece(const int dnaPieceIndex){
    int size = (int)dnaPieces.size();
    if (dnaPieceIndex >= 0 && dnaPieceIndex < size) {
        return dnaPieces[dnaPieceIndex];
    }
    return NULL;
}

void Plasmid::addDnaPiece(DNAPiece* newDNAPiece){
    dnaPieces.push_back(newDNAPiece);
}
    
    
void Plasmid::readMotifs(const int cellIndex){

	std::stringstream ss;
	ss << SAVES_PATH << "Cell_"<<cellIndex <<"_Motifs.txt";
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
	infile.close();
}

    
void Plasmid::findMotifsCandidates(
                             const int numRow,
                             const int numColumn,
                             const std::string* namesOfGenes,
                             const std::string* namesOfPromoters,
                             const int** database
                             ){

	//find single motif cadidates
	// j -> gene
	// k -> promoter
    for(int i = 0; i < numSingleMotifs; i ++){
		for(int j = 0; j < numColumn; j ++){
			for(int k = 0; k < numRow; k ++){
				if(database[k][j] == singleMotifsMatrice[i][0][0]){
					MotifCandidate* mc;
					GeneCandidate* gc;
					gc -> name = namesOfGenes[j];
					gc -> itsPromoterStrings.push_back(namesOfPromoters[k]);
					mc -> size = _1X1;
					mc -> geneString.push_back(namesOfGenes[j]);
					mc -> genes.push_back(gc);
					mc -> promoterStrings.push_back(namesOfPromoters[k]);
					motifCandidates.push_back(mc);
				}
			}
		}
	}



	//find double motif cadidats
	//i -> gene1
	//k -> gene2
	//j -> promt1
	//l -> promt2
	for(int n = 0; n < numDoubleMotifs; n ++){
		for(int i = 0; i < numColumn; i ++){
			for(int j = 0; j < numRow; j ++){
				// promt1 satisfy motifmatrice[n][0][0]
				if(database[j][i] != doubleMotifsMatrice[n][0][0]){ 
					continue;
				}
				for(int k = 0; k < numColumn; k ++){ 
					//gene2 != gene1  ||   gene2 satisfy motifmatrice[n][0][1]
					if(k == i || database[j][k] != doubleMotifsMatrice[n][0][1]){
						continue;
					}
					for(int l = 0; l < numRow; l++){
						// this will satisfy all condition
						if(database[l][i] == doubleMotifsMatrice[n][1][0] && database[k][l] == doubleMotifsMatrice[n][1][1]){
							GeneCandidate* gc1;
							GeneCandidate* gc2;
							MotifCandidate* mc;
							gc1 -> name = namesOfGenes[i];
							gc1 -> itsPromoterStrings.push_back(namesOfPromoters[j]);
							gc2 -> name = namesOfGenes[k];
							gc2 -> itsPromoterStrings.push_back(namesOfPromoters[l]);
							mc -> size = _2X2;
							mc -> geneString.push_back(gc1 -> name);
							mc -> geneString.push_back(gc2 -> name);
							mc -> promoterStrings.push_back(namesOfPromoters[j]);
							mc -> promoterStrings.push_back(namesOfPromoters[l]);
							mc -> genes.push_back(gc1);
							mc -> genes.push_back(gc2);
							motifCandidates.push_back(mc);
						}
					}
				}
			}
		}
	}
	


	// find trible motif candidates
	// i -> gene1
	// k -> gene2
	// p -> gene3
	// j -> promt1
	// l -> promt2
	// q -> promt3
	for(int n = 0; n < numTripleMotifs; n ++){
		for(int i = 0; i < numColumn; i ++){
			for(int j = 0; j < numRow; j ++){
				// promt1 satisfy motifmatrice[n][0][0]
				if(database[j][i] != tripleMotifsMatrice[n][0][0]){ 
					continue;
				}
				for(int k = 0; k < numColumn; k ++){ 
					//gene2 != gene1  ||   gene2 satisfy motifmatrice[n][0][1]
					if(k == i || database[j][k] != tripleMotifsMatrice[n][0][1]){
						continue;
					}
					for(int l = 0; l < numRow; l++){
						if(database[l][i] != tripleMotifsMatrice[n][1][0] || database[l][k] != tripleMotifsMatrice[n][1][1]){
							continue;
						}
						for(int p = 0; p < numColumn; p ++){
							if(p ==	k || p == i || database[j][p] != tripleMotifsMatrice[n][0][2] || database[l][p] != tripleMotifsMatrice[n][1][2]){
								continue;
							}
							for(int q  = 0; q < numRow; q ++){
								//satisfy all condition
								if(database[q][i] == tripleMotifsMatrice[n][2][0] && database[q][k] == tripleMotifsMatrice[n][2][1] && database[q][p] == tripleMotifsMatrice[n][2][2]){
									GeneCandidate* gc1;
									GeneCandidate* gc2;
									GeneCandidate* gc3;
									MotifCandidate* mc;
									gc1 -> name = namesOfGenes[i];
									gc1 -> itsPromoterStrings.push_back(namesOfPromoters[j]);
									gc2 -> name = namesOfGenes[k];
									gc2 -> itsPromoterStrings.push_back(namesOfPromoters[l]);
									gc3 -> name = namesOfGenes[p];
									gc3 -> itsPromoterStrings.push_back(namesOfPromoters[q]);
									mc -> size = _3X3;
									mc -> geneString.push_back(gc1 -> name);
									mc -> geneString.push_back(gc2 -> name);
									mc -> geneString.push_back(gc3 -> name);
									mc -> promoterStrings.push_back(namesOfPromoters[j]);
									mc -> promoterStrings.push_back(namesOfPromoters[l]);
									mc -> promoterStrings.push_back(namesOfPromoters[q]);
									mc -> genes.push_back(gc1);
									mc -> genes.push_back(gc2);
									mc -> genes.push_back(gc3);
									motifCandidates.push_back(mc);
								}
							}
						}
					}
				}
			}
		}
	}

}

    
void Plasmid::readCompleteMatrix(const int cellIndex){
    std::stringstream fileName;
    fileName << SAVES_PATH << "Cell_" << cellIndex << "_Complete.txt";
    std::fstream file;
    file.open(fileName.str().c_str());
    file >> numOfGenes;
    wholeRegulatoryMatrix = new int*[numOfGenes];
    for (int i = 0; i < numOfGenes; i++) {
        wholeRegulatoryMatrix[i] = new int[numOfGenes];
        for (int j = 0; j < numOfGenes; j++) {
            file >> wholeRegulatoryMatrix[i][j];
        }
    }
    file.close();
}


    
    
    //==============================================================================//
    //                                                                              //
    //  findMatrixRecurssion method                                                 //
    //  purpose: find small matrix in a big matrix                                  //
    //  it uses recurssion algorithm to find MxM matrix in NxN matrix, based on     //
    //      recussion of finding (M-1)x(M-1) matrix in (N-1)x(N-1) matrix           //
    //  parameters:                                                                 //
    //          databaseMatrix is the global giant matrix                           //
    //          targetMatrix is the current matrix to be find                       //
    //          choicesPool contains the indice of rows or coloumns can be chosen   //
    //          numberOfChoicesInPool is the size of choicesPool                    //
    //          numberOfChoicesToBeChosen is the size of targetMatrix               //
    //          numberOfPossibleChoicesSets is the number of possible solutions     //
    //  return value:                                                               //
    //          return value is the 2-d array that contains the possible solutions  //
    //                                                                              //
    //==============================================================================//
int** findMatrixRecurssion(
                           const int** databaseMatrix,
                           const int** targetMatrix,
                           const int* choicesPool,
                           const int numberOfChoicesInPool,
                           const int numberOfChoicesToBeChosen,
                           int* numberOfPossibleChoiceSets
                           )
{
    int** findMatrixRecurssion(
                               const int** databaseMatrix,
                               const int** targetMatrix,
                               const int* choicesPool,
                               const int numberOfChoicesInPool,
                               const int numberOfChoicesToBeChosen,
                               int* numberOfPossibleChoiceSets
                               );
   
    //
    //      for choosing 1 element from the choicePool,
    //      which is also the end of recurssion.
    //
    if(numberOfChoicesToBeChosen == 1){
        int numPossible = 0;
        std::vector<int*> choicesVector;
        for (int i = 0; i < numberOfChoicesInPool; i++ ) {
            if(databaseMatrix[choicesPool[i]][choicesPool[i]] == 2 ||
               databaseMatrix[choicesPool[i]][choicesPool[i]] == targetMatrix[0][0]){ //  2 means it can not only up-regulate, but also down-regulate
                numPossible++;
                int* aChoice = new int(choicesPool[i]);
                choicesVector.push_back(aChoice);
            }
        }
        int** choices = new int*[numPossible];
        for(int i = 0; i < numPossible; i++){
            choices[i] = choicesVector[i];
        }
        *numberOfPossibleChoiceSets = numPossible;
        return choices;
    }
    
    
    //
    //      general situation
    //
    
    //      variables for recurssion
    int** newTargetMatrix = new int*[numberOfChoicesToBeChosen - 1];
    int numOfWorkingChoices = 0;
    std::vector<int*> choicesVector;
    
    
    //
    //      constructing new targetMatrix, which is the same as the old without first row and first coloumn
    //
    for(int i = 0; i < numberOfChoicesToBeChosen - 1; i++){
        newTargetMatrix[i] = new int[numberOfChoicesToBeChosen - 1];
        for(int j = 0; j < numberOfChoicesToBeChosen - 1; j++){
            newTargetMatrix[i][j] = targetMatrix[i + 1][j + 1];
        }
    }
    
    
    //      choose which one is suitable for the first element
    for(int i = 0; i < numberOfChoicesInPool; i++){
        int firstElementIndex = choicesPool[i];
        if (databaseMatrix[firstElementIndex][firstElementIndex] != 2 &&
            databaseMatrix[firstElementIndex][firstElementIndex] != targetMatrix[0][0]) { //  2 means it can not only up-regulate, but also down-regulate
            continue;
        }
        
        //
        //  constructing new choicePool for recurssion
        //
        int* newPool = new int[numberOfChoicesInPool - 1];
        for(int j = 0, k = 0; j < numberOfChoicesInPool; j++){
            if( i != j){
                newPool[k] = choicesPool[j];
                k++;
            }
        }
        
        int numberOfNewSets;
        int** possibleSolutions = findMatrixRecurssion(
                                                       databaseMatrix,
                                                       (const int**)newTargetMatrix,
                                                       newPool,
                                                       numberOfChoicesInPool - 1,
                                                       numberOfChoicesToBeChosen - 1,
                                                       &numberOfNewSets
                                                       );
        //
        //  judge whether the solutions work
        //
        for (int i = 0 ; i < numberOfNewSets; i++) {
            bool work = true;
            for (int j = 0; j < numberOfChoicesToBeChosen - 1; j++) {
                if ((databaseMatrix[possibleSolutions[i][j]][firstElementIndex] != 2 &&
                    targetMatrix[j + 1][0] != databaseMatrix[possibleSolutions[i][j]][firstElementIndex]) ||
                    (databaseMatrix[firstElementIndex][possibleSolutions[i][j]] != 2 &&
                    targetMatrix[0][j + 1] != databaseMatrix[firstElementIndex][possibleSolutions[i][j]])) {
                    work = false;
                    break;
                }
            }
            if (work) {
                numOfWorkingChoices++;
                int* aChoice = new int[numberOfChoicesToBeChosen];
                aChoice[0] = firstElementIndex;
                for (int k = 1; k < numberOfChoicesToBeChosen; k++) {
                    aChoice[k] = possibleSolutions[i][k - 1];
                }
                choicesVector.push_back(aChoice);
            }
        }
    }
    *numberOfPossibleChoiceSets = numOfWorkingChoices;
    int **choices = new int*[numOfWorkingChoices];
    for (int i = 0; i < numOfWorkingChoices; i++) {
        choices[i] = choicesVector[i];
    }
    return choices;
}

    
    
    
    
    
void Plasmid::findCompleteCondidates(
                                     const int numRowOfDatabase,
                                     const int numColumnOfDatabase,
                                     const std::string* namesOfRegualters,
                                     const std::string* namesOfRegulatees,
                                     const int** database
                                     )
{

    //  find wholeRegulatoryMatrix in the database,
    //  and store the suitable choices in the cadidateIndice array,
    //  each 1-d array stores the choices, and number of 1-ds is the numberOfCandidates
    int** candidateIndice;
    int numberOfCandidates;
    int* choicePool = new int[numColumnOfDatabase];
    for (int i = 0; i < numColumnOfDatabase; i++) {
        choicePool[i] = i;
    }
    candidateIndice = findMatrixRecurssion(database, (const int**)wholeRegulatoryMatrix, choicePool, numColumnOfDatabase, numOfGenes, &numberOfCandidates);
    
    //  construct condidates based on the results from findMatrixRecussion methods
    for (int i = 0; i < numberOfCandidates; i++) {
        CompleteCandidate* aNewCandidate = new CompleteCandidate;
        aNewCandidate -> numOfRegulatees = numOfGenes;
        aNewCandidate -> numOfRegulators = numOfGenes;
        
        
        //  constructing its' regulators and regulateeStrings
        for (int j = 0; j < numOfGenes; j++) {
            
            //  push back a regulatee's string
            (aNewCandidate -> regulateeStrings).push_back(namesOfRegulatees[candidateIndice[i][j]]);
            
            //  constructing a regulator
            RegulatorCandidates* aRegulator = new RegulatorCandidates;
            aRegulator -> name = namesOfRegualters[candidateIndice[i][j]];
            //  find its regulatees
            for (int k = 0; k < numOfGenes; k++) {
                //  find regulation that is not 0
                if (wholeRegulatoryMatrix[k][j] != 0) {
                    Regulation* aRegulation = new Regulation;
                    aRegulation -> itsRegulateeNames = namesOfRegulatees[candidateIndice[i][k]];
                    aRegulation -> regulateDirection = database[candidateIndice[i][k]][candidateIndice[i][j]];
                    aRegulator -> itsRegulatees.push_back(aRegulation);
                }
            }
            
            aNewCandidate -> regulators.push_back(aRegulator);
        }
        
        
        //  add this candidate
        completeCandidates.push_back(aNewCandidate);
    }   //  finished canstructing the completeCandidates vector
    
}

    
    
void Plasmid::generatePlans(){
    
    //          clear previous plans
    plans.clear();
    
    //          determine the number of plans
    int numberOfPlans = (int)completeCandidates.size();
    
    //          constructing every plan
    for (int i = 0; i < numberOfPlans; i++) {
        CompleteCandidate* aCandidate = completeCandidates[i];        
        std::vector<Operon*> aPlan;
        
        //      constructing every DNA piece of the plan
        for (int j = 0; j < numOfGenes; j++) {
            std::string operonName = (aCandidate -> regulators)[j] -> name;
            
            //  for case all regulatees are operons, and no promoters
            Operon* anOperon = new Operon(operonName);
            aPlan.push_back(anOperon);
        }
        
        plans.push_back(aPlan);
    }
    
}


void Plasmid::generatePlanOutputs(const int plasmidIndex){

    std::vector<std::vector<Operon*> >::iterator iter_plan = plans.begin();
    std::vector<std::vector<Operon*> >::const_iterator iter_plan_end = plans.end();
    int planIndex = 1;
    while (iter_plan != iter_plan_end) {
        std::ofstream plasmidPlanFile;
        std::stringstream plasmidPlanFileName;
        plasmidPlanFileName << OUTPUT_PATH << "Plasmid_" << plasmidIndex <<"_Plan_" << planIndex << "_.txt";
        plasmidPlanFile.open(plasmidPlanFileName.str().c_str());
        
        //  write this plan into a file
        std::vector<Operon*>::iterator iter_operon = iter_plan -> begin();
        std::vector<Operon*>::const_iterator iter_operon_end = iter_plan -> end();
        while (iter_operon != iter_operon_end) {
            plasmidPlanFile << (*iter_operon) -> description();
                
            iter_operon++;
        }
        
        plasmidPlanFile.close();
        
        iter_plan++;
        planIndex++;
    }
ndl;
    
}
    
    
    
    
}// namespace ustc
