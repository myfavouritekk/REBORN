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
int** findMatrixRecursion(
                           const int** databaseMatrix,
                           const int** targetMatrix,
                           const int* choicesPool,
                           const int numberOfChoicesInPool,
                           const int numberOfChoicesToBeChosen,
                           int* numberOfPossibleChoiceSets
                           )
{
    int** findMatrixRecursion(
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
        int** possibleSolutions = findMatrixRecursion(
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

	//==============================================================================//
    //                                                                              //
    //  findMatrixRecurssion2 method                                                //
    //  purpose: find small matrix in a big matrix                                  //
    //  it uses recurssion algorithm to find MxM matrix in NxN matrix, based on     //
    //      recussion of finding (M-1)x(M-1) matrix in (N-1)x(N-1) matrix           //
    //  parameters:                                                                 //
    //          databaseMatrix is the global giant matrix                           //
    //          targetMatrix is the current matrix to be find						//
	//																				//
    //          choicesPoolOfRows contains the indice of rows can be chosen			//
	//			choicesPoolOfColumns contains the indice of columns can be chosen	//	
    //																				//
	//			numberOfRowChoicesInPool is the size of choicesPool of rows			//
	//			numberOfColumnChoicesInPool is the size of choicesPool of columns	//
    //																				//
	//			numberOfChoicesToBeChosen is the size of targetMatrix               //
    //          numberOfPossibleChoicesSets is the number of possible solutions     //
    //  return value:                                                               //
    //          return value is the 2-d array that contains the possible solutions  //
    //                                                                              //
    //==============================================================================//

int*** findMatrixRecursion2	(int** databaseMatrix,
                             int** targetMatrix,
                             int* choicesPoolOfRows,
                             int* choicesPoolOfColumns,
                             int numberOfRowChoicesInPool,
                             int numberOfColmunsChoicesInPool,
                             int numberOfChoicesToBeChosen,
                             int* numberOfPossibleChoices
                             )
{
	int*** findMatrixRecursion2(int** databaseMatrix,
								int** targetMatrix,
								int* choicesPoolOfRows,
								int* choicesPoolOfColumns,
								int numberOfRowChoicesInPool,
								int numberOfColmunsChoicesInPool,
								int numberOfChoicesToBeChosen,
								int* numberOfPossibleChoices
								);
	if(numberOfChoicesToBeChosen == 1){
		int numPossible = 0;
		std::vector<int**> choicesVector;
		


		for(int i = 0; i < numberOfRowChoicesInPool; i ++){
			for(int j = 0; j < numberOfColmunsChoicesInPool; j ++){
				if(databaseMatrix[choicesPoolOfRows[i]][choicesPoolOfColumns[j]] == 2 ||
                   databaseMatrix[choicesPoolOfRows[i]][choicesPoolOfColumns[j]] == targetMatrix[0][0]){
                    numPossible ++;
                    int** aChoice = new int*[2];
                    for(int i = 0; i < 2; i ++){
                        aChoice[i] = new int[1];
                    }
                    aChoice[0][0] = choicesPoolOfRows[i];
                    aChoice[1][0] = choicesPoolOfColumns[j];
                    choicesVector.push_back(aChoice);
				}
			}
		}
		int *** choices = new int** [numPossible];
		for(int i = 0; i < numPossible; i ++ ){
			choices[i] = choicesVector[i];
		}
        *numberOfPossibleChoices = numPossible;
		return choices;
	}
    
    //      variables for recurssion
    int** newTargetMatrix = new int*[numberOfChoicesToBeChosen - 1];
    int numOfWorkingChoices = 0;
    std::vector<int**> choicesVector;
    
    
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
	int firstElement [2][1];
    for(int i = 0; i < numberOfRowChoicesInPool; i++){
		for(int j = 0; j < numberOfColmunsChoicesInPool; j ++ ){
			if(databaseMatrix[choicesPoolOfRows[i]][choicesPoolOfColumns[j]] != 2 &&
               databaseMatrix[choicesPoolOfRows[i]][choicesPoolOfColumns[j]] != targetMatrix[0][0]){
				continue;
			}
			firstElement[0][0] = choicesPoolOfRows[i];
			firstElement[1][0] = choicesPoolOfColumns[j];
            
	        int* newPoolOfRow = new int[numberOfRowChoicesInPool - 1];
			for(int h = 0, k = 0; h < numberOfRowChoicesInPool; h ++){
				if( i != h){
					newPoolOfRow[k] = choicesPoolOfRows[h];
					k++;
				}
			}

			
			int* newPoolOfColumn = new int[numberOfColmunsChoicesInPool - 1];
			for(int g = 0, k = 0; g < numberOfColmunsChoicesInPool; g ++){
				if( j != g){
					newPoolOfColumn[k] = choicesPoolOfColumns[g];
					k++;
				}
			}
            
			
			int numberOfNewSets;
			int*** possibleSolutions = findMatrixRecursion2(
                                                            databaseMatrix,
                                                            newTargetMatrix,
                                                            newPoolOfRow,
                                                            newPoolOfColumn,
                                                            numberOfRowChoicesInPool - 1,
                                                            numberOfColmunsChoicesInPool - 1,
                                                            numberOfChoicesToBeChosen - 1,
                                                            &numberOfNewSets
                                                            );
            //  judge whether the solutions work
			for(int i = 0; i < numberOfNewSets; i ++){
				bool work = true;
				for(int j = 0; j < numberOfChoicesToBeChosen - 1; j ++){
					if(databaseMatrix[firstElement[0][0]][possibleSolutions[i][1][j]] != 2 &&
                       targetMatrix[0][j + 1] != databaseMatrix[firstElement[1][0]][possibleSolutions[i][1][j]] &&
                       databaseMatrix[possibleSolutions[i][0][j]][firstElement[1][0]] != 2 &&
                       targetMatrix[j + 1][0] != databaseMatrix[possibleSolutions[i][0][j]][firstElement[1][0]] )
					   {
						work = false;
						break;
					}
				}
				if(work){
					numOfWorkingChoices ++;
					int** aChoice = new int* [2];
					aChoice[0] = new int [numberOfChoicesToBeChosen];
					aChoice[1] = new int [numberOfChoicesToBeChosen];
					aChoice[0][0] = firstElement[0][0];
					aChoice[1][0] = firstElement[1][0];
					for(int k = 1; k < numberOfChoicesToBeChosen; k++ ){
						aChoice[0][k] = possibleSolutions[i][0][k - 1];
						aChoice[1][k] = possibleSolutions[i][1][k - 1];
					}
					choicesVector.push_back (aChoice);
				}
			}
		}
    }
    *numberOfPossibleChoices = numOfWorkingChoices;
    int ***choices = new int**[numOfWorkingChoices];
	for(int i = 0 ; i < numOfWorkingChoices; i ++){
		choices[i] = choicesVector[i];
	}
	return choices;
    
}
    
    
    
    
void Plasmid::findCompleteCandidates(
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
    candidateIndice = findMatrixRecursion(database, (const int**)wholeRegulatoryMatrix, choicePool, numColumnOfDatabase, numOfGenes, &numberOfCandidates);
    
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
   
void Plasmid::findCompleteCandidatesUsingBiobricks(
													const int numRow,
													const int numcolumn,
													const std::string* namesOfRegulators,
													const std::string* namesOfRegulatees,
													const int** database)
{
	int*** candidateIndice;
	int numberOfCandidates;
	int *choicePoolOfRows = new int[numRow];
	int *choicePoolOfColumns = new int[numcolumn];

	for(int i = 0 ; i < numRow; i++){
		choicePoolOfRows[i] = i;
	}
	for(int i = 0 ; i < numcolumn; i++){
		choicePoolOfColumns[i] = i;
	}

	candidateIndice = findMatrixRecursion2(database, wholeRegulatoryMatrix,choicePoolOfRows,choicePoolOfColumns,numRow,numcolumn,numOfGenes,&numberOfCandidates);
	biobrickPlans.clear();

	for(int i = 0;i < numberOfCandidates;i++){
		std::vector<BioBrick*> aPlan;
		for(int j = 0; j < numOfGenes; j++){ 
				BioBrick *aNewBiobrick = new BioBrick(namesOfRegulators[candidateIndice[i][1][j]],namesOfRegulatees[candidateIndice[i][0][j]]);
			    aPlan.push_back(aNewBiobrick);
		}
		biobrickPlans.push_back(aPlan);
	}
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
        plasmidPlanFileName << OUTPUT_PATH << "Plasmid_" << plasmidIndex + 1 <<"_Plan_" << planIndex << "_.txt";
        plasmidPlanFile.open(plasmidPlanFileName.str().c_str());
        
        //  write this plan into a file
        std::vector<Operon*>::iterator iter_operon = iter_plan -> begin();
        std::vector<Operon*>::const_iterator iter_operon_end = iter_plan -> end();
        while (iter_operon != iter_operon_end) {
            plasmidPlanFile << (*iter_operon) -> description();
                
            iter_operon++;
        }
        
        plasmidPlanFile.close();
        std::cout << "Finished plan " << planIndex << std::endl;
        iter_plan++;
        planIndex++;
    }
    
    //  print status  in the debug area
    std::cout
        << "Finished Generating Plasmid Plans for: Plasmid " << plasmidIndex + 1
        << "\t//\t" << plans.size() << " plans in all"
        << std::endl;
    
}
    
    
    
    
}// namespace ustc
