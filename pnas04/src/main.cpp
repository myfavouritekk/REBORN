#include "population.h"
#include <ctime>
#include "consts.h"

int main (void) {

    srand(1);
    //  initialization
    ustc::Population myPop(POPULATION, TOTAL_EVO);
    myPop.init ();
	// ask users if they can input more information about the cell
	void askInformation(ustc::Population* targetPop);
	askInformation(&myPop);
	
    int i = 1;
    int sum = 1;
    //  evolution
    while (!myPop.isTerminate ()) {
        if ((TOTAL_EVO + 1 - myPop.getEvolution()) % (50 * sum) == 0) {
            //sort
            myPop.sort();
            myPop.output();
            myPop.mutation();
            i++;
            sum += i;
        }
        myPop.mut_parameters_simAnneal();
        std::cout << "Finished Evolution: " << TOTAL_EVO - myPop.getEvolution() << std::endl;
    }

    
    //  output
    myPop.sort();
    myPop.output();
    myPop.genSBMLFormat();

    return 0;
}


void askInformation(ustc::Population* targetPop){
	int rType = 0;
	
	std::cout << "Do you know more information about the reaction? <y/n>" << std::endl;
	while(getchar() == 'y'){
		if(getchar() == 'y'){
			//print the types of reaction
			std::cout<< "1. inducer binds a protein" << std::endl;
			std::cout<< "2. protein binds protein" << std::endl;
			std::cout<< "3. protein binds a gene" << std::endl;
			std::cout<< "4. protein modification" << std::endl;
			ustc::Cell* aCell = (targetPop -> getCells())[0];
			//size of nodes
			int size = (aCell->getNodesVector()) -> size();
			// print the Nstrings of nodes
			for(int i = 0; i < size; i++){
				std::stringstream num;
				num << (i + 1);
				std::cout << "index: " << num.str() << (*(aCell -> getNodesVector()))[i] -> getNstring() << std::endl;
			}
			// ask the reaction type
			std::cout<< "please type in the type of reaction" << std::endl;
			if(getchar() == '4'){
				std::cout<< "please type in the index of reactans";
				int index = ((int()(getchar()) - 1));
				rType = 3;
				for(int i = 0; i < POPULATION; i ++){
				((targetPop -> getCells())[i]) -> addReaction(rType,index);
				}
			}
			else{
				if(getchar() == '1'){
					rType = 8;
				}
				else{
					if(getchar() == '2'){
						rType = 5;
					}
					else{
						if(getchar() == '3'){
							rType = 2;
						}
					}
				}
				int index1 = 0;
				int index2 = 0;
				std::cout<< "please type in the index of reactans";
				index1 = ((int()(getchar()) - 1));
				std::cout<< "please type in the index of another reactans";
				index2 = ((int()(getchar()) - 1));
				for(int i = 0; i < POPULATION; i ++){
				((targetPop -> getCells())[i]) -> addReaction(rType,index1,index2);
				}
			}
			
		}
		else{
			if(getchar() == 'n'){
				break;
			}
		}
		std::cout << "Do you know more information about the reaction? <y/n>" << std::endl;
	}
}
