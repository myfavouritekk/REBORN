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
	askInformation();
	
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
	std::cout << "Do you know more information about the reaction? <y/n>" << std::endl;
	if(getchar() == 'y'){
		std::cout << targetPop -> cells[0] -> 
	}
