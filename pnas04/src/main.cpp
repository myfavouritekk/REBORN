#include "consts.h"
#include "networkinference.h"
#include "buildplasmids.h"
#include <stdlib.h>

int main (int argc, char *argv[]) {
    
    //  global constants
    //  default values
    
    //      Total evolution number of the program
    TOTAL_EVO = 1000;
    
    //      Total number of cells
    POPULATION = 200;
    
    //      Number of cells that do not change topology
    CELLS_UNCHANGED = 5;
    
    //      Number of cells that are generated to SBMLModels
    NUM_SBMLMODEL = 10;

    
    for (int i = 1; i < argc; i++) {
        if ('-' == argv[i][0]) {
            switch (argv[i][1]) {
                
                case 'e':
                case 'E':
                    TOTAL_EVO = atoi(argv[++i]);
                    break;
                    
                case 'p':
                case 'P':
                    POPULATION = atoi(argv[++i]);
                    break;
                
                case 'u':
                case 'U':
                    CELLS_UNCHANGED = atoi(argv[++i]);
                    break;
                    
                case 'b':
                case 'B':
                    NUM_SBMLMODEL = atoi(argv[++i]);
                    break;
                    
                default:
                    break;
            }
        }
    }
    
    

    time_t start, end;
    time(&start);
    
    ustc::NetworkInference networkInference;
    networkInference.reverseEngineering();
    
    time(&end);
    double NIDuration = difftime(end, start);
    
    std::cout << "\nIt tooks " << NIDuration << " seconds.\n";
    
    time(&start);
    
    ustc::BuildPlasmids buildingPlasmids;
    buildingPlasmids.buildProcess();
    
    time(&end);
    double BPDuration = difftime(end, start);
    
    std::cout << "\nIt tooks " << BPDuration << " seconds.\n";
    
    return 0;
}



