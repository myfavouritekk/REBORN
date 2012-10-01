#include "consts.h"
#include "networkinference.h"
#include "buildplasmids.h"
#include <stdlib.h>


//      Total evolution number of the program
int total_evo = 1000;

//      Total number of cells
int population = 200;

//      Number of cells that do not change topology
int cells_unchanged = 5;

//      Number of cells that are generated to SBMLModels
int num_sbmlmodel = 10;


int main (int argc, char *argv[]) {
    
    //  global constants
    //  default values
    


    
    for (int i = 1; i < argc; i++) {
        if ('-' == argv[i][0]) {
            switch (argv[i][1]) {
                
                case 'e':
                case 'E':
                    total_evo = atoi(argv[++i]);
                    break;
                    
                case 'p':
                case 'P':
                    population = atoi(argv[++i]);
                    break;
                
                case 'u':
                case 'U':
                    cells_unchanged = atoi(argv[++i]);
                    break;
                    
                case 'b':
                case 'B':
                    num_sbmlmodel = atoi(argv[++i]);
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



