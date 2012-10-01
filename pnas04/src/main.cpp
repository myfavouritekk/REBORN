#include "consts.h"
#include "networkinference.h"
#include "buildplasmids.h"
#include <stdlib.h>


//  global constants
//  default values
//  
//      Total evolution number of the program
int total_evo = 1000;

//      Total number of cells
int population = 200;

//      Number of cells that do not change topology
int cells_unchanged = 5;

//      Number of cells that are generated to SBMLModels
int num_sbmlmodel = 10;


void printHelp() {
    std::cout << "Usage: reborn_cl [options]\n\n-f input_file_name\n-e total_evo\t\t-p population\n-u cells_unchanged\t-b num_sbmlmodel" << std::endl;
    std::cout << "Example: reborn_cl -e 1000 -p 200 -u 5 -b 10" << std::endl;
}




int main (int argc, char *argv[]) {
    
    std::string inputfilename;

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

                case 'f':
                case 'F':
                    inputfilename = string(argv[++i]);
                    break;

                default:
                    printHelp();
                    return 0;
            }
        }
        else {
            printHelp();
            return 0;
        }
    }
    
    if (total_evo == 0 || population == 0 || num_sbmlmodel == 0)
    {
        printHelp();
        return 0;
    }

    time_t start, end;
    time(&start);
    
    ustc::NetworkInference networkInference;
    networkInference.reverseEngineering(inputfilename);
    
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



