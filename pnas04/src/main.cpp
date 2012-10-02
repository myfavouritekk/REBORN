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



std::string output_path("output/");
std::string saves_path("saves/");
std::string html_saves_path("saves/html/");





void printHelp() {
    std::cout << "Usage: reborn_cl [options]\n\n-f input_file_name\t-o output folder\n-n no more information\n-e total_evo\t\t-p population\n-u cells_unchanged\t-b num_sbmlmodel" << std::endl;
    std::cout << "Example: reborn_cl -o Result/out -n -e 1000 -p 200 -u 5 -b 10" << std::endl;
}




int main (int argc, char *argv[]) {
    
    std::string inputfilename;
    bool isnoinfo = 0;

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

                case 'o':
                case 'O':
                    output_path = string(argv[++i]) + "/" + output_path;
                    saves_path = string(argv[i]) + "/" + saves_path;
                    html_saves_path = string(argv[i]) + "/" + html_saves_path;
                    break;

                case 'n':
                case 'N':
                    isnoinfo = 1;
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
    networkInference.reverseEngineering(inputfilename, isnoinfo);
    
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



