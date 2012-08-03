#include "consts.h"
#include "networkinference.h"
#include "buildplasmids.h"

int main (void) {

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



