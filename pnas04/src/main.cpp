#include "consts.h"
#include "networkinference.h"

int main (void) {

    time_t start, end;
    time(&start);
    
    ustc::NetworkInference networkInference;
    networkInference.reverseEngineering();
    
    time(&end);
    double duration = difftime(end, start);
    
    std::cout << "\nIt tooks " << duration << " seconds.\n";
    return 0;
}



