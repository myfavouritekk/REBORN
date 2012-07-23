#include "scorefunction.h"

namespace ustc{

//Constructor
ScoreFunc::ScoreFunc(int _type){

    type = _type;

}

ScoreFunc::~ScoreFunc(){
    
}



//return score based on the current "type"
double ScoreFunc::getScore(double* targetData, double* generatedData, int n){

    switch (type) {
    case 0:
        return getScore_0(targetData, generatedData, n);
    case 1:
        return getScore_1(targetData, generatedData, n);
            
    case 2:
        return getScore_2(targetData, generatedData, n);
    
    /*
     
    case 3:
        return getScore_3(data1, data2, n);
    case 4:
        return getScore_4(data1, data2, n);
    case 5:
        return getScore_5(data1, data2, n);
     *
     */
    }

    
    return 0;
}

//a score function based on the square of distance
double ScoreFunc::getScore_0(double* targetData, double* generatedData,int  n){   


    double sum = 0;
    for(int i = 0; i < n; i++){
        
        sum += (*targetData - *generatedData)*(*targetData - *generatedData);
        targetData++;
        generatedData++;

    }

    return sum;

}


//a score function based on the absolute distance
double ScoreFunc::getScore_1(double* targetData, double* generatedData,int  n){
    

    double sum = 0;
    for(int i = 0; i < n; i++){
        
        sum += fabs(*targetData-*generatedData);
        targetData++;
        generatedData++;

    }

    return sum;

}

#define RATIO 1.0
double ScoreFunc::getScore_2(double* targetData, double* generatedData,int  n){
    double targetMax = 0, targetMin = 0;
    double generatedMax = 0, generatedMin = 0;
    for (int i = 0; i < n; i++) {
        if (generatedData[i] > generatedMax) {
            generatedMax = generatedData[i];
        }
        if (generatedData[i] < generatedMin) {
            generatedMin = generatedData[i];
        }
        if (targetData[i] > targetMax) {
            targetMax = targetData[i];
        }
        if (targetData[i] < targetMin) {
            targetMin = targetData[i];
        }
    }
    
    //normalize generatedData so that it has same scale as targetData
    double sum = 0;
    double ratio = (targetMax - targetMin) / (generatedMax - generatedMin);
    double intercept = - (targetMax * generatedMin - targetMin * generatedMax) /  (generatedMax - generatedMin);
    for (int i = 0; i < n; i++) {
        double normal = generatedData[i] * ratio + intercept; // scale the generatedData
        sum += fabs(targetData[i] - normal);
        sum += RATIO * fabs(targetData[i] - generatedData[i]);
    }
    return sum;
}


double ScoreFunc::getScore_3(double* targetData, double* generatedData,int  n){
    
    double sum = 0;
    for (int i = 0; i < n - 1; i++) {
        double targetChange = targetData[i + 1] - targetData[i];
        double generateChange = generatedData[i + 1] - targetData[i];
        sum += RATIO * fabs(targetChange - generateChange);
        sum += fabs(targetData[i] - generatedData[i]);
    }
    return sum;
    
}

/*to be continued...
 *
double ScoreFunc::getScore_4(double* data1, double* data2,int  n){
    
}
double ScoreFunc::getScore_5(double* data1, double* data2,int  n){
}

*/


}   //namespace ustc