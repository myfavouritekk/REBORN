#include "scorefunction.h"


//Constructor
ScoreFunc::ScoreFunc(int _type){

    type = _type;

}

ScoreFunc::~ScoreFunc(){
    
}



//return score based on the current "type"
double ScoreFunc::getScore(double* data1, double* data2, int n){

    switch (type) {
    case 0:
        return getScore_0(data1, data2, n);
    case 1:
        return getScore_1(data1, data2, n);
    /*
     * to be contined...
     *
    case 2:
        return getScore_2(data1, data2, n);
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
double ScoreFunc::getScore_0(double* data1, double* data2,int  n){   


    double sum = 0;
    for(int i = 0; i < n; i++){
        
        sum += (*data1 - *data2)*(*data1 - *data2);
        data1++;
        data2++;

    }

    return sum;

}


//a score function based on the absolute distance
double ScoreFunc::getScore_1(double* data1, double* data2,int  n){
    

    double sum = 0;
    for(int i = 0; i < n; i++){
        
        sum += fabs(*data1-*data2);
        data1++;
        data2++;

    }

    return sum;

}

/*to be continued...
 *
double ScoreFunc::getScore_2(double* data1, double* data2,int  n){
    
}
double ScoreFunc::getScore_3(double* data1, double* data2,int  n){
    
}
double ScoreFunc::getScore_4(double* data1, double* data2,int  n){
    
}
double ScoreFunc::getScore_5(double* data1, double* data2,int  n){
}

*/
