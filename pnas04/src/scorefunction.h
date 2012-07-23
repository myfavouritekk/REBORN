#ifndef SCOREFUNC_H
#define SCOREFUNC_H


#include "consts.h"
#include <cmath>

namespace ustc{


class ScoreFunc{
    
    public:
        
        //Constructor
        ScoreFunc(int type);
       
        //Destructor
        ~ScoreFunc();


        //getScore function to calculate difference between data1 and data2
        double getScore(double* targetData, double* generatedData,int n);

    private:

        int type;//0-5 for different types

        //several Score function options
        double getScore_0(double* targetData, double* generatedData, int n);
        double getScore_1(double* targetData, double* generatedData, int n);
        double getScore_2(double* targetData, double* generatedData, int n);
        double getScore_3(double* targetData, double* generatedData, int n);
        double getScore_4(double* targetData, double* generatedData, int n);
        double getScore_5(double* targetData, double* generatedData, int n);

	};


}   //namespace ustc

#endif
