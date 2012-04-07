#ifndef SCOREFUNC_H
#define SCOREFUNC_H



class ScoreFunc{
    
    public:
        
        //Constructor
        ScoreFunc(int type);
       
        //Destructor
        ~ScoreFunc();


        //getScore function to calculate difference between data1 and data2
        double getScore(double* data1, double* data2,int n);

    private:

        int type;//0-5 for different types

        //several Score function options
        double getScore_0(double* data1, double* data2, int n);
        double getScore_1(double* data1, double* data2, int n);
        double getScore_2(double* data1, double* data2, int n);
        double getScore_3(double* data1, double* data2, int n);
        double getScore_4(double* data1, double* data2, int n);
        double getScore_5(double* data1, double* data2, int n);

	};




#endif
