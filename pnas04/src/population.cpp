#include "population.h"


#include "cell.h"
#define HALFTOTAL 100


//for growth phase
void Population::growth(){;
    Cell* currCell;
    for(int i=0; i < HALFTOTAL; i++){
        cells[HALFTOTAL+i]=cells[i];    
        currCell = &cells[TOTAL+i];
        currCell->mutation();
        currCell->getScore(&sfunc, targetData);
    } 
}


//for seclction phase
void Population::selection(){

    int i,j,k;
    /*
    double maxScore=0;
    int maxPosition;
    Cell temp;
    Cell* currCell;
    for(i=0; i<HALFTOTAL; i++){
        currCell = &cells[i];
        maxScore = 0;
        for(j=i; j<HALFTOTAL*2,j++){
            if(currCell->getCurrScore()>maxScore){
                maxScore = currCell->getCurrScore;
                maxPosition = j;
            }
            currCell++;
        }
        temp = cells[i];
        cells[i] = cells[maxPosition];
        cells[maxPosition] = temp;
    }
    */
    quickSort(cells,HALFTOTAL*2);

}

void quickSort(Cell* cells[],int num){
    
    int i;
    int num1=0,num2=0;//to store lengths of two subgroups
    Cell* less[num],greater[num];//to store pointers of subgroups
    if(num<=1)return;
    for(i=1;i<num;i++){
        if(cells[i]->getCurrScore()>cells[0]->getCurrScore()){
            greater[num2]=cells[i];
            num2++;
        }else{
            less[num1]=cells[i];
            num1++;
        }
    }


    //recursion
    quickSort(less,num1);
    quickSort(greater,num2);

    //final assignments
    cells[num1]=cells[0];
    for(i=0;i<num1;i++){
        cells[i]=less[i];
    }
    for(i=0;i+num1+1<num;i++){
        cells[i+num1+1]=greater[i];
    }
}


