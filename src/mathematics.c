/* ***************************************************************************************************************
 * This code was written as a part of LAMC, a program to perform molecular Monte Carlo simulations of linear
 * alkanes on multiple ensembles.
 * 
 * This module mathematics.c is a code file for defining mathematical functions
 * 
 * Developer: Rodolfo J Amancio (rodolfojamancio@gmail.com)
 * University of Campinas, School of Chemical Engineering, Campinas/SP -\sBrazil
 * 
 * **************************************************************************************************************/

// ---------------------------------------------------------------------------------------------------------------
// Headers loading
// ---------------------------------------------------------------------------------------------------------------

#include"mathematics.h"

// ---------------------------------------------------------------------------------------------------------------


/***************************************************************
 * Name       | LimMinMax
 * -------------------------------------------------------------
 * Function   | Compares a given value against a minimum and
 *              maximum value and returns the checked value.
 * Parameters |
 *              - val: The value to be checked against the 
 *                     minimum and maximum values.
 *              - min: The minimum value that val can be.
 *              - max: The maximum value that val can be.
 * Returns    | 
 *              The double value that has been checked against 
 *              the minimum and maximum values.
 ****************************************************************/
double LimMinMax(double val, double min, double max){
  if(val > max){
    return max;
  }else if(val < min){
    return min;
  }else{
    return val;
  }
}

/***************************************************************
 * Name       | GetRandomNumber
 * -------------------------------------------------------------
 * Function   | Generates a random number between 0 and 1.
 * Parameters | None
 * Returns    | The generated random number as a double.
 ****************************************************************/
double GetRandomNumber(){
  return ((double) rand())/((double) RAND_MAX);
}

/***************************************************************
 * Name       | GetRandomNumber
 * -------------------------------------------------------------
 * Function   | Generates a random number between 0 and 1.
 * Parameters | None
 * Returns    | The generated random number as a double.
 ****************************************************************/
double GetRandomDoubleInterval(double min, double max){
  return ((double) rand())/((double)(RAND_MAX/(max-min))) + min;
}

/***************************************************************
 * Name       | GetRandomIntegerInterval
 * -------------------------------------------------------------
 * Function   | Generates a random integer within a specified 
 *              interval.
 * Parameters |
 *              - lower: The lower bound of the interval.
 *              - upper: The upper bound of the interval.
 * Returns    | The generated random integer within the specified 
 *              interval.
 ****************************************************************/
int GetRandomIntegerInterval(int lower, int upper){
  return (rand() % (upper - lower + 1)) + lower;
}

/***************************************************************
 * Name       | GetRandomGaussianNumber
 * -------------------------------------------------------------
 * Function   | Generates a random number following a Gaussian 
 *              distribution.
 * Parameters | None.
 * Returns    | A random number following a Gaussian distribution.
 ****************************************************************/
double GetRandomGaussianNumber(){
  double ran1, ran2, r2 = 2;
  while(r2 >= 1 || r2 == 0){
    ran1 = GetRandomDoubleInterval(-1, 1);
    ran2 = GetRandomDoubleInterval(-1, 1);
    r2 = Squared(ran1) + Squared(ran2);
  }
  return ran1*sqrt(-2.0*log(r2)/r2);
}

/***************************************************************
 * Name       | Metropolis
 * -------------------------------------------------------------
 * Function   | Determines whether a movement should be accepted 
 *              using the Metropolis algorithm based on a given 
 *              probability ratio.
 * Parameters | ProbabilityRatio: The probability ratio used to 
 *              determine whether the movement is accepted.
 * Returns    | True if the movement is accepted, false otherwise.
 ****************************************************************/
bool Metropolis(double ProbabilityRatio){
  bool movementAccepted = false;
  if(ProbabilityRatio >= 1){
    movementAccepted =  true;
  }else if(ProbabilityRatio > 0){
    movementAccepted = GetRandomNumber() < ProbabilityRatio ? true : false;
  }
  return movementAccepted;
}