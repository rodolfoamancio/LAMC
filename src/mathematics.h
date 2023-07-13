/* ***************************************************************************************************************
 * This code was written as a part of LAMC, a program to perform molecular Monte Carlo simulations of linear
 * alkanes on multiple ensembles.
 * 
 * This module mathematics.h is a header file for defining mathematical functions and macros
 * 
 * Developer: Rodolfo J Amancio (rodolfojamancio@gmail.com)
 * University of Campinas, School of Chemical Engineering, Campinas/SP -	Brazil
 * 
 * **************************************************************************************************************/

// ---------------------------------------------------------------------------------------------------------------
// Headers loading
// ---------------------------------------------------------------------------------------------------------------

#ifndef MATHEMATICS_H
#define MATHEMATICS_H
#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<stdbool.h>
#include<math.h>

// ---------------------------------------------------------------------------------------------------------------
// Macros
// ---------------------------------------------------------------------------------------------------------------

#define Squared(x) ((x)*(x))
#define Cube(x) ((x)*(x)*(x))
#define Min2(x, y) (x < y ? x : y)
#define Max2(x, y) (x > y ? x : y)
#define Min3(x, y, z) (Min2(x, Min2(y, z)))
#define Max3(x, y, z) (Max2(x, Max2(y, z)))
#define Sign(x)   (x > 0 ? 1 : (x < 0 ? -1: 0))

// ---------------------------------------------------------------------------------------------------------------

double LimMinMax(double val, double min, double max);

double GetRandomNumber(void);
double GetRandomDoubleInterval(double min, double max);
double GetRandomGaussianNumber(void);
double GetRandomStandardNormalNumber(void);
int GetRandomIntegerInterval(int lower, int upper);

bool Metropolis(double ProbabilityRatio);

#endif
// ---------------------------------------------------------------------------------------------------------------