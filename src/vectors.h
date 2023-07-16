/* ***************************************************************************************************************
 * This code was written as a part of LAMC, a program to perform molecular Monte Carlo simulations of linear
 * alkanes on multiple ensembles.
 * 
 * This module vectors.h is a header file for vectors structs and functions
 * 
 * Developer: Rodolfo J Amancio (rodolfojamancio@gmail.com)
 * University of Campinas, School of Chemical Engineering, Campinas/SP -  Brazil
 * 
 * **************************************************************************************************************/

// ---------------------------------------------------------------------------------------------------------------
// Headers loading
// ---------------------------------------------------------------------------------------------------------------

#ifndef VECTORS_H
#define VECTORS_H
#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<stdbool.h>
#include<math.h>
#include"mathematics.h"

// ---------------------------------------------------------------------------------------------------------------
// Structures definitions
// ---------------------------------------------------------------------------------------------------------------

typedef struct{
  double x;
  double y;
  double z;
} VECTOR;

typedef struct
{
    double xx;
    double xy;
    double xz;
    double yx;
    double yy;
    double yz;
    double zx;
    double zy;
    double zz;
} DOUBLE_MATRIX_3X3;

typedef struct{
  double xMin;
  double yMin;
  double zMin;
  double xMax;
  double yMax;
  double zMax;
  double xSize;
  double ySize;
  double zSize;
  double volume;
  bool   ClosedBox;
} BOX;


// ---------------------------------------------------------------------------------------------------------------
// Vector functions
// ---------------------------------------------------------------------------------------------------------------

double Norm(VECTOR r);
double NormSquared(VECTOR r);
double CalculateDistance(VECTOR r1, VECTOR r2);
double DotProduct(VECTOR r1, VECTOR r2);
double InternalAngle(VECTOR r1, VECTOR r2);
VECTOR VectorSubtraction(VECTOR r1, VECTOR r2);
VECTOR VectorSum(VECTOR r1, VECTOR r2);
VECTOR CrossProduct(VECTOR r1, VECTOR r2);
VECTOR NormalizeVector(VECTOR r);
VECTOR MultiplyVectorScalar(VECTOR r, double a);
VECTOR RandomUnitVector();
VECTOR ApplyPeriodicBoundaryConditionsVector(VECTOR r);
VECTOR GetRandomPosition(void);

// ---------------------------------------------------------------------------------------------------------------
// Matrix functions
// ---------------------------------------------------------------------------------------------------------------

double MatrixTrace(DOUBLE_MATRIX_3X3 Matrix);

// ---------------------------------------------------------------------------------------------------------------
// Globals
// ---------------------------------------------------------------------------------------------------------------

extern DOUBLE_MATRIX_3X3 StrainDerivativeTensor;
extern BOX SimulationBox;

#endif
