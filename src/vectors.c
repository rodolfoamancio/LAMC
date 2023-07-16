/* ***************************************************************************************************************
 * This code was written as a part of LAMC, a program to perform molecular Monte Carlo simulations of linear
 * alkanes on multiple ensembles.
 * 
 * This module vectors.c is a code file for vectors functions
 * 
 * Developer: Rodolfo J Amancio (rodolfojamancio@gmail.com)
 * University of Campinas, School of Chemical Engineering, Campinas/SP -  Brazil
 * 
 * **************************************************************************************************************/

#include"vectors.h"

DOUBLE_MATRIX_3X3 StrainDerivativeTensor;
BOX SimulationBox;

/***************************************************************
 * Name       | Norm
 * -------------------------------------------------------------
 * Function   | Calculate the norm of a 3D vector.
 * Parameters | r: The input vector (type: VECTOR)
 * Returns    | The norm of the input vector (type: double)
 ****************************************************************/

double Norm(VECTOR r){
  return sqrt(Squared(r.x) + Squared(r.y) + Squared(r.z));
}

/***************************************************************
 * Name       | NormSquared
 * -------------------------------------------------------------
 * Function   | Calculate the square of the norm of a 3D vector.
 * Parameters | r: The input vector (type: VECTOR)
 * Returns    | The square of the norm of the input vector (type: double)
 ****************************************************************/
double NormSquared(VECTOR r){
  return (Squared(r.x) + Squared(r.y) + Squared(r.z));
}

/***************************************************************
 * Name       | CalculateDistance
 * -------------------------------------------------------------
 * Function   | Calculate the distance between two 3D vectors.
 * Parameters | r1: The first input vector (type: VECTOR)
 *              r2: The second input vector (type: VECTOR)
 * Returns    | The distance between the two input vectors (type: double)
 ****************************************************************/
double CalculateDistance(VECTOR r1, VECTOR r2){
  VECTOR rDifference;
  rDifference = VectorSubtraction(r1, r2);
  return Norm(rDifference);
}

/***************************************************************
 * Name       | 
 * -------------------------------------------------------------
 * Function   | 
 * Parameters | 
 * Returns    | 
 ****************************************************************/
double DotProduct(VECTOR r1, VECTOR r2){
  return r1.x*r2.x + r1.y*r2.y + r1.z*r2.z;
}

/***************************************************************
 * Name       | DotProduct
 * -------------------------------------------------------------
 * Function   | Calculate the dot product of two 3D vectors.
 * Parameters | r1: The first input vector (type: VECTOR)
 *              r2: The second input vector (type: VECTOR)
 * Returns    | The dot product of the two input vectors (type: double)
 ****************************************************************/
double InternalAngle(VECTOR r1, VECTOR r2){
  double r1Norm = Norm(r1);
  double r2Norm = Norm(r2);
  double angle = 0;

  if((r1Norm != 0) && (r2Norm != 0)){
    double dotProductPerNormProduct = DotProduct(r1, r2)/(r1Norm*r2Norm);

    if(dotProductPerNormProduct < -1){
      angle = M_PI;
    }else if(dotProductPerNormProduct > 1){
      angle = 0;
    }else{
      angle = acos(dotProductPerNormProduct);
    }

  }

  return angle;
}

/***************************************************************
 * Name       | VectorSubtraction
 * -------------------------------------------------------------
 * Function   | Subtract two 3D vectors element-wise.
 * Parameters | r1: The first vector (type: VECTOR)
 *              r2: The second vector (type: VECTOR)
 * Returns    | The resulting vector after subtracting r2 from r1 (type: VECTOR)
 ****************************************************************/

VECTOR VectorSubtraction(VECTOR r1, VECTOR r2){
  VECTOR rDifference;
  rDifference.x = r1.x - r2.x;
  rDifference.y = r1.y - r2.y;
  rDifference.z = r1.z - r2.z;
  return rDifference;
}

/***************************************************************
 * Name       | VectorSum
 * -------------------------------------------------------------
 * Function   | Add two 3D vectors element-wise.
 * Parameters | r1: The first vector (type: VECTOR)
 *              r2: The second vector (type: VECTOR)
 * Returns    | The resulting vector after adding r1 and r2 together (type: VECTOR)
 ****************************************************************/
VECTOR VectorSum(VECTOR r1, VECTOR r2){
  VECTOR rSum;
  rSum.x = r1.x + r2.x;
  rSum.y = r1.y + r2.y;
  rSum.z = r1.z + r2.z;
  return rSum;
}

/***************************************************************
 * Name       | CrossProduct
 * -------------------------------------------------------------
 * Function   | Calculate the cross product of two 3D vectors.
 * Parameters | r1: The first vector (type: VECTOR)
 *              r2: The second vector (type: VECTOR)
 * Returns    | The resulting vector after calculating the 
 *              cross product of r1 and r2 (type: VECTOR)
 ****************************************************************/

VECTOR CrossProduct(VECTOR r1, VECTOR r2){
  VECTOR VectorCrossProduct;
  VectorCrossProduct.x =  r1.y*r2.z - r1.z*r2.y;
  VectorCrossProduct.y = -r1.x*r2.z + r1.z*r2.x;
  VectorCrossProduct.z =  r1.x*r2.y - r1.y*r2.x;
  return VectorCrossProduct;
}

/***************************************************************
 * Name       | NormalizeVector
 * -------------------------------------------------------------
 * Function   | Normalize a 3D vector.
 * Parameters | r: The vector to be normalized (type: VECTOR)
 * Returns    | The normalized vector (type: VECTOR)
 ****************************************************************/
VECTOR NormalizeVector(VECTOR r){
  VECTOR Normalized;
  double rNorm = Norm(r);
  Normalized.x = r.x/rNorm;
  Normalized.y = r.y/rNorm;
  Normalized.z = r.z/rNorm;
  return Normalized;
}

/***************************************************************
 * Name       | MultiplyVectorScalar
 * -------------------------------------------------------------
 * Function   | Multiply a vector by a scalar.
 * Parameters | r: The vector to be multiplied (type: VECTOR)
 *            | a: The scalar to multiply with (type: double)
 * Returns    | The resulting vector (type: VECTOR)
 ****************************************************************/
VECTOR MultiplyVectorScalar(VECTOR r, double a){
  VECTOR Result;
  Result.x = a*r.x;
  Result.y = a*r.y;
  Result.z = a*r.z;
  return Result;
}

/***************************************************************
 * Name       | RandomUnitVector
 * -------------------------------------------------------------
 * Function   | Generate a random unit vector.
 * Parameters | None
 * Returns    | The generated random unit vector (type: VECTOR)
 ****************************************************************/
VECTOR RandomUnitVector(){
  VECTOR RandomVector;
  double ransq = 2., ran1, ran2;
  while (ransq >= 1){
    ran1 = GetRandomDoubleInterval(-1, 1);
    ran2 = GetRandomDoubleInterval(-1, 1);
    ransq = Squared(ran1) + Squared(ran2);
  }
  double ranh = 2.*sqrt(1.-ransq);
  RandomVector.x = ran1 * ranh;
  RandomVector.y = ran2 * ranh;
  RandomVector.z = (1. - 2.*ransq);
  return RandomVector;
}

/***************************************************************
 * Name       | MatrixTrace
 * -------------------------------------------------------------
 * Function   | Calculate the trace of a 3x3 matrix.
 * Parameters | Matrix: The input 3x3 matrix (type: DOUBLE_MATRIX_3X3)
 * Returns    | The sum of the diagonal elements of the matrix (type: double)
 ****************************************************************/
double MatrixTrace(DOUBLE_MATRIX_3X3 Matrix){
  return Matrix.xx + Matrix.yy + Matrix.zz;
}

/***************************************************************
 * Name       | ApplyPeriodicBoundaryConditionsVector
 * -------------------------------------------------------------
 * Function   | Apply periodic boundary conditions to a vector.
 * Parameters | r: The input vector (type: VECTOR)
 * Returns    | The transformed vector after applying the 
 *              periodic boundary conditions (type: VECTOR)
 ****************************************************************/
VECTOR ApplyPeriodicBoundaryConditionsVector(VECTOR r){
  VECTOR rTransformed = r;
  rTransformed.x = r.x - SimulationBox.xSize*round(r.x/SimulationBox.xSize);
  rTransformed.y = r.y - SimulationBox.ySize*round(r.y/SimulationBox.ySize);

  if(!SimulationBox.ClosedBox)
    rTransformed.z = r.z - SimulationBox.zSize*round(r.z/SimulationBox.zSize);
  
  return rTransformed;
}

/***************************************************************
 * Name       | GetRandomPosition
 * -------------------------------------------------------------
 * Function   | Generate a random position vector within the 
 *              simulation box.
 * Parameters | None
 * Returns    | The randomly generated position vector 
 *              (type: VECTOR)
 ****************************************************************/
VECTOR GetRandomPosition(void){
  VECTOR RandomPosition;
  RandomPosition.x = GetRandomDoubleInterval(SimulationBox.xMin, SimulationBox.xMax);
  RandomPosition.y = GetRandomDoubleInterval(SimulationBox.yMin, SimulationBox.yMax);
  RandomPosition.z = GetRandomDoubleInterval(SimulationBox.zMin, SimulationBox.zMax);
  return RandomPosition;
}