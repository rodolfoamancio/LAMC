/* ***************************************************************************************************************
 * This code was written as a part of LAMC, a program to perform molecular Monte Carlo simulations of linear
 * alkanes on multiple ensembles.
 * 
 * This module molecule.c is a code file for defining molecule function
 * 
 * Developer: Rodolfo J Amancio (rodolfojamancio@gmail.com)
 * University of Campinas, School of Chemical Engineering, Campinas/SP -\sBrazil
 * 
 * **************************************************************************************************************/

#include"molecule.h"

/* ***************************************************************************************************************
 * Name       | GetSigmaCombination
 * ---------------------------------------------------------------------------------------------------------------
 * Function   | Calculates the average of two given sigma values.
 * Parameters |
 *              - SigmaA: The first sigma value.
 *              - SigmaB: The second sigma value.
 * Returns    | The average of SigmaA and SigmaB.
 * **************************************************************************************************************/
double GetSigmaCombination(double SigmaA, double SigmaB){
  return (SigmaA != SigmaB) ? (SigmaA + SigmaB)/2 : SigmaA;
}

/* ***************************************************************************************************************
 * Name       | GetEpsilonCombination
 * ---------------------------------------------------------------------------------------------------------------
 * Function   | Calculates the square root of the product of two given epsilon values.
 * Parameters |
 *              - EpsilonA: The first epsilon value.
 *              - EpsilonB: The second epsilon value.
 * Returns    | The square root of the product of EpsilonA and EpsilonB.
 * **************************************************************************************************************/
double GetEpsilonCombination(double EpsilonA, double EpsilonB){
  return (EpsilonA != EpsilonB) ? sqrt(EpsilonA*EpsilonB) : EpsilonA;
}

/* ***************************************************************************************************************
 * Name       | GetExponentCombination
 * ---------------------------------------------------------------------------------------------------------------
 * Function   | Applies the combination rule from Lafitte (2013) for the exponents of the Mie potential
 * Parameters |
 *              - ExponentA: The first exponent value.
 *              - ExponentA: The second exponent value.
 * Returns    | The combination exponent
 * **************************************************************************************************************/
double GetExponentCombination(double ExponentA, double ExponentB){
  return (ExponentA != ExponentB) ? sqrt((ExponentA - 3)*(ExponentB - 3)) + 3 : ExponentA;
}

enum CarbonType GetCarbonType(int ChainLength, int Position){
  if(ChainLength == 1){
    return CH4;
  }else if(ChainLength == 2){
    return CH3e;
  }else if(ChainLength > 2){
    return (Position == 0 || Position == ChainLength - 1) ? CH3 : CH2;
  }
}

/* ***************************************************************************************************************
 * Name       | GetAlkaneEpsilon
 * ---------------------------------------------------------------------------------------------------------------
 * Function   | Retrieves the epsilon value based on the given CarbonType.
 * Parameters | - Type: The CarbonType enum representing the type of alkane pseudoatom.
 * Returns    | The epsilon value for the specified alkane type.
 * **************************************************************************************************************/
double GetAlkaneEpsilon(enum CarbonType Type){
  switch(Type){
    case CH4:
      return 148*BOLTZMANN_CONSTANT;
      break;
    case CH3e:
      return 130.780*BOLTZMANN_CONSTANT;
      break;
    case CH3:
      return 136.318*BOLTZMANN_CONSTANT;
      break;
    case CH2:
      return 52.9133*BOLTZMANN_CONSTANT;
  }
}

/* ***************************************************************************************************************
 * Name       | GetAlkaneSigma
 * ---------------------------------------------------------------------------------------------------------------
 * Function   | Retrieves the sigma value based on the given CarbonType.
 * Parameters | - Type: The CarbonType enum representing the type of alkane pseudoatom.
 * Returns    | The sigma value for the specified alkane type.
 * **************************************************************************************************************/
double GetAlkaneSigma(enum CarbonType Type){
  switch(Type){
    case CH4:
      return 3.73;
      break;
    case CH3e:
      return 3.6463;
      break;
    case CH3:
      return 3.6034;
      break;
    case CH2:
      return 4.0400;
      break;
  }
}

/* ***************************************************************************************************************
 * Name       | GetAlkaneRepulsiveExponent
 * ---------------------------------------------------------------------------------------------------------------
 * Function   | Retrieves the repulsive exponent value based on the given CarbonType.
 * Parameters | - Type: The CarbonType enum representing the type of alkane pseudoatom.
 * Returns    | The repulsive exponent value for the specified alkane type.
 * **************************************************************************************************************/
double GetAlkaneRepulsiveExponent(enum CarbonType Type){
  return (Type == CH4) ? 12 : 14;
}

/* ***************************************************************************************************************
 * Name       | GetAlkaneAttractiveExponent
 * ---------------------------------------------------------------------------------------------------------------
 * Function   | Retrieves the attractive exponent value based on the given CarbonType.
 * Parameters | - Type: The CarbonType enum representing the type of alkane pseudoatom.
 * Returns    | The attractive exponent value for the specified alkane type.
 * **************************************************************************************************************/
double GetAlkaneAttractiveExponent(enum CarbonType Type){
  return 6;
}

/* ***************************************************************************************************************
 * Name       | GetAlkaneAtomMolarMass
 * ---------------------------------------------------------------------------------------------------------------
 * Function   | Retrieves the molar mass of the alkane atom based on the given CarbonType.
 * Parameters | - Type: The CarbonType enum representing the type of alkane pseudoatom.
 * Returns    | The molar mass of the alkane atom for the specified alkane type.
 * **************************************************************************************************************/
double GetAlkaneAtomMolarMass(enum CarbonType Type){
  double MolarMass = 0.0;
  switch (Type){
    case CH4:
      MolarMass = 16.04206;
      break;
    case CH3e:
      MolarMass = 15.03422;
      break;
    case CH3:
      MolarMass = 15.03422;
      break;
    case CH2:
      MolarMass = 14.02638;
      break;

  }
  return MolarMass;
}

/* ***************************************************************************************************************
 * Name       | CalculateBendingAngle
 * ---------------------------------------------------------------------------------------------------------------
 * Function   | Calculates the bending angle between three vectors.
 * Parameters |
 *              - r1: The first vector.
 *              - r2: The second vector (center vector).
 *              - r3: The third vector.
 * Returns    | The bending angle in radians.
 * **************************************************************************************************************/
double CalculateBendingAngle(VECTOR r1, VECTOR r2, VECTOR r3){
  VECTOR r21, r23;
  r21 = VectorSubtraction(r1, r2);
  r23 = VectorSubtraction(r3, r2);
  return InternalAngle(r21, r23);
}

/* ***************************************************************************************************************
 * Name       | CalculateTorsionAngle
 * ---------------------------------------------------------------------------------------------------------------
 * Function   | Calculates the torsion angle between four vectors.
 * Parameters |
 *              - r1: The first vector.
 *              - r2: The second vector.
 *              - r3: The third vector.
 *              - r4: The fourth vector.
 * Returns    | The torsion angle in radians.
 * **************************************************************************************************************/
double CalculateTorsionAngle(VECTOR r1, VECTOR r2, VECTOR r3, VECTOR r4){
  VECTOR d12, d23, d34;
  VECTOR d23xd12, d34xd23;
  VECTOR r21, r23, r34;
  VECTOR r, s;
  VECTOR DihedralVector1, DihedralVector2;
  double CosPhi, CosPhi2, SinPhi, phi, r21r23, r34r23;

  d12 = VectorSubtraction(r2, r1); 
  d23 = VectorSubtraction(r3, r2); 
  d34 = VectorSubtraction(r4, r3); 

  d23xd12 = CrossProduct(d23, d12); 
  d23xd12 = NormalizeVector(d23xd12);
  d34xd23 = CrossProduct(d34, d23); 
  d34xd23 = NormalizeVector(d34xd23);

  CosPhi = DotProduct(d34xd23, d23xd12);

  phi = acos(CosPhi);
  
  return phi;
}