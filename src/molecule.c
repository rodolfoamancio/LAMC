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
 * Name       | GetSigma
 * ---------------------------------------------------------------------------------------------------------------
 * Function   | Calculates the average of two given sigma values.
 * Parameters |
 *              - SigmaA: The first sigma value.
 *              - SigmaB: The second sigma value.
 * Returns    | The average of SigmaA and SigmaB.
 * **************************************************************************************************************/
double GetSigma(double SigmaA, double SigmaB){
  return (SigmaA != SigmaB) ? (SigmaA + SigmaB)/2 : SigmaA;
}

/* ***************************************************************************************************************
 * Name       | GetEpsilon
 * ---------------------------------------------------------------------------------------------------------------
 * Function   | Calculates the square root of the product of two given epsilon values.
 * Parameters |
 *              - EpsilonA: The first epsilon value.
 *              - EpsilonB: The second epsilon value.
 * Returns    | The square root of the product of EpsilonA and EpsilonB.
 * **************************************************************************************************************/
double GetEpsilon(double EpsilonA, double EpsilonB){
  return (EpsilonA != EpsilonB) ? sqrt(EpsilonA*EpsilonB) : EpsilonA;
}

/* ***************************************************************************************************************
 * Name       | GetInteractionExponent
 * ---------------------------------------------------------------------------------------------------------------
 * Function   | Applies the combination rule from Lafitte (2013) for the exponents of the Mie potential
 * Parameters |
 *              - ExponentA: The first exponent value.
 *              - ExponentA: The second exponent value.
 * Returns    | The combination exponent
 * **************************************************************************************************************/
double GetInteractionExponent(double ExponentA, double ExponentB){
  return (ExponentA != ExponentB) ? sqrt((ExponentA - 3)*(ExponentB - 3)) + 3 : ExponentA;
}

/* ***************************************************************************************************************
 * Name       | GetAlkaneEpsilon
 * ---------------------------------------------------------------------------------------------------------------
 * Function   | Retrieves the epsilon value based on the given CarbonType.
 * Parameters | - Type: The CarbonType enum representing the type of alkane pseudoatom.
 * Returns    | The epsilon value for the specified alkane type.
 * **************************************************************************************************************/
double GetAlkaneEpsilon(enum CarbonType Type){
  double Epsilon = 0.0;
  switch(Type){
    case CH4:
      Epsilon = EpsilonAlkane[0];
      break;
    case CH3:
      Epsilon = EpsilonAlkane[1];
      break;
    case CH2:
      Epsilon = EpsilonAlkane[2];
  }
  return Epsilon;
}

/* ***************************************************************************************************************
 * Name       | GetAlkaneSigma
 * ---------------------------------------------------------------------------------------------------------------
 * Function   | Retrieves the sigma value based on the given CarbonType.
 * Parameters | - Type: The CarbonType enum representing the type of alkane pseudoatom.
 * Returns    | The sigma value for the specified alkane type.
 * **************************************************************************************************************/
double GetAlkaneSigma(enum CarbonType Type){
  double Sigma = 0.0;
  switch (Type){
    case CH4:
      Sigma = SigmaAlkane[0];
      break;
    case CH3:
      Sigma = SigmaAlkane[1];
      break;
    case CH2:
      Sigma = SigmaAlkane[2];
      break;

  }
  return Sigma;
}

/* ***************************************************************************************************************
 * Name       | GetAlkaneRepulsiveExponent
 * ---------------------------------------------------------------------------------------------------------------
 * Function   | Retrieves the repulsive exponent value based on the given CarbonType.
 * Parameters | - Type: The CarbonType enum representing the type of alkane pseudoatom.
 * Returns    | The repulsive exponent value for the specified alkane type.
 * **************************************************************************************************************/
double GetAlkaneRepulsiveExponent(enum CarbonType Type){
  double exponent = 12;
  switch (Type){
  case CH4:
    exponent = RepulsiveExponentAlkane[0];
    break;
  case CH3:
    exponent = RepulsiveExponentAlkane[1];
    break;
  case CH2:
    exponent = RepulsiveExponentAlkane[1];
    break;
  }
  return exponent;
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
      MolarMass = PseudoAtomMolarMass[0];
      break;
    case CH3:
      MolarMass = PseudoAtomMolarMass[1];
      break;
    case CH2:
      MolarMass = PseudoAtomMolarMass[2];
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