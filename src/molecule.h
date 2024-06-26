/* ***************************************************************************************************************
 * This code was written as a part of LAMC, a program to perform molecular Monte Carlo simulations of linear
 * alkanes on multiple ensembles.
 * 
 * This module molecule.h is a header file for defining molecular
 * 
 * Developer: Rodolfo J Amancio (rodolfojamancio@gmail.com)
 * University of Campinas, School of Chemical Engineering, Campinas/SP -	Brazil
 * 
 * **************************************************************************************************************/
#ifndef MOLECULE_H
#define MOLECULE_H
#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<stdbool.h>
#include<math.h>
#include"constants.h"
#include"mathematics.h"
#include"vectors.h"
#include"simulation_setup.h"


// ---------------------------------------------------------------------------------------------------------------
// Macros
// ---------------------------------------------------------------------------------------------------------------

#define NUMBER_PSEUDO_ATOMS_TYPES 4

// ---------------------------------------------------------------------------------------------------------------
// Force-field parameters
// ---------------------------------------------------------------------------------------------------------------

static const double CUTOFF_DISTANCE  = 14; // A
static const double CUTOFF_SQUARED = (Squared(CUTOFF_DISTANCE));
// Intramolecular interacions
// bond streaching (TraPPE)
static const double dEq = 1.529, kStreaching = 96500*BOLTZMANN_CONSTANT; // in A and J/A²
// bond bending (TraPPE)
static const double thetaEq = 114*M_PI/180., kBending = 62500*BOLTZMANN_CONSTANT; // in rad and J/rad²
// molecular torsion (TraPPE)
static const double cTorsion[3] = {
    335.03*BOLTZMANN_CONSTANT, 
    -68.19*BOLTZMANN_CONSTANT, 
    791.32*BOLTZMANN_CONSTANT
};  // in J

// ---------------------------------------------------------------------------------------------------------------
// Structures definitions
// ---------------------------------------------------------------------------------------------------------------

enum CarbonType {CH4, CH3e, CH3, CH2};

typedef struct{
    enum CarbonType Type;
    VECTOR Position;
    VECTOR Force;
    double Sigma;
    double Epsilon;
    double RepulsiveExponent;
    double AttractiveExponent;
    double MolarMass;
} ATOM;

typedef struct{
    ATOM* Atoms;
    int Size;
    double MolarMass;
    VECTOR CenterOfMass;
    double OrderParameter;
} MOLECULE;

typedef struct{
	int 	  NumberMolecules;
	MOLECULE* Molecules;
} CONFIGURATION;

// ---------------------------------------------------------------------------------------------------------------
// Molecules functions
// ---------------------------------------------------------------------------------------------------------------

// Combinig rules fucntions
double GetSigmaCombination(double SigmaA, double SigmaB);
double GetEpsilonCombination(double EpsilonA, double EpsilonB);
double GetExponentCombination(double ExponentA, double ExponentB);

// Functions for configuring molecules
enum CarbonType GetCarbonType(int ChainLength, int Position);
double GetAlkaneRepulsiveExponent(enum CarbonType Type);
double GetAlkaneAttractiveExponent(enum CarbonType Type);
double GetAlkaneEpsilon(enum CarbonType Type);
double GetAlkaneSigma(enum CarbonType Type);
double GetAlkaneAtomMolarMass(enum CarbonType Type);

// Molecules conformation functions
double CalculateBendingAngle(VECTOR r1, VECTOR r2, VECTOR r3);
double CalculateTorsionAngle(VECTOR r1, VECTOR r2, VECTOR r3, VECTOR r4);

#endif