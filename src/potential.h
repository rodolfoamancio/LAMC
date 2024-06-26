/* ***************************************************************************************************************
 * This code was written as a part of LAMC, a program to perform molecular Monte Carlo simulations of linear
 * alkanes on multiple ensembles.
 * 
 * This module potential.h is a header file for defining potential functions and related resources
 * 
 * Developer: Rodolfo J Amancio (rodolfojamancio@gmail.com)
 * University of Campinas, School of Chemical Engineering, Campinas/SP -	Brazil
 * 
 * **************************************************************************************************************/
#ifndef POTENTIAL_H
#define POTENTIAL_H

#include<stdio.h>
#include<stdbool.h>
#include<math.h>
#include"constants.h"
#include"molecule.h"
#include"simulation_setup.h"
#include"vectors.h"
#include"simulation_setup.h"
#include"structural.h"

typedef struct{
	double potential;
	bool   overlap;
} POTENTIAL;

enum PotentialType {
  MIE, 
  HARD_SPHERE, 
  BARKER_HENDERSON_REFERENCE, 
  BARKER_HENDERSON_PERTURBED, 
  WEEKS_CHANDLER_ANDERSEN_REFERENCE, 
  WEEKS_CHANDLER_ANDERSEN_PERTURBED
};
extern enum PotentialType ReferencePotential, PerturbationPotential;

double GetCMie(double RepulsiveExponent, double AttractiveExponent);
double GetPotentialMie(
  double RepulsiveExponent, 
  double AttractiveExponent, 
  double Sigma, 
  double Epsilon, 
  double Distance);
double GetBHPotentialReference(
  double RepulsiveExponent, 
  double AttractiveExponent, 
  double Sigma, 
  double Epsilon, 
  double Distance);
double GetBHPotentialPerturbed(
  double RepulsiveExponent, 
  double AttractiveExponent, 
  double Sigma, 
  double Epsilon, 
  double Distance);
double GetWCAPotentialReference(
  double RepulsiveExponent, 
  double AttractiveExponent, 
  double Sigma, 
  double Epsilon, 
  double Distance);
double GetWCAPotentialPerturbed(
  double RepulsiveExponent, 
  double AttractiveExponent, 
  double Sigma, 
  double Epsilon, 
  double Distance);
double GetNonbondedPotentialPair(
  double RepulsiveExponent, 
  double AttractiveExponent, 
  double Sigma, 
  double Epsilon, 
  double Distance,
  enum PotentialType Potential);
double GetPotentialStreaching(double d);
double GetPotentialBending(double theta);
double GetPotentialTorsion(double phi);
double GetPotentialBonded(CONFIGURATION Configuration);

double SampleBondLength(enum CarbonType TypeA, enum CarbonType TypeB);
VECTOR SampleBendingAngle(VECTOR r1, VECTOR r2);
VECTOR SampleBendingTorsionAngles(VECTOR r1, VECTOR r2, VECTOR r3);

POTENTIAL GetPartialExternalPotential(CONFIGURATION Configuration, int referenceMolecule, int referenceParticle);
POTENTIAL GetPotentialExternalField(ATOM Atom);
VECTOR GetPotentialGradient(VECTOR SeparationVector, ATOM AtomA, ATOM AtomB);
double GetPotentialSteele(ATOM Atom, double height);
double GetPotentialNonbonded(CONFIGURATION Configuration, enum PotentialType Potential);
double GetTotalPotentialExternal(CONFIGURATION Configuration);
double GetPotentialLongRangeCorrection(CONFIGURATION Configuration, enum PotentialType Potential);
void ComputeStrainDerivativeTensor(CONFIGURATION *Configuration);


#endif