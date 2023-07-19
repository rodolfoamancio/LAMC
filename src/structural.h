/* ***************************************************************************************************************
 * This code was written as a part of LAMC, a program to perform molecular Monte Carlo simulations of linear
 * alkanes on multiple ensembles.
 * 
 * This module structural.h is a header file for with struc and functions for structutal quantities
 * 
 * Developer: Rodolfo J Amancio (rodolfojamancio@gmail.com)
 * University of Campinas, School of Chemical Engineering, Campinas/SP -  Brazil
 * 
 * **************************************************************************************************************/
#ifndef STRUCTURAL_H
#define STRUCTURAL_H
#include<stdio.h>
#include<stdbool.h>
#include"constants.h"
#include"mathematics.h"
#include"vectors.h"
#include"molecule.h"
#include"simulation_setup.h"

typedef struct{
  int numberOfBins;
  double binHeight;
  double binsDelimiters[MAX_PROFILE_SIZE];
  double positions[MAX_PROFILE_SIZE];
  double values[MAX_PROFILE_SIZE];
} PROFILE;

VECTOR GetMoleculeCenterOfMass(MOLECULE Molecule);
void GetCenterOfMassAllMolecules(CONFIGURATION *Configuration);
void CalculateOrderParameter(CONFIGURATION *Configuraiton);
void DetermineDensityProfile(CONFIGURATION Configuration, PROFILE *DensityProfile);
void DetermineOrientationProfile(CONFIGURATION Configuration, PROFILE *OrientationProfile);
PROFILE ConfigureProfile(int numberOfBins);
PROFILE SumProfiles(PROFILE Profile1, PROFILE Profile2);

#endif