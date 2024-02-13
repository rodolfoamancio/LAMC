/* ***************************************************************************************************************
 * This code was written as a part of LAMC, a program to perform molecular Monte Carlo simulations of linear
 * alkanes on multiple ensembles.
 * 
 * This module monte_carlo.h is a header file for Monte Carlo procedures
 * 
 * Developer: Rodolfo J Amancio (rodolfojamancio@gmail.com)
 * University of Campinas, School of Chemical Engineering, Campinas/SP -  Brazil
 * 
 * **************************************************************************************************************/

#ifndef MONTE_CARLO_H
#define MONTE_CARLO_H
#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<stdbool.h>
#include<math.h>
#include"mathematics.h"
#include"vectors.h"
#include"molecule.h"
#include"potential.h"
#include"structural.h"
#include"simulation_setup.h"

#define NUMBER_TRIAL_ORIENTATIONS 30

enum ensemble {NVT, NPT, muVT};
extern enum ensemble SimulationEnsemble;
extern enum ensemble EquilibraionEnsemble;
extern enum ensemble Ensemble;

// ---------------------------------------------------------------------------------------------------------------
// Configurational bias monte carlo algorithm
// ---------------------------------------------------------------------------------------------------------------

// Memory handling functions
void CopyConfiguration(CONFIGURATION Source, CONFIGURATION* Dest);
void InitializeConfiguration(CONFIGURATION* Configuration);

// Initialaization function
double GenerateInitialConfiguration(CONFIGURATION* Configuration);

// Supplementary function
double GetRosenbluthWeightIdealChain(MOLECULE Molecule);
double GetRosenbluthWeightGhostMolecule(CONFIGURATION Configuration);

// Chain grown function
VECTOR SampleBeadPosition(MOLECULE Molecule, int Bead);
int SelectTrialOrientation(double *WeightList, double TotalWeight);

// Displacement function
void GetMoleculeDisplacement(CONFIGURATION OldConfiguration, CONFIGURATION* NewConfiguration);

// Volume change functions
void GetVolumeChange(CONFIGURATION OldConfiguration, CONFIGURATION* NewConfiguration);

// Moleculeculs inserion/deletion
void GetMoleculeInsertion(CONFIGURATION OldConfiguration, CONFIGURATION* NewConfiguration);
void GetMoleculeDeletion(CONFIGURATION OldConfiguration, CONFIGURATION* NewConfiguration);

// Ensemble functions
void GetNVTMovement(CONFIGURATION OldConfiguration, CONFIGURATION* NewConfiguration);
void GetNPTMovement(CONFIGURATION OldConfiguration, CONFIGURATION* NewConfiguration);
void GetmuVTMovement(CONFIGURATION OldConfiguration, CONFIGURATION* NewConfiguration);

// Final function
void GetMonteCarloMove(CONFIGURATION OldConfiguration, CONFIGURATION* NewConfiguration);

#endif

// ---------------------------------------------------------------------------------------------------------------