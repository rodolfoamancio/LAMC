/* ***************************************************************************************************************
 * This code was written as a part of LAMC, a program to perform molecular Monte Carlo simulations of linear
 * alkanes on multiple ensembles.
 * 
 * This module simulation_setup.h is a header file for defininf global variables
 * 
 * Developer: Rodolfo J Amancio (rodolfojamancio@gmail.com)
 * University of Campinas, School of Chemical Engineering, Campinas/SP -	Brazil
 * 
 * **************************************************************************************************************/

#ifndef SIMULATION_SETUP_H
#define SIMULATION_SETUP_H
#include<stdio.h>
#include<stdbool.h>
#include<math.h>
#include<time.h>
#include"constants.h"

extern double  Temperature;
extern double  SimulationTemperature;
extern double  SimulationPressure;
extern double  SimulationFugacity;
extern double  Beta;

extern double  DisplacementAcceptanceTarget;
extern double  VolumeAcceptanceTarget;
extern double  WeightIdealChain;
extern double* ProbabilityCavitySummed;
extern double* ProbabilityCavityArray;

extern double  AveragePressure;
extern double  AveragePressureExcess;
extern double  AveragePressureIdealGas;
extern double  AveragePressureLongRangeCorrection;

extern double  AveragePotential;
extern double  AveragePotentialNonbonded;
extern double  AveragePotentialBonded;
extern double  AveragePotentialLongRangeCorrection;
extern double  AveragePotentialWalls;

extern double  AverageWeightGhostMolecule;

extern double  AverageNumberMolecules;
extern double  AverageDensityMass;
extern double  AverageDensityMolar;

extern double  BlockPressure[MAX_NUMBER_BLOCKS];
extern double  BlockPressureExcess[MAX_NUMBER_BLOCKS];
extern double  BlockPressureIdealGas[MAX_NUMBER_BLOCKS];
extern double  BlockPressureLongRangeCorrection[MAX_NUMBER_BLOCKS];

extern double  BlockPotential[MAX_NUMBER_BLOCKS];
extern double  BlockPotentialNonbonded[MAX_NUMBER_BLOCKS];
extern double  BlockPotentialBonded[MAX_NUMBER_BLOCKS];
extern double  BlockPotentialLongRangeCorrection[MAX_NUMBER_BLOCKS];
extern double  BlockPotentialWalls[MAX_NUMBER_BLOCKS];

extern double  BlockWeightGhostMolecule[MAX_NUMBER_BLOCKS];

extern double  BlockNumberMolecules[MAX_NUMBER_BLOCKS];
extern double  BlockDensityMass[MAX_NUMBER_BLOCKS];
extern double  BlockDensityMolar[MAX_NUMBER_BLOCKS];

extern double  MaxTranslationDistance;
extern double  MaxLnVolumeChange;
extern double  DisplacementAcceptanceCalculated;
extern double  VolumeAcceptanceCalculated;
extern double  DeletionAttemptProbability;
extern double  InsertionAttemptProbability;
extern double  DisplacementAttemptProbability;

extern bool    MovementAccepted;
extern bool    TemperatureRamp;
extern bool    Debug;
extern time_t  StartTime;
extern time_t  EndTime;
extern int     DisplacementStepsAccepted;
extern int     DisplacementStepsAttempted;
extern int     InsertionStepsAccepted;
extern int     InsertionStepsAttempted;
extern int     DeletionStepsAccepted;
extern int     DeletionStepsAttempted;
extern int     VolumeStepsAccepted;
extern int     VolumeStepsAttempted;
extern int     InitialNumberMolecules;
extern int     NumberTotalSteps;
extern int     NumberEquilibrationSteps;
extern int     NumberTotalCycles;
extern int     NumberEquilibrationCycles;
extern int     StepsCalculateProperties;
extern int     StepsRecordConfiguration;
extern int     StepsGetProfiles;
extern int     CyclesCalculateProperties;
extern int     CyclesRecordConfiguration;
extern int     CyclesGetProfiles;
extern int     NumberBinsForProfiles;
extern int     ChainSize;
extern int     CountProductionStates;
extern int*    NumberCavityScan;

#endif