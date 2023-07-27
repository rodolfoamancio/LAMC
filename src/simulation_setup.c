/* ***************************************************************************************************************
 * This code was written as a part of LAMC, a program to perform molecular Monte Carlo simulations of linear
 * alkanes on multiple ensembles.
 * 
 * This module properties.c is a code file for defining auxliary global variables
 * 
 * Developer: Rodolfo J Amancio (rodolfojamancio@gmail.com)
 * University of Campinas, School of Chemical Engineering, Campinas/SP -	Brazil
 * 
 * **************************************************************************************************************/

#include"simulation_setup.h"

double  Temperature;
double  SimulationTemperature;
double  SimulationPressure;
double  SimulationFugacity;
double  Beta;
double  DisplacementAcceptanceTarget;
double  VolumeAcceptanceTarget;
double  WeightIdealChain;
double* ProbabilityCavitySummed;
double* ProbabilityCavityArray;
double  AveragePressure;
double  AveragePressureExcess;
double  AveragePressureIdealGas;
double  AveragePressureLongRangeCorrection;
double  AveragePotential;
double  AveragePotentialNonbonded;
double  AveragePotentialBonded;
double  AveragePotentialLongRangeCorrection;
double  AveragePotentialWalls;
double  AverageWeightGhostMolecule;
double  AverageNumberMolecules;
double  AverageDensityMass;
double  AverageDensityMolar;
double  MaxTranslationDistance;
double  MaxLnVolumeChange;
double  DisplacementAcceptanceCalculated;
double  VolumeAcceptanceCalculated;
double  DeletionAttemptProbability;
double  InsertionAttemptProbability;
double  DisplacementAttemptProbability;
bool    MovementAccepted;
bool    TemperatureRamp;
bool    Debug;
time_t  StartTime;
time_t  EndTime;
int     DisplacementStepsAccepted;
int     DisplacementStepsAttempted;
int     InsertionStepsAccepted;
int     InsertionStepsAttempted;
int     DeletionStepsAccepted;
int     DeletionStepsAttempted;
int     VolumeStepsAccepted;
int     VolumeStepsAttempted;
int     InitialNumberMolecules;
int     NumberTotalSteps;
int     NumberEquilibrationSteps;
int     NumberTotalCycles;
int     NumberEquilibrationCycles;
int     StepsCalculateProperties;
int     StepsRecordConfiguration;
int     StepsGetProfiles;
int     CyclesCalculateProperties;
int     CyclesRecordConfiguration;
int     CyclesGetProfiles;
int     NumberBinsForProfiles;
int     ChainSize;
int     CountProductionStates;
int*    NumberCavityScan;