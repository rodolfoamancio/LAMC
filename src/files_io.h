/* ***************************************************************************************************************
 * This code was written as a part of LAMC, a program to perform molecular Monte Carlo simulations of linear
 * alkanes on multiple ensembles.
 * 
 * This module files_io.h is a header file for defining files input output functions
 * 
 * Developer: Rodolfo J Amancio (rodolfojamancio@gmail.com)
 * University of Campinas, School of Chemical Engineering, Campinas/SP -	Brazil
 * 
 * **************************************************************************************************************/

#ifndef PARSERS_H
#define PARSERS_H
#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<stdbool.h>
#include<math.h>
#include<time.h>
#include"mathematics.h"
#include"potential.h"
#include"monte_carlo.h"
#include"vectors.h"
#include"simulation_setup.h"

// ---------------------------------------------------------------------------------------------------------------
// File parsing functions
// ---------------------------------------------------------------------------------------------------------------

void RemoveSubstring (char *string, char *sub);

void ReadInitialConfiguration(char InitialConfigurationFilePath[], CONFIGURATION* Configuration);
void ReadInputFile(char inputsFilePath[]);

FILE* InitializePropertiesDataFile(char BaseName[]);

void  RecordInitialConfiguraiton(char BaseName[], CONFIGURATION Configuration);
FILE* InitializeTrajectoryFile(char BaseName[]);
void  RecordConfiguration(FILE* IntermediateConfigurationFile, CONFIGURATION Configuration);

FILE* InitializeProfilesFile(char BaseName[]);
void  RecordProfileData(FILE* ProfileFile, CONFIGURATION Configuration);

void  RecordSimulationLog(char BaseName[]);

char* GetEnsembleLabel(enum ensemble InputEnsemble);
char* GetPotentialTypeLabel(enum PotentialType Potential);

#endif

