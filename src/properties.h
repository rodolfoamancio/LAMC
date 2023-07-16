/* ***************************************************************************************************************
 * This code was written as a part of LAMC, a program to perform molecular Monte Carlo simulations of linear
 * alkanes on multiple ensembles.
 * 
 * This module properties.h is a header file for calculating properties
 * 
 * Developer: Rodolfo J Amancio (rodolfojamancio@gmail.com)
 * University of Campinas, School of Chemical Engineering, Campinas/SP -	Brazil
 * 
 * **************************************************************************************************************/

#ifndef PROPERTIES_H
#define PROPERTIES_H
#include<stdio.h>
#include<stdbool.h>
#include"constants.h"
#include"molecule.h"
#include"vectors.h"
#include"simulation_setup.h"
#include"structural.h"
#include"potential.h"
#include"simulation_setup.h"

double GetPressureLongRangeCorrection(CONFIGURATION Configuration);
double GetPressureIdealGas(int numberOfMolecules, double volume);
double GetPressureExcess(CONFIGURATION Configuration);

double GetDensityMolar(CONFIGURATION Configuration);
double GetDensityMass(CONFIGURATION Configuration);

#endif