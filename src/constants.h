/* ***************************************************************************************************************
 * This code was written as a part of LAMC, a program to perform molecular Monte Carlo simulations of linear
 * alkanes on multiple ensembles.
 * 
 * This module constants.h is a header file with contants used durint the simulation
 * 
 * Developer: Rodolfo J Amancio (rodolfojamancio@gmail.com)
 * University of Campinas, School of Chemical Engineering, Campinas/SP -	Brazil
 * 
 * **************************************************************************************************************/

#ifndef CONSTANTS_H
#define CONSTANTS_H
#include<stdio.h>
#include<math.h>

// ---------------------------------------------------------------------------------------------------------------
// Constants definitions
// ---------------------------------------------------------------------------------------------------------------

// Physical constants
#define BOLTZMANN_CONSTANT 1.38064852E-23 // J/K 
#define AVOGADRO_NUMBER 6.02214086E23 // 1/mol
#define PLANCK_CONSTANT 6.62607015E-34 // J.s
#define GAS_CONSTANT 8.314 // J/mol.K

// Conversion constants
#define ANGSTRON 1.0E-10 // m
#define PICO_SECOND 1.0E-12 // s
#define ATOM_MASS 1.6605402E-27 // kg
#define METER_TO_ANGSTRON 1E10 // A/m 
#define KILOGRAM_TO_GRAM 1E3 // g/kg

// Arrays allocation limits and auxiliary constants
#define MAX_NUMBER_BLOCKS 10
#define MAX_NUMBER_MOLECULES 2000
#define MAX_CHAIN_SIZE 10
#define MAX_PROFILE_SIZE 200
#define NUMBER_SAMPLES_IDEAL_CHAIN 50000

#endif
