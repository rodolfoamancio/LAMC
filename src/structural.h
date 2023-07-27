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

VECTOR GetMoleculeCenterOfMass(MOLECULE Molecule);
void GetCenterOfMassAllMolecules(CONFIGURATION *Configuration);
void CalculateOrderParameter(CONFIGURATION *Configuraiton);

#endif