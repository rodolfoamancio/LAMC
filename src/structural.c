/* ***************************************************************************************************************
 * This code was written as a part of LAMC, a program to perform molecular Monte Carlo simulations of linear
 * alkanes on multiple ensembles.
 * 
 * This module structural.c is a header file for with struc and functions for structutal quantities
 * 
 * Developer: Rodolfo J Amancio (rodolfojamancio@gmail.com)
 * University of Campinas, School of Chemical Engineering, Campinas/SP -  Brazil
 * 
 * **************************************************************************************************************/
#include"structural.h"

/***************************************************************
 * Name       | GetMoleculeCenterOfMass
 * -------------------------------------------------------------
 * Function   | Calculates the center of mass position of the 
 *              molecule
 * Parameters | Molecule: a MOLECULE struct type variable
 * Returns    | The VECTOR struct witht the position of the 
 *              center of mass
 ****************************************************************/
VECTOR GetMoleculeCenterOfMass(MOLECULE Molecule) {
  VECTOR CenterOfMass;
  double totalMass = 0.0; // Track the total mass of the molecule

  CenterOfMass.x = CenterOfMass.y = CenterOfMass.z = 0;

  for (int j = 0; j < Molecule.Size; j++) {
    CenterOfMass.x += Molecule.Atoms[j].Position.x * Molecule.Atoms[j].MolarMass;
    CenterOfMass.y += Molecule.Atoms[j].Position.y * Molecule.Atoms[j].MolarMass;
    CenterOfMass.z += Molecule.Atoms[j].Position.z * Molecule.Atoms[j].MolarMass;
    totalMass += Molecule.Atoms[j].MolarMass;
  }

  if (totalMass > 0.0) {
    CenterOfMass.x /= totalMass;
    CenterOfMass.y /= totalMass;
    CenterOfMass.z /= totalMass;
  }

  return CenterOfMass;
}

/***************************************************************
 * Name       | GetCenterOfMassAllMolecules
 * -------------------------------------------------------------
 * Function   | Calculates the center of mass position of all 
 *              molecules in a given configuration
 * Parameters | Configuration: a pointer to a CONFIGURATION 
 *              struct type variable
 * Returns    | None
 ****************************************************************/
void GetCenterOfMassAllMolecules(CONFIGURATION *Configuration){
  for(int i = 0; i < Configuration->NumberMolecules; i++){
    Configuration->Molecules[i].CenterOfMass = GetMoleculeCenterOfMass(Configuration->Molecules[i]);
  }
}


/***************************************************************
 * Name       | CalculateOrderParameter
 * -------------------------------------------------------------
 * Function   | Calculates the order parameter for each molecule 
 *              in a given configuration
 * Parameters | Configuration: a pointer to a CONFIGURATION struct
 * Returns    | None
 ****************************************************************/
void CalculateOrderParameter(CONFIGURATION *Configuraiton){
  VECTOR HeadTailVector, zVersor;
  double OrientationAngle;
  zVersor.x = zVersor.y = 0;
  zVersor.z = 1;
  for(int i = 0; i < Configuraiton->NumberMolecules; i++){
    if(Configuraiton->Molecules[i].Size > 1){
      HeadTailVector = VectorSubtraction(
        Configuraiton->Molecules[i].Atoms[Configuraiton->Molecules[i].Size - 1].Position, 
        Configuraiton->Molecules[i].Atoms[0].Position
      );
      OrientationAngle = InternalAngle(HeadTailVector, zVersor);
      Configuraiton->Molecules[i].OrderParameter = 0.5*(3.*Squared(cos(OrientationAngle)) - 1);
    }else{
      Configuraiton->Molecules[i].OrderParameter = 0;
    }
    
  }
}

/***************************************************************
 * Name       | GetDensityMolar
 * -------------------------------------------------------------
 * Function   | Calculates and returns the molar density based on 
 *              the given configuration
 * Parameters | Configuration: the configuration of the system
 * Returns    | The molar density
 ****************************************************************/
double GetDensityMolar(CONFIGURATION Configuration){
  return Configuration.NumberMolecules/(SimulationBox.volume*Cube(ANGSTRON)*AVOGADRO_NUMBER);
}

/***************************************************************
 * Name       | GetDensityMass
 * -------------------------------------------------------------
 * Function   | Calculates and returns the mass density based on 
 *              the given configuration
 * Parameters | Configuration: the configuration of the system
 * Returns    | The mass density
 ****************************************************************/
double GetDensityMass(CONFIGURATION Configuration){
  double TotalMass = 0;
  for(int i = 0; i < Configuration.NumberMolecules; i++){
    TotalMass += Configuration.Molecules[i].MolarMass;
  }
  return (TotalMass/AVOGADRO_NUMBER)/(1000*SimulationBox.volume*Cube(ANGSTRON));
}