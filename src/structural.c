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
 * Name       | GetMoleculeOrderParameter
 * -------------------------------------------------------------
 * Function   | Calculates the order parameter of a molecule based 
 *              on its size and orientation
 * Parameters | Molecule: a MOLECULE struct type variable
 * Returns    | double: the calculated order parameter value
 ****************************************************************/
double GetMoleculeOrderParameter(MOLECULE Molecule){
  VECTOR HeadTailVector, zVersor;
  double Angle, OrderParameter;
  zVersor.x=zVersor.y=0.0;
  zVersor.z=1.0;
  if(Molecule.Size>1){
    HeadTailVector = VectorSubtraction(
      Molecule.Atoms[0].Position,
      Molecule.Atoms[Molecule.Size-1].Position
    );
    Angle = InternalAngle(HeadTailVector, zVersor);
    OrderParameter = 0.5*(3.0*Squared(cos(Angle)) - 1.0);
  }else{
    OrderParameter = 0.0;
  }
  return OrderParameter;
}

/***************************************************************
 * Name       | GetOrderParameterAllMolecules
 * -------------------------------------------------------------
 * Function   | Calculates the order parameter for all molecules 
 *              in a given configuration
 * Parameters | Configuration: a pointer to a CONFIGURATION struct
 * Returns    | None
 ****************************************************************/
void GetOrderParameterAllMolecules(CONFIGURATION *Configuration){
  for(int i=0; i<Configuration->NumberMolecules; i++){
    Configuration->Molecules[i].OrderParameter = GetMoleculeOrderParameter(Configuration->Molecules[i]);
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
    HeadTailVector = VectorSubtraction(Configuraiton->Molecules[i].Atoms[Configuraiton->Molecules[i].Size - 1].Position, Configuraiton->Molecules[i].Atoms[0].Position);
    OrientationAngle = InternalAngle(HeadTailVector, zVersor);
    Configuraiton->Molecules[i].OrderParameter = 0.5*(3.*Squared(cos(OrientationAngle)) - 1);
  }
}

/***************************************************************
 * Name       | ConfigureProfile
 * -------------------------------------------------------------
 * Function   | Configures a profile with the specified number of 
 *              bins and initializes its values
 * Parameters | numberOfBins: the number of bins in the profile
 * Returns    | PROFILE: the configured profile
 ****************************************************************/
PROFILE ConfigureProfile(int numberOfBins){
  PROFILE Profile;
  Profile.numberOfBins = numberOfBins;
  Profile.binHeight = (SimulationBox.zSize)/numberOfBins;
  Profile.binsDelimiters[0] = -SimulationBox.zSize/2;
  for(int i = 1; i < numberOfBins + 1; i++){
    Profile.binsDelimiters[i] = Profile.binsDelimiters[i-1] + Profile.binHeight;
    Profile.positions[i-1] = (Profile.binsDelimiters[i] + Profile.binsDelimiters[i-1])/2;
    Profile.values[i-1] = 0;
  }
  return Profile;
}

/***************************************************************
 * Name       | SumProfiles
 * -------------------------------------------------------------
 * Function   | Sums two profiles with matching dimensions and 
 *              returns the resulting profile
 * Parameters | Profile1: the first profile
 *             Profile2: the second profile
 * Returns    | PROFILE: the summed profile
 ****************************************************************/
PROFILE SumProfiles(PROFILE Profile1, PROFILE Profile2){
  PROFILE Profile = Profile1;
  if(Profile1.numberOfBins == Profile2.numberOfBins){
    for(int i = 0; i < Profile1.numberOfBins; i++){
      Profile.values[i] = Profile1.values[i] + Profile2.values[i];
    }
  }else{
    printf("ERROR: Profiles' dimensions do not match!\n");
  }
  return Profile;
}

/***************************************************************
 * Name       | DetermineDensityProfile
 * -------------------------------------------------------------
 * Function   | Determines the density profile based on a given 
 *              configuration and updates the provided DensityProfile
 * Parameters | Configuration: the configuration of the system
 *             DensityProfile: a pointer to the density profile
 * Returns    | (void)
 ****************************************************************/
void DetermineDensityProfile(CONFIGURATION Configuration, PROFILE *DensityProfile){
  double centerOfMassZPosition[MAX_NUMBER_MOLECULES] = {0};
  double binVolume = SimulationBox.xSize*SimulationBox.ySize*DensityProfile->binHeight/1E3;
  double AuxMass;

  for(int i = 0; i < (DensityProfile->numberOfBins); i++){
    int CountMoleculesInssideBin = 0;
    for(int j = 0; j < Configuration.NumberMolecules; j++){
      if(Configuration.Molecules[j].CenterOfMass.z > DensityProfile->binsDelimiters[i] && Configuration.Molecules[j].CenterOfMass.z <= DensityProfile->binsDelimiters[i+1]){
        CountMoleculesInssideBin++;
      }
    }
    DensityProfile->values[i] = CountMoleculesInssideBin/binVolume;
  }
}

/***************************************************************
 * Name       | DetermineOrientationProfile
 * -------------------------------------------------------------
 * Function   | Determines the orientation profile based on a 
 *              given configuration and updates the provided 
 *              OrientationProfile
 * Parameters | Configuration: the configuration of the system
 *             OrientationProfile: a pointer to the orientation
 *             profile
 * Returns    | (void)
 ****************************************************************/
void DetermineOrientationProfile(CONFIGURATION Configuration, PROFILE *OrientationProfile){
  if(Configuration.Molecules[0].Size > 1){
    for(int i = 0; i < OrientationProfile->numberOfBins; i++){
      int countMoleculesInsideBin = 0;
      double OrderParameterSum = 0;
      for(int j = 0; j < Configuration.NumberMolecules; j++){
        if(Configuration.Molecules[j].CenterOfMass.z > OrientationProfile->binsDelimiters[i] && Configuration.Molecules[j].CenterOfMass.z <= OrientationProfile->binsDelimiters[i+1]){
          countMoleculesInsideBin++;
          OrderParameterSum += Configuration.Molecules[j].OrderParameter;
        }
      }
      OrientationProfile->values[i] = countMoleculesInsideBin > 0 ? OrderParameterSum/countMoleculesInsideBin : 0.;
    }
  }else{
    printf("ERROR: Attemption to calculate orientation profile with single atom chain.\n");    
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