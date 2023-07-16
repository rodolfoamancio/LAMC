/* ***************************************************************************************************************
 * This code was written as a part of LAMC, a program to perform molecular Monte Carlo simulations of linear
 * alkanes on multiple ensembles.
 * 
 * This module properties.h is a code file for calculating properties
 * 
 * Developer: Rodolfo J Amancio (rodolfojamancio@gmail.com)
 * University of Campinas, School of Chemical Engineering, Campinas/SP -	Brazil
 * 
 * **************************************************************************************************************/
#include"properties.h"

/***************************************************************
 * Name       | GetPotentialLongRangeCorrection
 * -------------------------------------------------------------
 * Function   | Calculates the long-range correction to the total 
 *              potential energy between atoms in the system using 
 *              the specified potential field.
 * Parameters |
 *              - Configuration: The configuration of the system
 *                containing molecules and atoms.
 * Returns    | 
 *              The calculated long-range correction to the total 
 *              potential energy.
 ****************************************************************/
double GetPressureLongRangeCorrection(CONFIGURATION Configuration, double* EpsilonAlkane, double* SigmaAlkane) {
  ATOM   AtomA, AtomB;
  int    NumberPesudoAtoms[3] = {0, 0, 0};
  double VolumeCubicMeters = SimulationBox.volume * Cube(ANGSTRON);
  double AuxInteractions = 0;

  // Count the number of pseudo atoms of each type
  for (int i = 0; i < Configuration.NumberMolecules; i++) {
    for (int j = 0; j < Configuration.Molecules[i].Size; j++) {
      switch (Configuration.Molecules[i].Atoms[j].Type) {
        case CH4:
          NumberPesudoAtoms[0]++;
          break;

        case CH3:
          NumberPesudoAtoms[1]++;
          break;

        case CH2:
          NumberPesudoAtoms[2]++;
          break;
      }
    }
  }

  for (int i = 0; i < NUMBER_PSEUDO_ATOMS_TYPES; i++) {
    AtomA.Epsilon = EpsilonAlkane[i];
    AtomA.Sigma = SigmaAlkane[i];
    double EpsilonI = AtomA.Epsilon * Cube(ANGSTRON);
    double SigmaI = AtomA.Sigma / CUTOFF_DISTANCE;
    double SigmaOverCutoffI3 = Cube(SigmaI) / Cube(CUTOFF_DISTANCE);

    for (int j = 0; j < NUMBER_PSEUDO_ATOMS_TYPES; j++) {
      AtomB.Epsilon = EpsilonAlkane[j];
      AtomB.Sigma = SigmaAlkane[j];
      double EpsilonJ = AtomB.Epsilon * Cube(ANGSTRON);
      double SigmaJ = AtomB.Sigma / CUTOFF_DISTANCE;
      double SigmaOverCutoffJ3 = Cube(SigmaJ) / Cube(CUTOFF_DISTANCE);

      AuxInteractions += (NumberPesudoAtoms[i] * NumberPesudoAtoms[j] * GetEpsilon(EpsilonI, EpsilonJ) *
          Cube(AtomA.Sigma * ANGSTRON) * ((2.0/3.0) * SigmaOverCutoffJ3 - SigmaOverCutoffI3));
    }
  }

  return 16.0 * M_PI * AuxInteractions / (3.0 * Squared(VolumeCubicMeters));
}

/***************************************************************
 * Name       | GetPressureIdealGas
 * -------------------------------------------------------------
 * Function   | Calculates the pressure of an ideal gas using the 
 *              given number of molecules and volume.
 * Parameters |
 *              - numberOfMolecules: The number of molecules in the 
 *                system.
 *              - volume: The volume of the system in cubic meters.
 * Returns    | 
 *              The calculated pressure of the ideal gas.
 ****************************************************************/
double GetPressureIdealGas(int numberOfMolecules, double volume){
  volume = volume*Cube(ANGSTRON);
  return numberOfMolecules*BOLTZMANN_CONSTANT*Temperature/volume;
}

/***************************************************************
 * Name       | GetPressureExcess
 * -------------------------------------------------------------
 * Function   | Calculates the excess pressure of the system
 * Parameters | Configuration: CONFIGURATION struc withe the
 *              system's configuration
 * Returns    | The excess pressure
 ****************************************************************/
double GetPressureExcess(CONFIGURATION Configuration) {
  VECTOR Position, CenterOfMass, Force;
  double Average = 0.0;
  double VolumeMeters = SimulationBox.volume * Cube(ANGSTRON);
  double PressureExcess = 0.0;

  if (ReferencePotential == LENNARD_JONES) {
    StrainDerivativeTensor.xx = 0.0;
    StrainDerivativeTensor.yx = 0.0;
    StrainDerivativeTensor.zx = 0.0;

    StrainDerivativeTensor.xy = 0.0;
    StrainDerivativeTensor.yy = 0.0;
    StrainDerivativeTensor.zy = 0.0;

    StrainDerivativeTensor.xz = 0.0;
    StrainDerivativeTensor.yz = 0.0;
    StrainDerivativeTensor.zz = 0.0;
    
    ComputeNonbondedForces(&Configuration);

    for (int i = 0; i < Configuration.NumberMolecules; i++) {
      CenterOfMass = GetMoleculeCenterOfMass(Configuration.Molecules[i]);

      for (int j = 0; j < Configuration.Molecules[i].Size; j++) {
        Position = Configuration.Molecules[i].Atoms[j].Position;
        Force = Configuration.Molecules[i].Atoms[j].Force;

        double dx = (Position.x - CenterOfMass.x) * ANGSTRON;
        double dy = (Position.y - CenterOfMass.y) * ANGSTRON;
        double dz = (Position.z - CenterOfMass.z) * ANGSTRON;

        StrainDerivativeTensor.xx += Force.x * dx;
        StrainDerivativeTensor.yx += Force.x * dy;
        StrainDerivativeTensor.zx += Force.x * dz;

        StrainDerivativeTensor.xy += Force.y * dx;
        StrainDerivativeTensor.yy += Force.y * dy;
        StrainDerivativeTensor.zy += Force.y * dz;

        StrainDerivativeTensor.xz += Force.z * dx;
        StrainDerivativeTensor.yz += Force.z * dy;
        StrainDerivativeTensor.zz += Force.z * dz;
      }
    }

    Average = 0.5 * (StrainDerivativeTensor.xy + StrainDerivativeTensor.yx);
    StrainDerivativeTensor.xy = StrainDerivativeTensor.yx = Average;

    Average = 0.5 * (StrainDerivativeTensor.xz + StrainDerivativeTensor.zx);
    StrainDerivativeTensor.xz = StrainDerivativeTensor.zx = Average;

    Average = 0.5 * (StrainDerivativeTensor.yz + StrainDerivativeTensor.zy);
    StrainDerivativeTensor.yz = StrainDerivativeTensor.zy = Average;

    PressureExcess = -MatrixTrace(StrainDerivativeTensor) / (3.0 * VolumeMeters);
  }

  return PressureExcess;
}