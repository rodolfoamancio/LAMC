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

/* ***************************************************************************************************************
 * Name       | GetPressureLongRangeCorrection
 * ---------------------------------------------------------------------------------------------------------------
 * Function   | Calculates the long-range correction to the total pressure
 * Parameters | - Configuration: The configuration of the system containing molecules and atoms.
 * Returns    | The calculated long-range correction to the total potential energy.
 * **************************************************************************************************************/
double GetPressureLongRangeCorrection(CONFIGURATION Configuration){
  double PressureLongRangeCorrection = 0.0;
  enum CarbonType TypeA, TypeB;
  if ((ReferencePotential == MIE || ReferencePotential == BARKER_HENDERSON_PERTURBED)) {
    int NumberPesudoAtoms[NUMBER_PSEUDO_ATOMS_TYPES] = {0, 0, 0, 0};
    double VolumeCubicMeters = SimulationBox.volume / Cube(METER_TO_ANGSTRON);
    double AuxInteractions = 0;

    for (int i = 0; i < Configuration.NumberMolecules; i++) {
      for (int j = 0; j < Configuration.Molecules[i].Size; j++) {
        switch (Configuration.Molecules[i].Atoms[j].Type) {
          case CH4:
            NumberPesudoAtoms[0]++;
            break;
          
          case CH3e:
            NumberPesudoAtoms[1]++;
            break;
          
          case CH3:
            NumberPesudoAtoms[2]++;
            break;
          
          case CH2:
            NumberPesudoAtoms[3]++;
            break;
        }
      }
    }

    for (int i = 0; i < NUMBER_PSEUDO_ATOMS_TYPES; i++) {
      for (int j = i; j < NUMBER_PSEUDO_ATOMS_TYPES; j++) {
        if(i == 0){
          TypeA = CH4;
        }else if(i == 1){
          TypeA = CH3e;
        }else if(i == 2){
          TypeA = CH3;
        }else{
          TypeA = CH2;
        }

        if(j == 0){
          TypeB = CH4;
        }else if(j == 1){
          TypeB = CH3e;
        }else if(j == 2){
          TypeB = CH3;
        }else{
          TypeB = CH2;
        }

        double Sigma = GetSigmaCombination(GetAlkaneSigma(TypeA), GetAlkaneSigma(TypeB));
        double Epsilon = GetEpsilonCombination(GetAlkaneEpsilon(TypeA), GetAlkaneEpsilon(TypeB));
        double RepulsiveExponent = GetExponentCombination(
          GetAlkaneRepulsiveExponent(TypeA), 
          GetAlkaneRepulsiveExponent(TypeB)
        );
        double AttractiveExponent = GetExponentCombination(
          GetAlkaneAttractiveExponent(TypeA), 
          GetAlkaneAttractiveExponent(TypeB)
        );
        double C = GetCMie(RepulsiveExponent, AttractiveExponent);
        double SigmaCutoffN = pow(Sigma/CUTOFF_DISTANCE, RepulsiveExponent);
        double SigmaCutoffM = pow(Sigma/CUTOFF_DISTANCE, AttractiveExponent);
        AuxInteractions += (
          NumberPesudoAtoms[i]
          *NumberPesudoAtoms[j]
          *C
          *Epsilon
          *(
            (AttractiveExponent/(3-AttractiveExponent))*SigmaCutoffM
            -(RepulsiveExponent/(3-RepulsiveExponent))*SigmaCutoffN
          )
        );
      }
    }

    PressureLongRangeCorrection = 2.0*M_PI*AuxInteractions*Cube(CUTOFF_DISTANCE*ANGSTRON)/(3.*Squared(VolumeCubicMeters));
  }

  return PressureLongRangeCorrection;
}

/* ***************************************************************************************************************
 * Name       | GetPressureIdealGas
 * ---------------------------------------------------------------------------------------------------------------
 * Function   | Calculates the pressure of an ideal gas using the 
 *              given number of molecules and volume.
 * Parameters |
 *              - numberOfMolecules: The number of molecules in the 
 *                system.
 *              - volume: The volume of the system in cubic meters.
 * Returns    | 
 *              The calculated pressure of the ideal gas.
 * **************************************************************************************************************/
double GetPressureIdealGas(int numberOfMolecules, double volume){
  volume = volume*Cube(ANGSTRON);
  return numberOfMolecules*BOLTZMANN_CONSTANT*Temperature/volume;
}

/* ***************************************************************************************************************
 * Name       | GetPressureExcess
 * ---------------------------------------------------------------------------------------------------------------
 * Function   | Calculates the excess pressure of the system
 * Parameters | Configuration: CONFIGURATION struc withe the
 *              system's configuration
 * Returns    | The excess pressure
 * **************************************************************************************************************/
double GetPressureExcess(CONFIGURATION Configuration){
  double PressureExcess = 0.0;

  if(ReferencePotential == MIE || ReferencePotential == BARKER_HENDERSON_PERTURBED){
    StrainDerivativeTensor.xx=StrainDerivativeTensor.xy=StrainDerivativeTensor.xz=0.0;
    StrainDerivativeTensor.yx=StrainDerivativeTensor.yy=StrainDerivativeTensor.yz=0.0;
    StrainDerivativeTensor.zx=StrainDerivativeTensor.zy=StrainDerivativeTensor.zz=0.0;
    
    ComputeStrainDerivativeTensor(&Configuration);

    PressureExcess = -1.0*MatrixTrace(StrainDerivativeTensor);
  }
  return PressureExcess;
}