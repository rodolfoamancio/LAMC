/* ***************************************************************************************************************
 * This code was written as a part of LAMC, a program to perform molecular Monte Carlo simulations of linear
 * alkanes on multiple ensembles.
 * 
 * This module potential.c is a code file for potential functions
 * 
 * Developer: Rodolfo J Amancio (rodolfojamancio@gmail.com)
 * University of Campinas, School of Chemical Engineering, Campinas/SP -\sBrazil
 * 
 * **************************************************************************************************************/
#include"potential.h"

enum PotentialType ReferencePotential, PerturbationPotential;

/* ***************************************************************************************************************
 * Name       | GetCMie
 * ---------------------------------------------------------------------------------------------------------------
 * Function   | Calculates C pre-factor of the Mie potential
 * Parameters | 
 *              - RepulsiveExponent: The exponent for the repulsive part
 *              - AttractiveExponent: The exponent for the repulsive part
 * Returns    | The value of C
 * **************************************************************************************************************/
double GetCMie(double RepulsiveExponent, double AttractiveExponent){
  double C;
  C = (
    (RepulsiveExponent/(RepulsiveExponent - AttractiveExponent))
    *pow(RepulsiveExponent/AttractiveExponent, AttractiveExponent/(RepulsiveExponent - AttractiveExponent))
  );
  return C;
}

/* ***************************************************************************************************************
 * Name       | GetPotentialMie
 * ---------------------------------------------------------------------------------------------------------------
 * Function   | Calculates the Mie potential
 * Parameters | 
 *              - RepulsiveExponent: The exponent for the repulsive part
 *              - AttractiveExponent: The exponent for the repulsive part
 *              - Sigma: The segment diameter in A
 *              - Epsilon: the potential depth in J
 *              - Distance: the distance in A
 * Returns    | The potential energy in atomic units.
 * **************************************************************************************************************/
double GetPotentialMie(
  double RepulsiveExponent, 
  double AttractiveExponent, 
  double Sigma, 
  double Epsilon, 
  double Distance){
    double C;
    double SigmaDistance, SigmaDistanceN, SigmaDistanceM;
    C = GetCMie(RepulsiveExponent, AttractiveExponent);
    SigmaDistance = Sigma/Distance;
    SigmaDistanceN = pow(SigmaDistance, RepulsiveExponent);
    SigmaDistanceM = pow(SigmaDistance, AttractiveExponent);
    return C*Epsilon*(SigmaDistanceN - SigmaDistanceM);
}


/* ***************************************************************************************************************
 * Name       | GetPotentialStretching
 * ---------------------------------------------------------------------------------------------------------------
 * Function   | Calculates the potential energy due to stretching.
 * Parameters | - d: The displacement from the equilibrium distance.
 * Returns    | The potential energy in atomic units.
 * **************************************************************************************************************/
double GetPotentialStreaching(double d){
  return (kStreaching/2.0)*pow(d - dEq, 2);
}

/* ***************************************************************************************************************
 * Name       | GetPotentialBending
 * ---------------------------------------------------------------------------------------------------------------
 * Function   | Calculates the potential energy due to bending.
 * Parameters | - theta: The angular displacement from equilibrium.
 * Returns    | The potential energy in atomic units.
 * **************************************************************************************************************/
double GetPotentialBending(double theta){
  return (kBending/2.0)*pow(theta - thetaEq, 2);
}

/* ***************************************************************************************************************
 * Name       | GetPotentialTorsion
 * ---------------------------------------------------------------------------------------------------------------
 * Function   | Calculates the potential energy due to torsion.
 * Parameters | - phi: Torsinal angle.
 * Returns    | The potential energy in atomic units.
 * **************************************************************************************************************/
double GetPotentialTorsion(double phi){
  return (cTorsion[0]*(1 + cos(phi)) + cTorsion[1]*(1 - cos(2*phi)) + cTorsion[2]*(1 + cos(3*phi)));
}

/* ***************************************************************************************************************
 * Name       | GetPotentialBonded
 * ---------------------------------------------------------------------------------------------------------------
 * Function   | Calculates the potential energy of a bonded configuration.
 * Parameters | - Configuration: The configuration containing the molecules and their atoms.
 * Returns    | The total potential energy of the bonded configuration.
 * **************************************************************************************************************/
double GetPotentialBonded(CONFIGURATION Configuration) {
  double total = 0;
  
  for(int i=0; i<Configuration.NumberMolecules; ++i) {
    MOLECULE molecule = Configuration.Molecules[i];
    // loop for molecular bending energy potential
    for(int j=0; j<ChainSize-2; ++j) {
      ATOM atom1 = molecule.Atoms[j];
      ATOM atom2 = molecule.Atoms[j+1];
      ATOM atom3 = molecule.Atoms[j+2];
      
      double BendingAngle = CalculateBendingAngle(atom1.Position, atom2.Position, atom3.Position);
      total += GetPotentialBending(BendingAngle);
    }
    
    // loop for molecular torsion energy potential
    for(int j=0; j<ChainSize-3; ++j) {
      ATOM atom1 = molecule.Atoms[j];
      ATOM atom2 = molecule.Atoms[j+1];
      ATOM atom3 = molecule.Atoms[j+2];
      ATOM atom4 = molecule.Atoms[j+3];
      
      double TorsionAngle = CalculateTorsionAngle(atom1.Position, atom2.Position, atom3.Position, atom4.Position);
      total += GetPotentialTorsion(TorsionAngle);
    }
  }
  
  return total;
}

/* ***************************************************************************************************************
 * Name       | SampleBondLength
 * ---------------------------------------------------------------------------------------------------------------
 * Function   | Samples a bond length from the corresponding carbon types.
 * Returns    | The sampled bond length.
 * **************************************************************************************************************/
double SampleBondLength(enum CarbonType TypeA, enum CarbonType TypeB){
  if(TypeA == CH2 && TypeB == CH2){
    return 1.54;
  }else if (
    (TypeA == CH2 && TypeB == CH3)
    || (TypeA == CH3 && TypeB == CH2)
  ){
    return 1.74;  
  }else if(TypeA == CH3e && TypeB == CH3e){
    return 1.94;
  }
}

/* ***************************************************************************************************************
 * Name       | SampleBendingAngle
 * ---------------------------------------------------------------------------------------------------------------
 * Function   | Samples an orientation vector with a bending angle following a certain probability distribution.
 * Parameters | - r1: Vector representing the position of atom 1.
 *              - r2: Vector representing the position of atom 2.
 * Returns    | The sampled orientation vector.
 * **************************************************************************************************************/
VECTOR SampleBendingAngle(VECTOR r1, VECTOR r2){
  bool ready = false;
  VECTOR OrientationVector;
  double theta; 
  while(!ready){
    OrientationVector = RandomUnitVector();
    theta = CalculateBendingAngle(r1, r2, VectorSum(r2, OrientationVector));
    ready = GetRandomNumber() < exp(-Beta*GetPotentialBending(theta)) ? true : false;
  }
  return OrientationVector;
}

/* ***************************************************************************************************************
 * Name       | SampleBendingTorsionAngles
 * ---------------------------------------------------------------------------------------------------------------
 * Function   | Samples an orientation vector with bending and torsion angles following a certain probability 
 *              distribution.
 * Parameters | - r1: Vector representing the position of atom 1.
 *              - r2: Vector representing the position of atom 2.
 *              - r3: Vector representing the position of atom 3.
 * Returns    | The sampled orientation vector.
 * **************************************************************************************************************/
VECTOR SampleBendingTorsionAngles(VECTOR r1, VECTOR r2, VECTOR r3){
  bool ready = false;
  VECTOR OrientationVector, r4;
  double theta, phi, BendingPotential, TorsionPotential, TotalPotential; 
  while(!ready){
    OrientationVector = RandomUnitVector();
    r4 = VectorSum(r3, OrientationVector);
    theta = CalculateBendingAngle(r2, r3, r4);
    phi = CalculateTorsionAngle(r1, r2, r3, r4);
    ready = GetRandomNumber() < exp(-Beta*(GetPotentialBending(theta) + GetPotentialTorsion(phi))) ? true : false;
  }
  return OrientationVector;
}

/* ***************************************************************************************************************
 * Name       | GetPartialExternalPotential
 * ---------------------------------------------------------------------------------------------------------------
 * Function   | Calculates the partial external potential for a given reference molecule and reference particle.
 * Parameters | - Configuration: The configuration of molecules.
 *              - referenceMolecule: The index of the reference molecule.
 *              - referenceParticle: The index of the reference particle.
 * Returns    | The calculated partial external potential.
 * **************************************************************************************************************/
POTENTIAL GetPartialExternalPotential(CONFIGURATION Configuration, int referenceMolecule, int referenceParticle){
  POTENTIAL PartialPotential;
  VECTOR SeparationVector;
  double sigma, epsilon, RepulsiveExponent, AttractiveExponent;
  double Distance;
  double potential;
  bool ExcludedInterations;

  PartialPotential.overlap = false;
  PartialPotential.potential = 0.0;

  PartialPotential = GetPotentialExternalField(Configuration.Molecules[referenceMolecule].Atoms[referenceParticle]);
  if(PartialPotential.overlap) return PartialPotential;

  for(int i=0; i<Configuration.NumberMolecules; i++){
    for(int j=0; j<Configuration.Molecules[i].Size; j++){
      ExcludedInterations = (referenceMolecule==i) && (referenceParticle-j<4);
      if((!ExcludedInterations) && (!PartialPotential.overlap)){
        SeparationVector = VectorSubtraction(
          Configuration.Molecules[referenceMolecule].Atoms[referenceParticle].Position, 
          Configuration.Molecules[i].Atoms[j].Position
        );
        SeparationVector = ApplyPeriodicBoundaryConditionsVector(SeparationVector);
        Distance = Norm(SeparationVector);
        sigma = GetSigmaCombination(
          Configuration.Molecules[referenceMolecule].Atoms[referenceParticle].Sigma, 
          Configuration.Molecules[i].Atoms[j].Sigma
        );
        if((ReferencePotential == HARD_SPHERE) && (Distance < sigma)){
          PartialPotential.overlap = true;
          PartialPotential.potential = 1E30;
          return PartialPotential;
        }else if((ReferencePotential == MIE) && (Distance < CUTOFF_DISTANCE)){
          epsilon = GetEpsilonCombination(
              Configuration.Molecules[referenceMolecule].Atoms[referenceParticle].Epsilon, 
              Configuration.Molecules[i].Atoms[j].Epsilon
          );
          RepulsiveExponent = GetExponentCombination(
            Configuration.Molecules[referenceMolecule].Atoms[referenceParticle].RepulsiveExponent, 
            Configuration.Molecules[i].Atoms[j].RepulsiveExponent
          );
          AttractiveExponent = GetExponentCombination(
            Configuration.Molecules[referenceMolecule].Atoms[referenceParticle].AttractiveExponent, 
            Configuration.Molecules[i].Atoms[j].AttractiveExponent
          );
          potential = GetPotentialMie(
            RepulsiveExponent,
            AttractiveExponent,
            sigma,
            epsilon,
            Distance
          );
          if(Beta*potential > 10){
            PartialPotential.overlap = true;
            PartialPotential.potential = 1E30;
            return PartialPotential;
          }
          PartialPotential.potential += potential;
        }
      }
    }  
  }

  return PartialPotential;
}

/* ***************************************************************************************************************
 * Name       | GetPotentialExternalField
 * ---------------------------------------------------------------------------------------------------------------
 * Function   | Calculates the external potential for a given atom.
 * Parameters | - Atom: The atom for which to calculate the potential.
 * Returns    | The calculated external potential.
 * **************************************************************************************************************/
POTENTIAL GetPotentialExternalField(ATOM Atom) {
  POTENTIAL potentialField;
  potentialField.potential = 0.0;
  potentialField.overlap = false;

  if (SimulationBox.ClosedBox) {
    double z = Atom.Position.z;
    
    if (z <= SimulationBox.zMin || z >= SimulationBox.zMax){
      potentialField.overlap = true;
      potentialField.potential = 1E5;
      return potentialField;
    }

    double distanceLowerWall = z - SimulationBox.zMin;
    double distanceUpperWall = SimulationBox.zMax - z;
    double potentialLowerWall = GetPotentialSteele(Atom, distanceLowerWall);
    double potentialUpperWall = GetPotentialSteele(Atom, distanceUpperWall);

    potentialField.potential = potentialLowerWall + potentialUpperWall;
    potentialField.overlap = (Beta * potentialField.potential > 1000);
  }
  
  return potentialField;
}

/* ***************************************************************************************************************
 * Name       | GetPotentialGradient
 * ---------------------------------------------------------------------------------------------------------------
 * Function   | Calculates the gradient of the potential between two atoms.
 * Parameters | - SeparationVector: Vector representing the distance between AtomA and AtomB.
 *              - AtomA: First atom.
 *              - AtomB: Second atom.
 * Returns    | The calculated potential gradient vector.
 * **************************************************************************************************************/
VECTOR GetPotentialGradient(VECTOR SeparationVector, ATOM AtomA, ATOM AtomB) {
  VECTOR PotentialGradient;
  double Sigma, Epsilon, AttractiveExponent, RepulsiveExponent, C;
  double Distance;
  double SigmaOverDistanceM, SigmaOverDistanceN;

  Distance = Norm(SeparationVector);

  if (Distance < CUTOFF_DISTANCE) {
    Sigma = GetSigmaCombination(AtomA.Sigma, AtomB.Sigma);
    Epsilon = GetEpsilonCombination(AtomA.Epsilon, AtomB.Epsilon);
    AttractiveExponent = GetExponentCombination(AtomA.AttractiveExponent, AtomB.AttractiveExponent);
    RepulsiveExponent = GetExponentCombination(AtomA.RepulsiveExponent, AtomB.RepulsiveExponent);
    C = GetCMie(RepulsiveExponent, AttractiveExponent);
    SigmaOverDistanceM = pow(Sigma/Distance, AttractiveExponent);
    SigmaOverDistanceN = pow(Sigma/Distance, RepulsiveExponent);

    double PotentialGradientModule = (
      C
      *Epsilon
      *(1/(Distance*ANGSTRON))
      *(
        AttractiveExponent*SigmaOverDistanceM
        -RepulsiveExponent*SigmaOverDistanceN
      )
    );
    PotentialGradient.x = PotentialGradientModule*SeparationVector.x/Distance; 
    PotentialGradient.y = PotentialGradientModule*SeparationVector.y/Distance; 
    PotentialGradient.z = PotentialGradientModule*SeparationVector.z/Distance;
  } else {
    PotentialGradient.x = 0.0;
    PotentialGradient.y = 0.0;
    PotentialGradient.z = 0.0;
  }

  return PotentialGradient;
}


/* ***************************************************************************************************************
 * Name       | GetPotentialSteele
 * ---------------------------------------------------------------------------------------------------------------
 * Function   | Calculates the Steele potential for an atom.
 * Parameters | - Atom: Atom for which the Steele potential is calculated.
 *              - height: Height parameter for the calculation.
 * Returns    | The calculated Steele potential.
 * **************************************************************************************************************/
double GetPotentialSteele(ATOM Atom, double height) {
  double carbonDensity = 114.0E27;
  double interlayerSpacing = 3.35 * ANGSTRON;
  double sigmaCarbon = 3.4, epsilonCarbon = 28 * BOLTZMANN_CONSTANT;
  double sigma, epsilon, binaryInteractionParameter;

  binaryInteractionParameter = -0.035;

  sigma = GetSigmaCombination(Atom.Sigma, sigmaCarbon)*ANGSTRON;
  epsilon = (1-binaryInteractionParameter)*GetEpsilonCombination(Atom.Epsilon, epsilonCarbon);

  height = height*ANGSTRON;

  double term1 = 0.4*pow(sigma/height, 10);
  double term2 = pow(sigma/height, 4);
  double term3 = pow(sigma, 4)/(3*interlayerSpacing*pow(0.61*interlayerSpacing + height, 3));

  double SteelePotential = 2*M_PI*carbonDensity*epsilon*pow(sigma, 2)*interlayerSpacing*(term1 - term2 - term3);

  return SteelePotential;
}

/* ***************************************************************************************************************
 * Name       | GetPotentialNonbonded
 * ---------------------------------------------------------------------------------------------------------------
 * Function   | Calculates the potential energy between non-bonded atoms using the Lennard-Jones potential.
 * Parameters | - Configuration: The configuration of the system containing molecules and atoms.
 *              - Potential: The type of potential to be calculated (Lennard-Jones in this case).
 * Returns    | The total potential energy calculated.
 * **************************************************************************************************************/
double GetPotentialNonbonded(CONFIGURATION Configuration, enum PotentialType Potential) {
  VECTOR SeparationVector;
  double Sum = 0.0;
  double RepulsiveExponent, AttractiveExponent;
  double potential;
  double Distance;
  double sigma, epsilon;
  int i, j, k, l, k0, l0;

  for (i = 0; i < Configuration.NumberMolecules; i++) {
    for (j = 0; j < Configuration.Molecules[i].Size; j++) {
      k0 = j < Configuration.Molecules[i].Size - 4 ? i : i + 1;
      for (k = k0; k < Configuration.NumberMolecules; k++) {
        l0 = k == i ? j + 4 : 0;
        for (l = l0; l < Configuration.Molecules[k].Size; l++) {
          SeparationVector = VectorSubtraction(
            Configuration.Molecules[i].Atoms[j].Position, 
            Configuration.Molecules[k].Atoms[l].Position
          );
          SeparationVector = ApplyPeriodicBoundaryConditionsVector(SeparationVector);
          Distance = Norm(SeparationVector);
          sigma = GetSigmaCombination(Configuration.Molecules[i].Atoms[j].Sigma, Configuration.Molecules[k].Atoms[l].Sigma);

          if (Potential == MIE && Distance < CUTOFF_DISTANCE){
            epsilon = GetEpsilonCombination(
              Configuration.Molecules[i].Atoms[j].Epsilon, 
              Configuration.Molecules[k].Atoms[l].Epsilon
            );
            RepulsiveExponent = GetExponentCombination(
              Configuration.Molecules[i].Atoms[j].RepulsiveExponent, 
              Configuration.Molecules[k].Atoms[l].RepulsiveExponent
            );
            AttractiveExponent = GetExponentCombination(
              Configuration.Molecules[i].Atoms[j].AttractiveExponent, 
              Configuration.Molecules[k].Atoms[l].AttractiveExponent
            );
            potential = GetPotentialMie(
              RepulsiveExponent,
              AttractiveExponent,
              sigma,
              epsilon,
              Distance
            );
            Sum += potential;
            if(Beta*potential > 1000) return 1E6;
          }else if (Potential == HARD_SPHERE && Distance < sigma){
            return 1E6;
          }
        }
      }
    }
  }
  return Sum;
}

/* ***************************************************************************************************************
 * Name       | GetTotalPotentialExternal
 * ---------------------------------------------------------------------------------------------------------------
 * Function   | Calculates the total potential energy between external atoms using the specified potential field.
 * Parameters | - Configuration: The configuration of the systemcontaining molecules and atoms.
 * Returns    | The total potential energy calculated.
 * **************************************************************************************************************/
double GetTotalPotentialExternal(CONFIGURATION Configuration){
  double total = 0;
  bool overlap = false;
  if(SimulationBox.ClosedBox){
    for(int i = 0; i < Configuration.NumberMolecules; i++){
      for(int j = 0; j < ChainSize; j++)
      total += GetPotentialExternalField(Configuration.Molecules[i].Atoms[j]).potential;
    }
  }  
  return total;
}

/* ***************************************************************************************************************
 * Name       | GetPotentialLongRangeCorrection
 * ---------------------------------------------------------------------------------------------------------------
 * Function   | Calculates the long-range correction to the total potential energy between atoms in the system 
 *              using the specified potential field.
 * Parameters | - Configuration: The configuration of the system containing molecules and atoms.
 *            | - Potential: PotentialType for calculation
 * Returns    | The calculated long-range correction to the total potential energy.
 * **************************************************************************************************************/
double GetPotentialLongRangeCorrection(CONFIGURATION Configuration, enum PotentialType Potential) {
  double PotentialLongRangeCorrection = 0.0;
  enum CarbonType TypeA, TypeB;

  if (Potential == MIE) {
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
        double RepulsiveExponent = GetExponentCombination(GetAlkaneRepulsiveExponent(TypeA), GetAlkaneRepulsiveExponent(TypeB));
        double AttractiveExponent = GetExponentCombination(GetAlkaneAttractiveExponent(TypeA), GetAlkaneAttractiveExponent(TypeB));
        double C = GetCMie(RepulsiveExponent, AttractiveExponent);
        double SigmaCutoffN = pow(Sigma/CUTOFF_DISTANCE, RepulsiveExponent);
        double SigmaCutoffM = pow(Sigma/CUTOFF_DISTANCE, AttractiveExponent);
        AuxInteractions += (
          NumberPesudoAtoms[i]
          *NumberPesudoAtoms[j]
          *Epsilon
          *C
          *(
            (1.0/(3 - RepulsiveExponent))*SigmaCutoffN
            -(1.0/(3 - AttractiveExponent))*SigmaCutoffM
          )
        );
      }
    }

    PotentialLongRangeCorrection = -2.0*M_PI*AuxInteractions*Cube(CUTOFF_DISTANCE*ANGSTRON)/VolumeCubicMeters;
  }
  
  return PotentialLongRangeCorrection;
}

/* ***************************************************************************************************************
 * Name       | ComputeStrainDerivativeTensor
 * ---------------------------------------------------------------------------------------------------------------
 * Function   | Calculates the strain derivative tensor with the appropriate forces based on the reference
 *            | potential for nonbonded interactions
 * Parameters | - Configuration: The configuration of the system containing molecules and atoms.
 * Returns    | None
 * **************************************************************************************************************/
void ComputeStrainDerivativeTensor(CONFIGURATION *Configuration){
  VECTOR Forcejl, Force;
  VECTOR SeparationVectorlj, CenterOfMass, Position;
  double VolumeMeters = SimulationBox.volume*Cube(ANGSTRON), Average;

  for(int i=0; i<Configuration->NumberMolecules; i++){
    for(int j=0; j<Configuration->Molecules[i].Size; j++){
      Configuration->Molecules[i].Atoms[j].Force.x = 0;
      Configuration->Molecules[i].Atoms[j].Force.y = 0;
      Configuration->Molecules[i].Atoms[j].Force.z = 0;
    }
  }

  // atomic contribution
  for(int i=0; i<Configuration->NumberMolecules; i++){
    for(int j=0; j<Configuration->Molecules[i].Size; j++){
      for(int k=i+1; k<Configuration->NumberMolecules; k++){
        for(int l=0; l<Configuration->Molecules[k].Size; l++){
          SeparationVectorlj = VectorSubtraction(
            Configuration->Molecules[i].Atoms[j].Position,
            Configuration->Molecules[k].Atoms[l].Position
          );
          SeparationVectorlj = ApplyPeriodicBoundaryConditionsVector(SeparationVectorlj);
          Forcejl = GetPotentialGradient(
            SeparationVectorlj, 
            Configuration->Molecules[i].Atoms[j], 
            Configuration->Molecules[k].Atoms[l]
          );
          Configuration->Molecules[i].Atoms[j].Force = VectorSubtraction(
            Configuration->Molecules[i].Atoms[j].Force, Forcejl
          );
          Configuration->Molecules[k].Atoms[l].Force = VectorSum(
            Configuration->Molecules[k].Atoms[l].Force, Forcejl
          );

          StrainDerivativeTensor.xx += Forcejl.x*SeparationVectorlj.x*ANGSTRON;
          StrainDerivativeTensor.yx += Forcejl.y*SeparationVectorlj.x*ANGSTRON;
          StrainDerivativeTensor.zx += Forcejl.z*SeparationVectorlj.x*ANGSTRON;

          StrainDerivativeTensor.xy += Forcejl.x*SeparationVectorlj.y*ANGSTRON;
          StrainDerivativeTensor.yy += Forcejl.y*SeparationVectorlj.y*ANGSTRON;
          StrainDerivativeTensor.zy += Forcejl.z*SeparationVectorlj.y*ANGSTRON;

          StrainDerivativeTensor.xz += Forcejl.x*SeparationVectorlj.z*ANGSTRON;
          StrainDerivativeTensor.yz += Forcejl.y*SeparationVectorlj.z*ANGSTRON;
          StrainDerivativeTensor.zz += Forcejl.z*SeparationVectorlj.z*ANGSTRON;
        }
      }
    }
  }

  // molecular correction
  for(int i = 0; i<Configuration->NumberMolecules; i++){
    CenterOfMass = GetMoleculeCenterOfMass(Configuration->Molecules[i]);
    for(int j = 0; j<Configuration->Molecules[i].Size; j++){
      Position = Configuration->Molecules[i].Atoms[j].Position;
      Force = Configuration->Molecules[i].Atoms[j].Force;
      StrainDerivativeTensor.xx += Force.x*(Position.x - CenterOfMass.x)*ANGSTRON;
      StrainDerivativeTensor.yx += Force.y*(Position.x - CenterOfMass.x)*ANGSTRON;
      StrainDerivativeTensor.zx += Force.z*(Position.x - CenterOfMass.x)*ANGSTRON;

      StrainDerivativeTensor.xy += Force.x*(Position.y - CenterOfMass.y)*ANGSTRON;
      StrainDerivativeTensor.yy += Force.y*(Position.y - CenterOfMass.y)*ANGSTRON;
      StrainDerivativeTensor.zy += Force.z*(Position.y - CenterOfMass.y)*ANGSTRON;

      StrainDerivativeTensor.xz += Force.x*(Position.z - CenterOfMass.z)*ANGSTRON;
      StrainDerivativeTensor.yz += Force.y*(Position.z - CenterOfMass.z)*ANGSTRON;
      StrainDerivativeTensor.zz += Force.z*(Position.z - CenterOfMass.z)*ANGSTRON;
    }
  }

  StrainDerivativeTensor.xx /= 3*VolumeMeters;
  StrainDerivativeTensor.yx /= 3*VolumeMeters;
  StrainDerivativeTensor.zx /= 3*VolumeMeters;

  StrainDerivativeTensor.xy /= 3*VolumeMeters;
  StrainDerivativeTensor.yy /= 3*VolumeMeters;
  StrainDerivativeTensor.zy /= 3*VolumeMeters;

  StrainDerivativeTensor.xz /= 3*VolumeMeters;
  StrainDerivativeTensor.yz /= 3*VolumeMeters;
  StrainDerivativeTensor.zz /= 3*VolumeMeters;

  Average = 0.5*(StrainDerivativeTensor.xy+StrainDerivativeTensor.yx);
  StrainDerivativeTensor.xy = StrainDerivativeTensor.yx = Average;
  Average = 0.5*(StrainDerivativeTensor.xz+StrainDerivativeTensor.zx);
  StrainDerivativeTensor.xz = StrainDerivativeTensor.zx = Average;
  Average = 0.5*(StrainDerivativeTensor.yz+StrainDerivativeTensor.zy);
  StrainDerivativeTensor.yz = StrainDerivativeTensor.zy = Average;

}
