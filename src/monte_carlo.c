/* ***************************************************************************************************************
 * This code was written as a part of LAMC, a program to perform molecular Monte Carlo simulations of linear
 * alkanes on multiple ensembles.
 * 
 * This module files_io.c is a code file for defining files input output functions
 * 
 * Developer: Rodolfo J Amancio (rodolfojamancio@gmail.com)
 * University of Campinas, School of Chemical Engineering, Campinas/SP -\sBrazil
 * 
 * **************************************************************************************************************/

#include"monte_carlo.h"

enum ensemble SimulationEnsemble;
enum ensemble EquilibraionEnsemble;
enum ensemble Ensemble;

/* ***************************************************************************************************************
 * Name       | InitializeConfiguration
 * ---------------------------------------------------------------------------------------------------------------
 * Function   | Initializes the configuration with values for molecules
 * Parameters | Configuration: a pointer to the configuration struct
 * Returns    | None
 * **************************************************************************************************************/
void InitializeConfiguration(CONFIGURATION* Configuration){
  Configuration->NumberMolecules  = InitialNumberMolecules;
  Configuration->Molecules = (MOLECULE*) calloc(InitialNumberMolecules, sizeof(MOLECULE));
  Configuration->Molecules[0].Size = ChainSize;
  Configuration->Molecules[0].Atoms = (ATOM*) calloc(ChainSize, sizeof(ATOM));
  Configuration->Molecules[0].MolarMass = 0;

  for(int j=0; j<ChainSize; j++){
    if(ChainSize==1){
      Configuration->Molecules[0].Atoms[j].Type = CH4;
    }else{
      Configuration->Molecules[0].Atoms[j].Type = (
        (j==0 || j==Configuration->Molecules[0].Size-1) ? 
        CH3 : 
        CH2
      );
    }
    Configuration->Molecules[0].Atoms[j].Epsilon = GetAlkaneEpsilon(Configuration->Molecules[0].Atoms[j].Type);
    Configuration->Molecules[0].Atoms[j].Sigma = GetAlkaneSigma(Configuration->Molecules[0].Atoms[j].Type);
    Configuration->Molecules[0].Atoms[j].MolarMass = GetAlkaneAtomMolarMass(Configuration->Molecules[0].Atoms[j].Type);
    Configuration->Molecules[0].Atoms[j].AttractiveExponent = GetAlkaneAttractiveExponent(Configuration->Molecules[0].Atoms[j].Type);
    Configuration->Molecules[0].Atoms[j].RepulsiveExponent = GetAlkaneRepulsiveExponent(Configuration->Molecules[0].Atoms[j].Type);
    Configuration->Molecules[0].MolarMass += Configuration->Molecules[0].Atoms[j].MolarMass;
  }

  for(int i=1; i<InitialNumberMolecules; i++){
    Configuration->Molecules[i].Size = ChainSize;
    Configuration->Molecules[i].Atoms = (ATOM*) calloc(ChainSize, sizeof(ATOM));
    Configuration->Molecules[i].MolarMass = Configuration->Molecules[0].MolarMass;
    for(int j=0; j<ChainSize; j++) Configuration->Molecules[i].Atoms[j] = Configuration->Molecules[0].Atoms[j];
  }
}

/* ***************************************************************************************************************
 * Name       | CopyConfiguration
 * ---------------------------------------------------------------------------------------------------------------
 * Function   | Copies the content of the source configuration struct to the destination configuration struct, 
 *              including allocated memory
 * Parameters | - Source: the source configuration struct
 *            | - Dest: a pointer to the destination configuration struct
 * Returns    | None
 * **************************************************************************************************************/
void CopyConfiguration(CONFIGURATION Source, CONFIGURATION* Dest){
  for(int i=0; i<Dest->NumberMolecules; i++)
    free(Dest->Molecules[i].Atoms);
    free(Dest->Molecules);
    Dest->NumberMolecules = Source.NumberMolecules;
    Dest->Molecules = (MOLECULE*) calloc(Dest->NumberMolecules, sizeof(MOLECULE));
    for(int i = 0; i < Dest->NumberMolecules; i++){
      Dest->Molecules[i].Size = Source.Molecules[i].Size;
      Dest->Molecules[i].Atoms = (ATOM*) calloc(Dest->Molecules[i].Size, sizeof(ATOM));
      Dest->Molecules[i].MolarMass = Source.Molecules[i].MolarMass;
      for(int j = 0; j < ChainSize; j++){
        Dest->Molecules[i].Atoms[j] = Source.Molecules[i].Atoms[j];
      }
    }
}

/* ***************************************************************************************************************
 * Name       | SampleBeadPosition
 * ---------------------------------------------------------------------------------------------------------------
 * Function   | Calculate the position of a bead in a molecule based on its index
 * Parameters | - Molecule: The molecule containing the atoms and their positions
 *            | - Bead: The index of the bead for which to calculate the position
 * Returns    | The vector representing the position of the specified bead
 * **************************************************************************************************************/
VECTOR SampleBeadPosition(MOLECULE Molecule, int Bead) {
    double BondLength;
    VECTOR DirectionalUnitVector;

    // Sample bond length
    BondLength = SampleBondLength(
      Molecule.Atoms[Bead-1].Type,
      Molecule.Atoms[Bead].Type
    );

    if(Bead==1){
        // Random unit vector for first bead
        DirectionalUnitVector = RandomUnitVector();
    }else if(Bead == 2){
        // Calculate directional unit vector using bending angle
        DirectionalUnitVector = SampleBendingAngle(
          Molecule.Atoms[Bead-2].Position, 
          Molecule.Atoms[Bead-1].Position
        );
    }else if (Bead >= 3){
        // Calculate directional unit vector using bending torsion angles
        DirectionalUnitVector = SampleBendingTorsionAngles(
          Molecule.Atoms[Bead-3].Position, 
          Molecule.Atoms[Bead-2].Position, 
          Molecule.Atoms[Bead-1].Position
        );
    }

    // Multiply directional unit vector by bond length
    DirectionalUnitVector = MultiplyVectorScalar(DirectionalUnitVector, BondLength);

    // Return the sum of directional unit vector and position of the current bead
    return VectorSum(DirectionalUnitVector, Molecule.Atoms[Bead - 1].Position);
}


/* ***************************************************************************************************************
 * Name       | SelectTrialOrientation
 * ---------------------------------------------------------------------------------------------------------------
 * Function   | Selects a trial orientation based on the 
 *              weights provided and the total weight
 * Parameters | - WeightList: An array of weights for each 
 *              orientation
 *              - TotalWeight: The sum of all the weights in 
 *              WeightList
 * Returns    | The index of the selected trial orientation
 * **************************************************************************************************************/
int SelectTrialOrientation(double *WeightList, double TotalWeight){
  int i=0;
  double SelectedAccumulatedWeight, AccumulatedSum;
  SelectedAccumulatedWeight = GetRandomNumber()*TotalWeight;
  AccumulatedSum = WeightList[0];
  while(AccumulatedSum<SelectedAccumulatedWeight){
    i++;
    AccumulatedSum += WeightList[i];
  }
  return i;
}

/* ***************************************************************************************************************
 * Name       | GetRosenbluthWeightIdealChain
 * ---------------------------------------------------------------------------------------------------------------
 * Function   | Computes the Rsenbluth weight for an ideal chain
 * Parameters | Molecule (type MOLECULE) with molecular structure
 * Returns    | Roserbluth weight as double
 * **************************************************************************************************************/
double GetRosenbluthWeightIdealChain(MOLECULE Molecule){
  CONFIGURATION Configuration;
  POTENTIAL PartialPotential[MAX_CHAIN_SIZE][NUMBER_TRIAL_ORIENTATIONS];
  VECTOR DirectionalUnitVector;
  VECTOR PositionCandidates[MAX_CHAIN_SIZE][NUMBER_TRIAL_ORIENTATIONS];
  double WeightBeadTrials[MAX_CHAIN_SIZE][NUMBER_TRIAL_ORIENTATIONS] = {0};
  double WeightBeadTotal[NUMBER_TRIAL_ORIENTATIONS] = {0};
  double WeightChain;
  double BondLength;
  int SelectedTrialOrientation;
  bool SimulationBoxOriginalState;
  
  SimulationBoxOriginalState = SimulationBox.ClosedBox;
  SimulationBox.ClosedBox = false;

  Configuration.NumberMolecules = 1;
  Configuration.Molecules = (MOLECULE*) calloc(Configuration.NumberMolecules, sizeof(MOLECULE));
  Configuration.Molecules[0].Size = Molecule.Size;
  Configuration.Molecules[0].MolarMass = Molecule.MolarMass;
  Configuration.Molecules[0].Atoms = (ATOM*) calloc(Configuration.Molecules[0].Size, sizeof(ATOM));

  for(int i = 0; i <Configuration.Molecules[0].Size; i++)
    Configuration.Molecules[0].Atoms[i] = Molecule.Atoms[i];

  Configuration.Molecules[0].Atoms[0].Position = GetRandomPosition();
  PartialPotential[0][0] = GetPartialExternalPotential(Configuration, 0, 0);
  WeightBeadTrials[0][0] = !PartialPotential[0][0].overlap ? exp(-Beta*PartialPotential[0][0].potential) : 0.0;
  WeightBeadTotal[0] = WeightBeadTrials[0][0];
  WeightChain = WeightBeadTotal[0];

  for(int j=1; j<Configuration.Molecules[0].Size; j++){
    for(int k=0; k<NUMBER_TRIAL_ORIENTATIONS; k++){
      PositionCandidates[j][k] = SampleBeadPosition(Configuration.Molecules[0], j);
      Configuration.Molecules[0].Atoms[j].Position = PositionCandidates[j][k];
      PartialPotential[j][k] = GetPartialExternalPotential(Configuration, 0, j);
      WeightBeadTrials[j][k] = !PartialPotential[j][k].overlap ? exp(-Beta*PartialPotential[j][k].potential) : 0.0;
      WeightBeadTotal[j] += WeightBeadTrials[j][k];
    }
    SelectedTrialOrientation = SelectTrialOrientation(WeightBeadTrials[j], WeightBeadTotal[j]);
    Configuration.Molecules[0].Atoms[j].Position = PositionCandidates[j][SelectedTrialOrientation];
    WeightChain *= WeightBeadTotal[j]/NUMBER_TRIAL_ORIENTATIONS;
  }
  for(int i=0; i<Configuration.NumberMolecules; i++)
    free(Configuration.Molecules[i].Atoms);
  free(Configuration.Molecules);

  SimulationBox.ClosedBox = SimulationBoxOriginalState;

  return WeightChain;
}

/* ***************************************************************************************************************
 * Name       | GetRosenbluthWeightGhostMolecule
 * ---------------------------------------------------------------------------------------------------------------
 * Function   | Computers the Rosenbluth weight of a ghost 
 *              molecule in a given configuration
 * Parameters | Configuration (type CONFIGURATION)
 * Returns    | Rosenbluth weight as double
 * **************************************************************************************************************/
double GetRosenbluthWeightGhostMolecule(CONFIGURATION Configuration){
  CONFIGURATION TestConfiguration;
  POTENTIAL Potential;
  POTENTIAL	PartialPotential[MAX_CHAIN_SIZE][NUMBER_TRIAL_ORIENTATIONS];
  MOLECULE TestMolecule;
  VECTOR FirstBeadPosition;
  VECTOR DirectionalUnitVector;
  VECTOR PositionCandidates[MAX_CHAIN_SIZE][NUMBER_TRIAL_ORIENTATIONS];
  double WeightBeadTrials[MAX_CHAIN_SIZE][NUMBER_TRIAL_ORIENTATIONS] = {0};
  double WeightBeadTotal[NUMBER_TRIAL_ORIENTATIONS] = {0};
  double WeightChain;
  double BondLength;
  double DeltaPotentialLongRangeCorrection;
  bool FirstBeadPlaced = false;
  int SelectedTrialOrientation;
  int MoleculeIndex;


  // initializing test configuration
  TestConfiguration.NumberMolecules = Configuration.NumberMolecules+1;
  TestConfiguration.Molecules = (MOLECULE*) calloc(TestConfiguration.NumberMolecules, sizeof(MOLECULE));
  for(int i = 0; i < Configuration.NumberMolecules; i++){
    TestConfiguration.Molecules[i].Size  = Configuration.Molecules[i].Size;
    TestConfiguration.Molecules[i].MolarMass  = Configuration.Molecules[i].MolarMass;
    TestConfiguration.Molecules[i].Atoms = (ATOM*) calloc(TestConfiguration.Molecules[i].Size, sizeof(ATOM));
    for(int j = 0; j < TestConfiguration.Molecules[i].Size; j++)
      TestConfiguration.Molecules[i].Atoms[j] = Configuration.Molecules[i].Atoms[j];
  }
  
  MoleculeIndex = TestConfiguration.NumberMolecules-1;
  TestConfiguration.Molecules[MoleculeIndex].Size = Configuration.Molecules[0].Size;
  TestConfiguration.Molecules[MoleculeIndex].MolarMass = Configuration.Molecules[0].MolarMass;
  TestConfiguration.Molecules[MoleculeIndex].Atoms = (ATOM*) calloc(TestConfiguration.Molecules[0].Size, sizeof(ATOM));
  for(int j=0; j<TestConfiguration.Molecules[MoleculeIndex].Size; j++)
    TestConfiguration.Molecules[MoleculeIndex].Atoms[j] = Configuration.Molecules[0].Atoms[j];
    
  // insert molecule
  TestConfiguration.Molecules[MoleculeIndex].Atoms[0].Position = GetRandomPosition();
  Potential = GetPartialExternalPotential(TestConfiguration, MoleculeIndex, 0);
  WeightChain = !Potential.overlap ? exp(-Beta*Potential.potential) : 0.0;

  for(int j=1; j<TestConfiguration.Molecules[MoleculeIndex].Size; j++){
    for(int k=0; k<NUMBER_TRIAL_ORIENTATIONS; k++){
      PositionCandidates[j][k] = SampleBeadPosition(TestConfiguration.Molecules[MoleculeIndex], j);		
      TestConfiguration.Molecules[MoleculeIndex].Atoms[j].Position = PositionCandidates[j][k];
      PartialPotential[j][k] = GetPartialExternalPotential(TestConfiguration, MoleculeIndex, j);
      WeightBeadTrials[j][k] = !PartialPotential[j][k].overlap ? exp(-PartialPotential[j][k].potential*Beta) : 0.0;
      WeightBeadTotal[j] += WeightBeadTrials[j][k];
    }
    SelectedTrialOrientation = SelectTrialOrientation(WeightBeadTrials[j], WeightBeadTotal[j]);
    TestConfiguration.Molecules[MoleculeIndex].Atoms[j].Position = PositionCandidates[j][SelectedTrialOrientation];
    WeightChain *= WeightBeadTotal[j]/NUMBER_TRIAL_ORIENTATIONS;
  }
  if(!SimulationBox.ClosedBox){
    DeltaPotentialLongRangeCorrection = (
      GetPotentialLongRangeCorrection(TestConfiguration) - 
      GetPotentialLongRangeCorrection(Configuration)
    );
    WeightChain *= exp(-Beta*DeltaPotentialLongRangeCorrection);
  }

  return WeightChain;	
}

/* ***************************************************************************************************************
 * Name       | GenerateInitialConfiguration
 * ---------------------------------------------------------------------------------------------------------------
 * Function   | Generates an inition configuration
 * Parameters | Configuration (CONFIGURATION pointer) with 
 *              number of molecules and chain size already filled
 * Returns    | None
 * **************************************************************************************************************/
double GenerateInitialConfiguration(CONFIGURATION* Configuration){
  CONFIGURATION AuxConfiguration;
  POTENTIAL PartialPotential[MAX_CHAIN_SIZE][NUMBER_TRIAL_ORIENTATIONS];
  VECTOR PositionTrialOrientations[MAX_CHAIN_SIZE][NUMBER_TRIAL_ORIENTATIONS];
  VECTOR DirectionalVector;
  VECTOR BasePosition;
  double WeightBeadTrials[MAX_CHAIN_SIZE][NUMBER_TRIAL_ORIENTATIONS]={0};
  double WeightBeadTotal[MAX_CHAIN_SIZE]={0};
  double CellLength;
  double BondLength;
  double xSize, ySize, zSize;
  double xStart, xEnd, yStart, yEnd, zStart, zEnd;
  int SelectedTrialOrientation;

  xSize = SimulationBox.xSize - SigmaAlkane[2];
  ySize = SimulationBox.ySize - SigmaAlkane[2];
  zSize = SimulationBox.zSize - SigmaAlkane[2];

  CellLength = cbrt(xSize*ySize*zSize/Configuration->NumberMolecules);

  xStart = SimulationBox.xMin;
  yStart = SimulationBox.yMin;
  zStart = SimulationBox.zMin + SigmaAlkane[2]/2;

  xEnd = SimulationBox.xMax - SigmaAlkane[2];
  yEnd = SimulationBox.yMax - SigmaAlkane[2];
  zEnd = SimulationBox.zMax - SigmaAlkane[2]/2;

  BasePosition.x = xStart;
  BasePosition.y = yStart;
  BasePosition.z = zStart;

  AuxConfiguration.Molecules = (MOLECULE*) calloc(Configuration->NumberMolecules, sizeof(MOLECULE));
  for(int i=0; i<Configuration->NumberMolecules; i++){
    AuxConfiguration.NumberMolecules = i+1;
    AuxConfiguration.Molecules[i].Size = Configuration->Molecules[i].Size;
    AuxConfiguration.Molecules[i].MolarMass = Configuration->Molecules[i].MolarMass;
    AuxConfiguration.Molecules[i].Atoms = (ATOM*) calloc(AuxConfiguration.Molecules[i].Size, sizeof(ATOM));

    for(int j=0; j<AuxConfiguration.Molecules[i].Size; j++) 
      AuxConfiguration.Molecules[i].Atoms[j] = Configuration->Molecules[i].Atoms[j];

    AuxConfiguration.Molecules[i].Atoms[0].Position = BasePosition;

    for(int j=1; j<AuxConfiguration.Molecules[i].Size; j++){
      WeightBeadTotal[j]=0;
      for(int k=0; k<NUMBER_TRIAL_ORIENTATIONS; k++){
        PositionTrialOrientations[j][k] = SampleBeadPosition(AuxConfiguration.Molecules[i], j);
        AuxConfiguration.Molecules[i].Atoms[j].Position = PositionTrialOrientations[j][k];
        PartialPotential[j][k] = GetPartialExternalPotential(AuxConfiguration, i, j);
        WeightBeadTrials[j][k] = !PartialPotential[j][k].overlap ? exp(-Beta*PartialPotential[j][k].potential) : 0.0;
        WeightBeadTotal[j] += WeightBeadTrials[j][k];
      }
      SelectedTrialOrientation = SelectTrialOrientation(WeightBeadTrials[j], WeightBeadTotal[j]);
      AuxConfiguration.Molecules[i].Atoms[j].Position = PositionTrialOrientations[j][SelectedTrialOrientation];    
    }

    BasePosition.x += CellLength;
    if(BasePosition.x > xEnd){
      BasePosition.y += CellLength;
      BasePosition.x = xStart;
    }

    if(BasePosition.y > yEnd){
      BasePosition.z += CellLength;
      BasePosition.y = yStart;
    }
  }
  CopyConfiguration(AuxConfiguration, Configuration);
  return CellLength;
}

/* ***************************************************************************************************************
 * Name       | GetMoleculeDisplacement
 * ---------------------------------------------------------------------------------------------------------------
 * Function   | Generates displacement attempt
 * Parameters | - OldConfiguration (type CONFIGURATION)
 *              - NewConfiguration (type CONFIGURATION pointer)
 * Returns    | None
 * **************************************************************************************************************/
void GetMoleculeDisplacement(CONFIGURATION OldConfiguration, CONFIGURATION* NewConfiguration){
  POTENTIAL PartialPotential;
  POTENTIAL NewPartialPotential[MAX_CHAIN_SIZE][NUMBER_TRIAL_ORIENTATIONS];
  VECTOR NewCandidatePositions[MAX_CHAIN_SIZE][NUMBER_TRIAL_ORIENTATIONS]; 
  VECTOR OldChainPosition[MAX_CHAIN_SIZE];
  VECTOR DirectionalUnitVector;
  VECTOR NewBeadPosition;
  VECTOR Displacement;
  double NewWeightBeadTrials[MAX_CHAIN_SIZE][NUMBER_TRIAL_ORIENTATIONS] = {0}; 
  double NewAccumulatedWeight[MAX_CHAIN_SIZE][NUMBER_TRIAL_ORIENTATIONS] = {0};
  double OldWeightBeadTrials[MAX_CHAIN_SIZE][NUMBER_TRIAL_ORIENTATIONS] = {0}; 
  double NewWeightBeadTotal[MAX_CHAIN_SIZE] = {0}; 
  double OldWeightBeadTotal[MAX_CHAIN_SIZE] = {0}; 
  double NewWeightChain;
  double OldWeightChain;
  double BondLength;
  int ChosenMolecule = GetRandomIntegerInterval(0, OldConfiguration.NumberMolecules-1);
  int SelectedTrialOrientation;

  Displacement.x = GetRandomDoubleInterval(-1, 1)*MaxTranslationDistance;
  Displacement.y = GetRandomDoubleInterval(-1, 1)*MaxTranslationDistance;
  Displacement.z = GetRandomDoubleInterval(-1, 1)*MaxTranslationDistance;
  NewConfiguration->Molecules[ChosenMolecule].Atoms[0].Position = VectorSum(
    OldConfiguration.Molecules[ChosenMolecule].Atoms[0].Position, 
    Displacement
  );
  NewConfiguration->Molecules[ChosenMolecule].Atoms[0].Position = ApplyPeriodicBoundaryConditionsVector(
    NewConfiguration->Molecules[ChosenMolecule].Atoms[0].Position
  );
  PartialPotential = GetPartialExternalPotential(*NewConfiguration, ChosenMolecule, 0);
  NewWeightBeadTrials[0][0] = !PartialPotential.overlap ? exp(-PartialPotential.potential*Beta) : 0.0;
  NewWeightBeadTotal[0] = ((double) NUMBER_TRIAL_ORIENTATIONS)*NewWeightBeadTrials[0][0];
  NewWeightChain = NewWeightBeadTotal[0];

  for(int j=1; j<ChainSize; j++){
    for(int k=0; k<NUMBER_TRIAL_ORIENTATIONS; k++){
      NewCandidatePositions[j][k] = SampleBeadPosition(NewConfiguration->Molecules[ChosenMolecule], j);
      NewConfiguration->Molecules[ChosenMolecule].Atoms[j].Position = NewCandidatePositions[j][k];
      NewPartialPotential[j][k] = GetPartialExternalPotential(*NewConfiguration, ChosenMolecule, j);
      NewWeightBeadTrials[j][k] = (
        !NewPartialPotential[j][k].overlap ? exp(-NewPartialPotential[j][k].potential*Beta) : 0.0
      );
      NewWeightBeadTotal[j] += NewWeightBeadTrials[j][k];
    }
    SelectedTrialOrientation = SelectTrialOrientation(NewWeightBeadTrials[j], NewWeightBeadTotal[j]);
    NewConfiguration->Molecules[ChosenMolecule].Atoms[j].Position = NewCandidatePositions[j][SelectedTrialOrientation];
    NewWeightChain *= NewWeightBeadTotal[j];
  }

  for(int j = 0; j < ChainSize; j++)
    OldChainPosition[j] = OldConfiguration.Molecules[ChosenMolecule].Atoms[j].Position;

  PartialPotential = GetPartialExternalPotential(OldConfiguration, ChosenMolecule, 0);
  OldWeightBeadTrials[0][0] = !PartialPotential.overlap ? exp(-PartialPotential.potential*Beta) : 0.0;
  OldWeightBeadTotal[0] = ((double) NUMBER_TRIAL_ORIENTATIONS)*OldWeightBeadTrials[0][0];
  OldWeightChain = OldWeightBeadTotal[0];

  for(int j=1; j<ChainSize; j++){
    PartialPotential = GetPartialExternalPotential(OldConfiguration, ChosenMolecule, j);
    OldWeightBeadTrials[j][0] = !PartialPotential.overlap ? exp(-PartialPotential.potential*Beta) : 0.0;
    OldWeightBeadTotal[j] = OldWeightBeadTrials[j][0];
    for(int k=1; k<NUMBER_TRIAL_ORIENTATIONS; k++){
      OldConfiguration.Molecules[ChosenMolecule].Atoms[j].Position = SampleBeadPosition(
        OldConfiguration.Molecules[ChosenMolecule], 
        j
      );
      PartialPotential = GetPartialExternalPotential(OldConfiguration, ChosenMolecule, j);
      OldWeightBeadTrials[j][k] = !PartialPotential.overlap ? exp(-PartialPotential.potential*Beta) : 0.0;
      OldWeightBeadTotal[j] += OldWeightBeadTrials[j][k];
    }
    OldConfiguration.Molecules[ChosenMolecule].Atoms[j].Position = OldChainPosition[j];
    OldWeightChain *= OldWeightBeadTotal[j];
  }

  MovementAccepted = NewWeightChain > 0 ? Metropolis(NewWeightChain/OldWeightChain) : false;
  if(MovementAccepted) DisplacementStepsAccepted++;
}

/* ***************************************************************************************************************
 * Name       | GetMoleculeInsertion
 * ---------------------------------------------------------------------------------------------------------------
 * Function   | Generates an insertion attempt
 * Parameters | - OldConfiguration (type CONFIGURATION)
 *              - NewConfiguration (type CONFIGURATION pointer)
 * Returns    | None
 * **************************************************************************************************************/
void GetMoleculeInsertion(CONFIGURATION OldConfiguration, CONFIGURATION *NewConfiguration){
  POTENTIAL PartialPotential[MAX_CHAIN_SIZE][NUMBER_TRIAL_ORIENTATIONS];
  VECTOR DirectionalUnitVector;
  VECTOR PositionCandidates[MAX_CHAIN_SIZE][NUMBER_TRIAL_ORIENTATIONS];
  VECTOR* TrialPositionFirstBead;
  double WeightBeadTrials[MAX_CHAIN_SIZE][NUMBER_TRIAL_ORIENTATIONS] = {0};
  double WeightBeadTotal[NUMBER_TRIAL_ORIENTATIONS] = {0};
  double WeightChain;
  double BondLength;
  double VolumeMeters = SimulationBox.volume/Cube(METER_TO_ANGSTRON);
  double AcceptanceCriteria;
  double DeltaPotentialLongRangeCorrection;
  int SelectedTrialOrientation;
  int	MoleculeIndex;

  InsertionStepsAttempted++;

  NewConfiguration->NumberMolecules = OldConfiguration.NumberMolecules + 1;
  NewConfiguration->Molecules = (MOLECULE*) realloc(NewConfiguration->Molecules, NewConfiguration->NumberMolecules*sizeof(MOLECULE));
  for(int i = 0; i < OldConfiguration.NumberMolecules; i++){
    NewConfiguration->Molecules[i].Size = OldConfiguration.Molecules[i].Size;
    NewConfiguration->Molecules[i].MolarMass = OldConfiguration.Molecules[i].MolarMass;
    if(NewConfiguration->Molecules[i].Atoms==NULL){
      NewConfiguration->Molecules[i].Atoms = (ATOM*) calloc(NewConfiguration->Molecules[i].Size, sizeof(ATOM));
    }else{
      NewConfiguration->Molecules[i].Atoms = (ATOM*) realloc(NewConfiguration->Molecules[i].Atoms, NewConfiguration->Molecules[i].Size*sizeof(ATOM));
    }
    for(int j = 0; j < ChainSize; j++) 
      NewConfiguration->Molecules[i].Atoms[j] = OldConfiguration.Molecules[i].Atoms[j];
  }

  MoleculeIndex = NewConfiguration->NumberMolecules - 1;
  NewConfiguration->Molecules[MoleculeIndex].Size = OldConfiguration.Molecules[0].Size;
  NewConfiguration->Molecules[MoleculeIndex].MolarMass = OldConfiguration.Molecules[0].MolarMass;
  NewConfiguration->Molecules[MoleculeIndex].Atoms = (ATOM*) calloc(NewConfiguration->Molecules[MoleculeIndex].Size, sizeof(ATOM));

  for(int j = 0; j < NewConfiguration->Molecules[MoleculeIndex].Size; j++)
    NewConfiguration->Molecules[MoleculeIndex].Atoms[j] = OldConfiguration.Molecules[0].Atoms[j];
  
  NewConfiguration->Molecules[MoleculeIndex].Atoms[0].Position = GetRandomPosition();
  PartialPotential[0][0] = GetPartialExternalPotential(*NewConfiguration, MoleculeIndex, 0);
  WeightBeadTrials[0][0] = !PartialPotential[0][0].overlap ? exp(-PartialPotential[0][0].potential*Beta) : 0.0;
  WeightBeadTotal[0] = WeightBeadTrials[0][0];
  WeightChain = WeightBeadTotal[0];
  for(int j=1; j<ChainSize; j++){
    for(int k=0; k<NUMBER_TRIAL_ORIENTATIONS; k++){
      PositionCandidates[j][k] = SampleBeadPosition(NewConfiguration->Molecules[MoleculeIndex], j);			
      NewConfiguration->Molecules[MoleculeIndex].Atoms[j].Position = PositionCandidates[j][k];
      PartialPotential[j][k] = GetPartialExternalPotential(*NewConfiguration, MoleculeIndex, j);
      WeightBeadTrials[j][k] = !PartialPotential[j][k].overlap ? exp(-PartialPotential[j][k].potential*Beta) : 0.0;
      WeightBeadTotal[j] += WeightBeadTrials[j][k];
    }
    SelectedTrialOrientation = SelectTrialOrientation(WeightBeadTrials[j], WeightBeadTotal[j]);
    NewConfiguration->Molecules[MoleculeIndex].Atoms[j].Position = PositionCandidates[j][SelectedTrialOrientation];
    WeightChain *= WeightBeadTotal[j]/NUMBER_TRIAL_ORIENTATIONS;
  }

  if(!SimulationBox.ClosedBox)
    WeightChain *= exp(
      -Beta*(
        GetPotentialLongRangeCorrection(*NewConfiguration) 
        - GetPotentialLongRangeCorrection(OldConfiguration)
      )
    );
  
  AcceptanceCriteria = Beta*SimulationFugacity*VolumeMeters*WeightChain/(NewConfiguration->NumberMolecules*WeightIdealChain);
  
  MovementAccepted  = Metropolis(AcceptanceCriteria);
  if(MovementAccepted)
    InsertionStepsAccepted++;
}


/* ***************************************************************************************************************
 * Name       | GetMoleculeDeletion
 * ---------------------------------------------------------------------------------------------------------------
 * Function   | Computes a deleption attempt
 * Parameters | - OldConfiguration (type CONFIGURATION)
 *              - NewConfiguration (type CONFIGURATION pointer)
 * Returns    | None
 * **************************************************************************************************************/
void GetMoleculeDeletion(CONFIGURATION OldConfiguration, CONFIGURATION *NewConfiguration){
  POTENTIAL	PartialPotential;
  VECTOR DirectionalUnitVector;
  VECTOR PositionCandidates[MAX_CHAIN_SIZE][NUMBER_TRIAL_ORIENTATIONS];
  VECTOR* TrialPositionFirstBead;
  double WeightBeadTrials[MAX_CHAIN_SIZE][NUMBER_TRIAL_ORIENTATIONS];
  double WeightBeadTotal[NUMBER_TRIAL_ORIENTATIONS] = {0};
  double AccumulatedWeight[NUMBER_TRIAL_ORIENTATIONS] = {0};
  double WeightChain;
  double VolumeMeters = SimulationBox.volume*Cube(ANGSTRON);
  double AcceptanceCriteria;
  double DeltaPotentialLongRangeCorrection;
  double VolumeCavity;
  double ProbabilityCavity;
  bool* IsCavity;
  int SelectedTrialOrientation;
  int countTrials;
  int ChosenMolecule = GetRandomIntegerInterval(0, OldConfiguration.NumberMolecules - 1);
  int	CountCavity=0;

  DeletionStepsAttempted++;

  CopyConfiguration(OldConfiguration, NewConfiguration);

  PartialPotential = GetPartialExternalPotential(*NewConfiguration, ChosenMolecule, 0);
  WeightBeadTotal[0] = !PartialPotential.overlap ? exp(-Beta*PartialPotential.potential) : 0.0;
  WeightChain = WeightBeadTotal[0];

  for(int j=1; j<ChainSize; j++){
    PartialPotential = GetPartialExternalPotential(*NewConfiguration, ChosenMolecule, j);
    WeightBeadTotal[j] = !PartialPotential.overlap ? exp(-Beta*PartialPotential.potential) : 0.0;
    for(int k=1; k<NUMBER_TRIAL_ORIENTATIONS; k++){
      NewConfiguration->Molecules[ChosenMolecule].Atoms[j].Position = SampleBeadPosition(NewConfiguration->Molecules[ChosenMolecule], j);
      PartialPotential = GetPartialExternalPotential(*NewConfiguration, ChosenMolecule, j);
      WeightBeadTotal[j] += !PartialPotential.overlap ? exp(-Beta*PartialPotential.potential) : 0.0;
    }
    NewConfiguration->Molecules[ChosenMolecule].Atoms[j].Position = OldConfiguration.Molecules[ChosenMolecule].Atoms[j].Position;
    WeightChain *= WeightBeadTotal[j]/NUMBER_TRIAL_ORIENTATIONS;
  }

  if(!SimulationBox.ClosedBox)
    WeightChain *= exp(
      -Beta*(
        GetPotentialLongRangeCorrection(*NewConfiguration) 
        - GetPotentialLongRangeCorrection(OldConfiguration)
      )
    );

  if(WeightChain > 0){
    AcceptanceCriteria = OldConfiguration.NumberMolecules*WeightIdealChain/(Beta*VolumeMeters*SimulationFugacity*WeightChain);
    MovementAccepted = Metropolis(AcceptanceCriteria);
  }else{
    MovementAccepted = true;
  }
  
  if(MovementAccepted){
    int CountNew = 0;
    DeletionStepsAccepted++;
    NewConfiguration->NumberMolecules = OldConfiguration.NumberMolecules - 1;
    NewConfiguration->Molecules = (MOLECULE*) realloc(NewConfiguration->Molecules, NewConfiguration->NumberMolecules*sizeof(MOLECULE));
    for(int i = 0; i < OldConfiguration.NumberMolecules; i++){
      if(i != ChosenMolecule){
        NewConfiguration->Molecules[CountNew].Size = OldConfiguration.Molecules[i].Size;
        NewConfiguration->Molecules[CountNew].MolarMass = OldConfiguration.Molecules[i].MolarMass;
        if(NewConfiguration->Molecules[CountNew].Atoms==NULL){
          NewConfiguration->Molecules[CountNew].Atoms = (ATOM*) calloc(NewConfiguration->Molecules[CountNew].Size, sizeof(ATOM));
        }else{
          NewConfiguration->Molecules[CountNew].Atoms = (ATOM*) realloc(NewConfiguration->Molecules[CountNew].Atoms, NewConfiguration->Molecules[CountNew].Size*sizeof(ATOM));
        }

        for(int j = 0; j < NewConfiguration->Molecules[CountNew].Size; j++)
          NewConfiguration->Molecules[CountNew].Atoms[j] = OldConfiguration.Molecules[i].Atoms[j];
        CountNew++;
      }

    }
  }
}

/* ***************************************************************************************************************
 * Name       | GetVolumeChange
 * ---------------------------------------------------------------------------------------------------------------
 * Function   | Generates a volume change attempt
 * Parameters | - OldConfiguration (type CONFIGURATION)
 *              - NewConfiguration (type CONFIGURATION pointer)
 * Returns    | None
 * **************************************************************************************************************/
void GetVolumeChange(CONFIGURATION OldConfiguration, CONFIGURATION* NewConfiguration){
  BOX    OldBox;
  BOX	   NewBox;
  double PotentialOld, PotentialNew;
  double ScalingFactor;
  double DeltaH;
  double AcceptanceCriteria;
  VECTOR CenterOfMassDisplacement;
  VECTOR CenterOfMassNew;
  CONFIGURATION OldConfigurationCopy;
  VolumeStepsAttempted++;

  NewBox=OldBox=SimulationBox;

  PotentialOld = GetPotentialNonbonded(OldConfiguration, ReferencePotential)+GetPotentialBonded(OldConfiguration)+GetPotentialLongRangeCorrection(OldConfiguration);
  GetCenterOfMassAllMolecules(&OldConfiguration);

  NewBox.volume = SimulationBox.volume*exp(GetRandomDoubleInterval(-1, 1)*MaxLnVolumeChange);
  ScalingFactor = cbrt(NewBox.volume/OldBox.volume);
  NewBox.xSize = SimulationBox.xSize*ScalingFactor;
  NewBox.ySize = SimulationBox.ySize*ScalingFactor;
  NewBox.zSize = SimulationBox.zSize*ScalingFactor;
  NewBox.xMin = SimulationBox.xMin*ScalingFactor;
  NewBox.xMax = SimulationBox.xMax*ScalingFactor;
  NewBox.yMin = SimulationBox.yMin*ScalingFactor;
  NewBox.yMax = SimulationBox.yMax*ScalingFactor;
  NewBox.zMin = SimulationBox.zMin*ScalingFactor;
  NewBox.zMax = SimulationBox.zMax*ScalingFactor;

  SimulationBox = NewBox;

  CopyConfiguration(OldConfiguration, NewConfiguration);

  for(int i=0; i<NewConfiguration->NumberMolecules; i++){
    CenterOfMassNew = MultiplyVectorScalar(OldConfiguration.Molecules[i].CenterOfMass, ScalingFactor);
    CenterOfMassDisplacement = VectorSubtraction(CenterOfMassNew, OldConfiguration.Molecules[i].CenterOfMass);
    for(int j=0; j<NewConfiguration->Molecules[i].Size; j++){
      NewConfiguration->Molecules[i].Atoms[j].Position = VectorSum(NewConfiguration->Molecules[i].Atoms[j].Position, CenterOfMassDisplacement);
    }
  }

  PotentialNew = GetPotentialNonbonded(*NewConfiguration, ReferencePotential)+GetPotentialBonded(*NewConfiguration)+GetPotentialLongRangeCorrection(*NewConfiguration);
  DeltaH = ((PotentialNew - PotentialOld) + 
  SimulationPressure*(NewBox.volume - OldBox.volume)/Cube(METER_TO_ANGSTRON) - 
  ((NewConfiguration->NumberMolecules+1)/Beta)*log(NewBox.volume/OldBox.volume));
  AcceptanceCriteria = exp(-Beta*DeltaH);

  MovementAccepted = Metropolis(AcceptanceCriteria);

  if(MovementAccepted){
    VolumeStepsAccepted++;
  }else{
    SimulationBox = OldBox;
  }
}


/* ***************************************************************************************************************
 * Name       | GetNVTMovement
 * ---------------------------------------------------------------------------------------------------------------
 * Function   | Generates a NVT movement attempt
 * Parameters | - OldConfiguration (type CONFIGURATION)
 *              - NewConfiguration (type CONFIGURATION pointer)
 * Returns    | None
 * **************************************************************************************************************/
void GetNVTMovement(CONFIGURATION OldConfiguration, CONFIGURATION* NewConfiguration){
  DisplacementStepsAttempted++;
  CopyConfiguration(OldConfiguration, NewConfiguration);
  GetMoleculeDisplacement(OldConfiguration, NewConfiguration);
}

/* ***************************************************************************************************************
 * Name       | GetmuVTMovement
 * ---------------------------------------------------------------------------------------------------------------
 * Function   | Generates a muVT movement attempt
 * Parameters | - OldConfiguration (type CONFIGURATION)
 *              - NewConfiguration (type CONFIGURATION pointer)
 * Returns    | None
 * **************************************************************************************************************/
void GetmuVTMovement(CONFIGURATION OldConfiguration, CONFIGURATION* NewConfiguration){
  double RandomNumber = GetRandomNumber();
  double ProbabilitySum = DeletionAttemptProbability + InsertionAttemptProbability + DisplacementAttemptProbability;
  if(RandomNumber <= DeletionAttemptProbability/ProbabilitySum){
    GetMoleculeDeletion(OldConfiguration, NewConfiguration);
  }else if(RandomNumber < (DeletionAttemptProbability + InsertionAttemptProbability)/ProbabilitySum){
    GetMoleculeInsertion(OldConfiguration, NewConfiguration);
  }else{
    GetNVTMovement(OldConfiguration, NewConfiguration);
  }
}

/* ***************************************************************************************************************
 * Name       | GetNPTMovement
 * ---------------------------------------------------------------------------------------------------------------
 * Function   | Generates a NPT movement attempt
 * Parameters | - OldConfiguration (type CONFIGURATION)
 *              - NewConfiguration (type CONFIGURATION pointer)
 * Returns    | None
 * **************************************************************************************************************/
void GetNPTMovement(CONFIGURATION OldConfiguration, CONFIGURATION* NewConfiguration){
  if(GetRandomNumber() <= 1.0/((double) (OldConfiguration.NumberMolecules + 1))){
    GetVolumeChange(OldConfiguration, NewConfiguration);
  }else{
    GetNVTMovement(OldConfiguration, NewConfiguration);
  }
}

/* ***************************************************************************************************************
 * Name       | GetMonteCarloMove
 * ---------------------------------------------------------------------------------------------------------------
 * Function   | Generates a Monte Carlo movement attempt with the
 *              appropriate simulation ensemble
 * Parameters | - OldConfiguration (type CONFIGURATION)
 *              - NewConfiguration (type CONFIGURATION pointer)
 * Returns    | None
 * **************************************************************************************************************/
void GetMonteCarloMove(CONFIGURATION OldConfiguration, CONFIGURATION* NewConfiguration){
  MovementAccepted = false;
  switch (Ensemble){
    case NVT:
      GetNVTMovement(OldConfiguration, NewConfiguration);
      break;
    case NPT:
      GetNPTMovement(OldConfiguration, NewConfiguration);
      break;
    case muVT:
      GetmuVTMovement(OldConfiguration, NewConfiguration);
      break;
  }
}

// ---------------------------------------------------------------------------------------------------------------