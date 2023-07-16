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
void InitializeConfiguration(CONFIGURATION* Configuration) {
  Configuration->NumberMolecules = InitialNumberMolecules;

  // Calculate values for Molecules[0]
  MOLECULE moleculeZero;
  moleculeZero.Size = ChainSize;
  moleculeZero.Atoms = (ATOM*) calloc(ChainSize, sizeof(ATOM));
  moleculeZero.MolarMass = 0;

  for(int j = 0; j < ChainSize; j++) {
    if(ChainSize == 1) {
      moleculeZero.Atoms[j].Type = CH4;
    } else {
      moleculeZero.Atoms[j].Type = (
        (j == 0 || j == moleculeZero.Size - 1) ? 
        CH3 : 
        CH2);
    }
    moleculeZero.Atoms[j].Epsilon = GetAlkaneEpsilon(moleculeZero.Atoms[j].Type);
    moleculeZero.Atoms[j].Sigma = GetAlkaneSigma(moleculeZero.Atoms[j].Type);
    moleculeZero.Atoms[j].MolarMass = GetAlkaneAtomMolarMass(moleculeZero.Atoms[j].Type);
    moleculeZero.MolarMass += moleculeZero.Atoms[j].MolarMass;
  }

  for(int i = 0; i < InitialNumberMolecules; i++) {
    Configuration->Molecules[i].Size = moleculeZero.Size;
    Configuration->Molecules[i].Atoms = (ATOM*) calloc(ChainSize, sizeof(ATOM));
    Configuration->Molecules[i].MolarMass = moleculeZero.MolarMass;
    memcpy(Configuration->Molecules[i].Atoms, moleculeZero.Atoms, ChainSize * sizeof(ATOM));
  }
}

/* ***************************************************************************************************************
 * Name       | FreeConfiguration
 * ---------------------------------------------------------------------------------------------------------------
 * Function   | Frees the memory allocated for the configuration struct
 * Parameters | Configuration: a pointer to the configuration struct
 * Returns    | None
 * **************************************************************************************************************/
void FreeConfiguration(CONFIGURATION* Configuration){
  for(int i = 0; i < Configuration->NumberMolecules; i++)
  free(Configuration->Molecules[i].Atoms);
  free(Configuration->Molecules);
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
    BondLength = SampleBondLength();

    if (Bead == 1) {
        // Random unit vector for first bead
        DirectionalUnitVector = RandomUnitVector();
    } else if (Bead == 2) {
        // Calculate directional unit vector using bending angle
        DirectionalUnitVector = SampleBendingAngle(
          Molecule.Atoms[Bead-2].Position, 
          Molecule.Atoms[Bead-1].Position
        );
    } else if (Bead >= 3) {
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
  double AccumulatedWeight[NUMBER_TRIAL_ORIENTATIONS] = {0};
  double WeightChain = 1.0;
  int SelectedTrialOrientation;
  bool SimulationBoxOriginalState;

  SimulationBoxOriginalState = SimulationBox.ClosedBox;
  SimulationBox.ClosedBox = false;

  Configuration.NumberMolecules = 1;
  Configuration.Molecules = (MOLECULE*) calloc(Configuration.NumberMolecules, sizeof(MOLECULE));
  Configuration.Molecules[0].Size = Molecule.Size;
  Configuration.Molecules[0].MolarMass = Molecule.MolarMass;
  Configuration.Molecules[0].Atoms = (ATOM*) calloc(Configuration.Molecules[0].Size, sizeof(ATOM));

  for(int i = 0; i < Configuration.Molecules[0].Size; i++)
    Configuration.Molecules[0].Atoms[i] = Molecule.Atoms[i];

  for(int j = 0; j < Configuration.Molecules[0].Size; j++){
    WeightBeadTotal[j] = 0.0;
    if(j == 0){
      Configuration.Molecules[0].Atoms[j].Position = GetRandomPosition();
      PartialPotential[j][0] = GetPartialExternalPotential(Configuration, 0, j);
      WeightBeadTrials[j][0] = !PartialPotential[j][0].overlap ? exp(-Beta*PartialPotential[j][0].potential) : 0.0;
      WeightBeadTotal[j] = WeightBeadTrials[j][0];
      WeightChain *= WeightBeadTotal[j];
    }else{
      for(int k=0; k<NUMBER_TRIAL_ORIENTATIONS; k++){
        PositionCandidates[j][k] = SampleBeadPosition(Configuration.Molecules[0], j);
        Configuration.Molecules[0].Atoms[j].Position = PositionCandidates[j][k];
        PartialPotential[j][k] = GetPartialExternalPotential(Configuration, 0, j);
        WeightBeadTrials[j][k] = !PartialPotential[j][k].overlap ? exp(-Beta*PartialPotential[j][k].potential) : 0.0;
        WeightBeadTotal[j] += WeightBeadTrials[j][k];
        AccumulatedWeight[k] = WeightBeadTotal[j];
      }
      SelectedTrialOrientation = SelectTrialOrientation(WeightBeadTrials[j], WeightBeadTotal[j]);
      Configuration.Molecules[0].Atoms[j].Position = PositionCandidates[j][SelectedTrialOrientation];
      WeightChain *= WeightBeadTotal[j]/NUMBER_TRIAL_ORIENTATIONS;
    }
    free(Configuration.Molecules[0].Atoms);
  }

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
double GetRosenbluthWeightGhostMolecule(CONFIGURATION Configuration) {
    CONFIGURATION TestConfiguration;
    POTENTIAL Potential;
    POTENTIAL PartialPotential[MAX_CHAIN_SIZE][NUMBER_TRIAL_ORIENTATIONS];
    MOLECULE TestMolecule;
    VECTOR FirstBeadPosition;
    VECTOR DirectionalUnitVector;
    VECTOR PositionCandidates[MAX_CHAIN_SIZE][NUMBER_TRIAL_ORIENTATIONS];
    double WeightBeadTrials[MAX_CHAIN_SIZE][NUMBER_TRIAL_ORIENTATIONS] = {0};
    double AccumulatedWeight[NUMBER_TRIAL_ORIENTATIONS] = {0};
    double WeightChain = 0.0;
    double BondLength;
    double DeltaPotentialLongRangeCorrection;
    int SelectedTrialOrientation;
    int MoleculeIndex;

    // initializing test configuration
    TestConfiguration.NumberMolecules = Configuration.NumberMolecules + 1;
    TestConfiguration.Molecules = (MOLECULE*)malloc(TestConfiguration.NumberMolecules * sizeof(MOLECULE));
    for (int i = 0; i < Configuration.NumberMolecules; i++) {
        TestConfiguration.Molecules[i].Size = Configuration.Molecules[i].Size;
        TestConfiguration.Molecules[i].MolarMass = Configuration.Molecules[i].MolarMass;
        TestConfiguration.Molecules[i].Atoms = (ATOM*)calloc(TestConfiguration.Molecules[i].Size, sizeof(ATOM));
        memcpy(TestConfiguration.Molecules[i].Atoms, Configuration.Molecules[i].Atoms, TestConfiguration.Molecules[i].Size * sizeof(ATOM));
    }

    MoleculeIndex = TestConfiguration.NumberMolecules - 1;
    TestConfiguration.Molecules[MoleculeIndex].Size = Configuration.Molecules[0].Size;
    TestConfiguration.Molecules[MoleculeIndex].MolarMass = Configuration.Molecules[0].MolarMass;
    TestConfiguration.Molecules[MoleculeIndex].Atoms = (ATOM*)calloc(TestConfiguration.Molecules[0].Size, sizeof(ATOM));

    // insert molecule
    for (int j = 0; j < TestConfiguration.Molecules[MoleculeIndex].Size; j++) {
        if (j == 0) {
            TestConfiguration.Molecules[MoleculeIndex].Atoms[j].Position = GetRandomPosition();
            Potential = GetPartialExternalPotential(TestConfiguration, MoleculeIndex, j);
            WeightChain = !Potential.overlap ? exp(-Beta * Potential.potential) : 0.0;
        } else {
            for (int k = 0; k < NUMBER_TRIAL_ORIENTATIONS; k++) {
                PositionCandidates[j][k] = SampleBeadPosition(TestConfiguration.Molecules[MoleculeIndex], j);
                TestConfiguration.Molecules[MoleculeIndex].Atoms[j].Position = PositionCandidates[j][k];
                PartialPotential[j][k] = GetPartialExternalPotential(TestConfiguration, MoleculeIndex, j);
                WeightBeadTrials[j][k] = !PartialPotential[j][k].overlap ? exp(-PartialPotential[j][k].potential * Beta) : 0.0;
                AccumulatedWeight[k] = WeightBeadTrials[j][k] + (k > 0 ? AccumulatedWeight[k - 1] : 0.0);
            }
            SelectedTrialOrientation = SelectTrialOrientation(WeightBeadTrials[j], AccumulatedWeight[NUMBER_TRIAL_ORIENTATIONS - 1]);
            TestConfiguration.Molecules[MoleculeIndex].Atoms[j].Position = PositionCandidates[j][SelectedTrialOrientation];
            WeightChain *= AccumulatedWeight[NUMBER_TRIAL_ORIENTATIONS - 1] / NUMBER_TRIAL_ORIENTATIONS;
        }
    }

    if (!SimulationBox.ClosedBox) {
        DeltaPotentialLongRangeCorrection = GetPotentialLongRangeCorrection(TestConfiguration) - GetPotentialLongRangeCorrection(Configuration);
        WeightChain *= exp(-Beta * DeltaPotentialLongRangeCorrection);
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
void GenerateInitialConfiguration(CONFIGURATION* Configuration) {
  CONFIGURATION AuxConfiguration;
  POTENTIAL PartialPotential[MAX_CHAIN_SIZE][NUMBER_TRIAL_ORIENTATIONS];
  VECTOR PositionTrialOrientations[MAX_CHAIN_SIZE][NUMBER_TRIAL_ORIENTATIONS];
  VECTOR DirectionalVector;
  VECTOR BasePosition;
  double WeightBeadTrials[MAX_CHAIN_SIZE][NUMBER_TRIAL_ORIENTATIONS]={0};
  double WeightBeadTotal[MAX_CHAIN_SIZE]={0};
  double CellLength;
  double BondLength;
  double xStart, xEnd, yStart, yEnd, zStart, zEnd;

  CellLength = cbrt((SimulationBox.xSize * 0.9) * (SimulationBox.ySize * 0.9) * (SimulationBox.zSize * 0.9) / Configuration->NumberMolecules);

  xStart = SimulationBox.xMin + 0.05 * SimulationBox.xSize;
  yStart = SimulationBox.yMin + 0.05 * SimulationBox.ySize;
  zStart = SimulationBox.zMin + 0.05 * SimulationBox.zSize;

  xEnd = SimulationBox.xMax - 0.05 * SimulationBox.xSize;
  yEnd = SimulationBox.yMax - 0.05 * SimulationBox.ySize;
  zEnd = SimulationBox.zMax - 0.05 * SimulationBox.zSize;

  BasePosition.x = xStart;
  BasePosition.y = yStart;
  BasePosition.z = zStart;

  AuxConfiguration.Molecules = (MOLECULE*) malloc(Configuration->NumberMolecules * sizeof(MOLECULE));
  for (int i = 0; i < Configuration->NumberMolecules; i++) {
    int moleculeSize = Configuration->Molecules[i].Size;
    AuxConfiguration.NumberMolecules = i + 1;
    AuxConfiguration.Molecules[i].Size = moleculeSize;
    AuxConfiguration.Molecules[i].MolarMass = Configuration->Molecules[i].MolarMass;
    AuxConfiguration.Molecules[i].Atoms = (ATOM*) malloc(moleculeSize * sizeof(ATOM));

    memcpy(AuxConfiguration.Molecules[i].Atoms, Configuration->Molecules[i].Atoms, moleculeSize * sizeof(ATOM));

    for (int j = 0; j < moleculeSize; j++) {
      WeightBeadTotal[j] = 0;
      if (j == 0) {
        AuxConfiguration.Molecules[i].Atoms[j].Position = BasePosition;
      } else {
        for (int k = 0; k < NUMBER_TRIAL_ORIENTATIONS; k++) {
          PositionTrialOrientations[j][k] = SampleBeadPosition(AuxConfiguration.Molecules[i], j);
          AuxConfiguration.Molecules[i].Atoms[j].Position = PositionTrialOrientations[j][k];
          PartialPotential[j][k] = GetPartialExternalPotential(AuxConfiguration, i, j);
          WeightBeadTrials[j][k] = !PartialPotential[j][k].overlap ? exp(-Beta * PartialPotential[j][k].potential) : 0.0;
          WeightBeadTotal[j] += WeightBeadTrials[j][k];
        }
        AuxConfiguration.Molecules[i].Atoms[j].Position = PositionTrialOrientations[j][SelectTrialOrientation(WeightBeadTrials[j], WeightBeadTotal[j])];
      }
    }

    BasePosition.x += CellLength;
    if (BasePosition.x > xEnd) {
      BasePosition.y += CellLength;
      BasePosition.x = xStart;
    }

    if (BasePosition.y > yEnd) {
      BasePosition.z += CellLength;
      BasePosition.y = yStart;
    }
  }

  CopyConfiguration(AuxConfiguration, Configuration);

  // Free memory allocated for AuxConfiguration
  for (int i = 0; i < Configuration->NumberMolecules; i++) {
    free(AuxConfiguration.Molecules[i].Atoms);
  }
  free(AuxConfiguration.Molecules);
}

/* ***************************************************************************************************************
 * Name       | GetMoleculeDisplacement
 * ---------------------------------------------------------------------------------------------------------------
 * Function   | Generates displacement attempt
 * Parameters | - OldConfiguration (type CONFIGURATION)
 *              - NewConfiguration (type CONFIGURATION pointer)
 * Returns    | None
 * **************************************************************************************************************/
void GetMoleculeDisplacement(CONFIGURATION OldConfiguration, CONFIGURATION* NewConfiguration) {
  POTENTIAL PartialPotential;
  POTENTIAL NewPartialPotential[MAX_CHAIN_SIZE][NUMBER_TRIAL_ORIENTATIONS];
  VECTOR NewCandidatePositions[MAX_CHAIN_SIZE][NUMBER_TRIAL_ORIENTATIONS];
  VECTOR OldChainPosition[MAX_CHAIN_SIZE];
  VECTOR Displacement;
  double NewWeightBeadTrials[MAX_CHAIN_SIZE][NUMBER_TRIAL_ORIENTATIONS] = {0};
  double NewAccumulatedWeight[MAX_CHAIN_SIZE][NUMBER_TRIAL_ORIENTATIONS] = {0};
  double OldWeightBeadTrials[MAX_CHAIN_SIZE][NUMBER_TRIAL_ORIENTATIONS] = {0};
  double NewWeightBeadTotal[MAX_CHAIN_SIZE] = {0};
  double OldWeightBeadTotal[MAX_CHAIN_SIZE] = {0};
  double NewWeightChain;
  double OldWeightChain;
  double BondLength;
  int ChosenMolecule = GetRandomIntegerInterval(0, OldConfiguration.NumberMolecules - 1);
  int SelectedTrialOrientation;

  // Calculate Displacement outside the loop
  Displacement.x = GetRandomDoubleInterval(-1, 1) * MaxTranslationDistance;
  Displacement.y = GetRandomDoubleInterval(-1, 1) * MaxTranslationDistance;
  Displacement.z = GetRandomDoubleInterval(-1, 1) * MaxTranslationDistance;

  // Calculate PartialPotential and NewWeightBeadTrials for the first atom
  NewConfiguration->Molecules[ChosenMolecule].Atoms[0].Position = VectorSum(
    OldConfiguration.Molecules[ChosenMolecule].Atoms[0].Position,
    Displacement
  );
  VECTOR NewPosition = ApplyPeriodicBoundaryConditionsVector(NewConfiguration->Molecules[ChosenMolecule].Atoms[0].Position);
  NewConfiguration->Molecules[ChosenMolecule].Atoms[0].Position = NewPosition;
  PartialPotential = GetPartialExternalPotential(*NewConfiguration, ChosenMolecule, 0);
  NewWeightBeadTrials[0][0] = !PartialPotential.overlap ? exp(-PartialPotential.potential * Beta) : 0.0;
  NewWeightBeadTotal[0] = ((double)NUMBER_TRIAL_ORIENTATIONS) * NewWeightBeadTrials[0][0];
  NewWeightChain = NewWeightBeadTotal[0];

  // Calculate NewCandidatePositions, NewPartialPotential, and NewWeightBeadTrials for the remaining atoms
  for (int j=1; j<ChainSize; j++) {
    for (int k=0; k<NUMBER_TRIAL_ORIENTATIONS; k++) {
      NewCandidatePositions[j][k] = SampleBeadPosition(NewConfiguration->Molecules[ChosenMolecule], j);
      NewConfiguration->Molecules[ChosenMolecule].Atoms[j].Position = NewCandidatePositions[j][k];
      NewPartialPotential[j][k] = GetPartialExternalPotential(*NewConfiguration, ChosenMolecule, j);
      NewWeightBeadTrials[j][k] = (
        !NewPartialPotential[j][k].overlap ? exp(-NewPartialPotential[j][k].potential * Beta) : 0.0
      );
    }

    // Calculate NewWeightBeadTotal and SelectedTrialOrientation
    double weightSum = 0.0;
    for (int k = 0; k < NUMBER_TRIAL_ORIENTATIONS; k++) {
      weightSum += NewWeightBeadTrials[j][k];
    }
    NewWeightBeadTotal[j] = weightSum;
    SelectedTrialOrientation = SelectTrialOrientation(NewWeightBeadTrials[j], NewWeightBeadTotal[j]);

    // Update NewConfiguration and NewWeightChain
    NewConfiguration->Molecules[ChosenMolecule].Atoms[j].Position = NewCandidatePositions[j][SelectedTrialOrientation];
    NewWeightChain *= NewWeightBeadTotal[j];
  }

  // Store OldChainPosition before modifying it
  for (int j = 0; j < ChainSize; j++) {
    OldChainPosition[j] = OldConfiguration.Molecules[ChosenMolecule].Atoms[j].Position;
  }

  // Calculate OldWeightBeadTrials and OldWeightBeadTotal for each atom
  for (int j = 0; j < ChainSize; j++) {
    PartialPotential = GetPartialExternalPotential(OldConfiguration, ChosenMolecule, j);
    OldWeightBeadTrials[j][0] = !PartialPotential.overlap ? exp(-PartialPotential.potential * Beta) : 0.0;
    OldWeightBeadTotal[j] = OldWeightBeadTrials[j][0];

    if (j > 0) {
      OldConfiguration.Molecules[ChosenMolecule].Atoms[j].Position = SampleBeadPosition(
        OldConfiguration.Molecules[ChosenMolecule],
        j
      );
      PartialPotential = GetPartialExternalPotential(OldConfiguration, ChosenMolecule, j);
      OldWeightBeadTrials[j][1] = !PartialPotential.overlap ? exp(-PartialPotential.potential * Beta) : 0.0;
      OldWeightBeadTotal[j] += OldWeightBeadTrials[j][1];

      for (int k = 2; k < NUMBER_TRIAL_ORIENTATIONS; k++) {
        OldConfiguration.Molecules[ChosenMolecule].Atoms[j].Position = SampleBeadPosition(
          OldConfiguration.Molecules[ChosenMolecule],
          j
        );
        PartialPotential = GetPartialExternalPotential(OldConfiguration, ChosenMolecule, j);
        OldWeightBeadTrials[j][k] = !PartialPotential.overlap ? exp(-PartialPotential.potential * Beta) : 0.0;
        OldWeightBeadTotal[j] += OldWeightBeadTrials[j][k];
      }
      OldConfiguration.Molecules[ChosenMolecule].Atoms[j].Position = OldChainPosition[j];
      OldWeightChain *= OldWeightBeadTotal[j];
    }else{
      OldWeightChain = OldWeightBeadTotal[j];
    }

  }

  MovementAccepted = NewWeightChain > 0 ? Metropolis(NewWeightChain / OldWeightChain) : false;
  if (MovementAccepted)
    DisplacementStepsAccepted++;
}

/* ***************************************************************************************************************
 * Name       | GetMoleculeInsertion
 * ---------------------------------------------------------------------------------------------------------------
 * Function   | Generates an insertion attempt
 * Parameters | - OldConfiguration (type CONFIGURATION)
 *              - NewConfiguration (type CONFIGURATION pointer)
 * Returns    | None
 * **************************************************************************************************************/
void GetMoleculeInsertion(CONFIGURATION OldConfiguration, CONFIGURATION *NewConfiguration) {
  POTENTIAL PartialPotential[MAX_CHAIN_SIZE * NUMBER_TRIAL_ORIENTATIONS];
  VECTOR DirectionalUnitVector;
  VECTOR PositionCandidates[MAX_CHAIN_SIZE * NUMBER_TRIAL_ORIENTATIONS];
  VECTOR* TrialPositionFirstBead;
  double WeightBeadTrials[MAX_CHAIN_SIZE * NUMBER_TRIAL_ORIENTATIONS] = {0};
  double WeightBeadTotal[NUMBER_TRIAL_ORIENTATIONS] = {0};
  double WeightChain;
  double BondLength;
  double VolumeMeters = SimulationBox.volume / Cube(METER_TO_ANGSTRON);
  double VolumeCavities;
  double AcceptanceCriteria;
  double DeltaPotentialLongRangeCorrection;
  double ProbabilityCavity;
  bool* IsCavity;
  int SelectedTrialOrientation;
  int MoleculeIndex;
  int CountCavity = 0;
  int SelectedCavity;
  int CavityList[NUMBER_TRIAL_CAVITY_SEARCH];

  InsertionStepsAttempted++;

  NewConfiguration->NumberMolecules = OldConfiguration.NumberMolecules + 1;
  NewConfiguration->Molecules = realloc(NewConfiguration->Molecules, NewConfiguration->NumberMolecules * sizeof(MOLECULE));
  memcpy(
    NewConfiguration->Molecules, 
    OldConfiguration.Molecules, 
    OldConfiguration.NumberMolecules * sizeof(MOLECULE)
  );

  MoleculeIndex = NewConfiguration->NumberMolecules - 1;
  NewConfiguration->Molecules[MoleculeIndex].Size = OldConfiguration.Molecules[0].Size;
  NewConfiguration->Molecules[MoleculeIndex].MolarMass = OldConfiguration.Molecules[0].MolarMass;
  NewConfiguration->Molecules[MoleculeIndex].Atoms = calloc(NewConfiguration->Molecules[MoleculeIndex].Size, sizeof(ATOM));
  memcpy(
    NewConfiguration->Molecules[MoleculeIndex].Atoms, 
    OldConfiguration.Molecules[0].Atoms, 
    NewConfiguration->Molecules[MoleculeIndex].Size * sizeof(ATOM)
  );

  for (int j = 0; j < NewConfiguration->Molecules[MoleculeIndex].Size; j++)
    NewConfiguration->Molecules[MoleculeIndex].Atoms[j] = OldConfiguration.Molecules[0].Atoms[j];

  CountCavity = 0;

  NewConfiguration->Molecules[MoleculeIndex].Atoms[0].Position = GetRandomPosition();

  WeightChain = 0.0;
  for (int j = 0; j < ChainSize; j++) {
    if (j == 0) {
      PartialPotential[j] = GetPartialExternalPotential(*NewConfiguration, MoleculeIndex, j);
      WeightBeadTrials[j] = !PartialPotential[j].overlap ? exp(-PartialPotential[j].potential * Beta) : 0.0;
      WeightBeadTotal[j] = WeightBeadTrials[j];
      WeightChain = WeightBeadTotal[j];
    } else {
      for (int k = 0; k < NUMBER_TRIAL_ORIENTATIONS; k++) {
        int index = j * NUMBER_TRIAL_ORIENTATIONS + k;
        PositionCandidates[index] = SampleBeadPosition(NewConfiguration->Molecules[MoleculeIndex], j);
        NewConfiguration->Molecules[MoleculeIndex].Atoms[j].Position = PositionCandidates[index];
        PartialPotential[index] = GetPartialExternalPotential(*NewConfiguration, MoleculeIndex, j);
        WeightBeadTrials[index] = !PartialPotential[index].overlap ? exp(-PartialPotential[index].potential * Beta) : 0.0;
        WeightBeadTotal[j] += WeightBeadTrials[index];
      }
      SelectedTrialOrientation = SelectTrialOrientation(&WeightBeadTrials[j * NUMBER_TRIAL_ORIENTATIONS], WeightBeadTotal[j]);
      int selectedIndex = j * NUMBER_TRIAL_ORIENTATIONS + SelectedTrialOrientation;
      NewConfiguration->Molecules[MoleculeIndex].Atoms[j].Position = PositionCandidates[selectedIndex];
      WeightChain *= WeightBeadTotal[j] / NUMBER_TRIAL_ORIENTATIONS;
    }
  }

  if (!SimulationBox.ClosedBox)
    WeightChain *= exp(
      -Beta * (
        GetPotentialLongRangeCorrection(*NewConfiguration)
        - GetPotentialLongRangeCorrection(OldConfiguration)
      )
    );

  AcceptanceCriteria = Beta * SimulationFugacity * VolumeMeters * WeightChain / (NewConfiguration->NumberMolecules * WeightIdealChain);

  MovementAccepted = Metropolis(AcceptanceCriteria);
  if (MovementAccepted)
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
void GetMoleculeDeletion(CONFIGURATION OldConfiguration, CONFIGURATION *NewConfiguration) {
  POTENTIAL PartialPotential;
  VECTOR DirectionalUnitVector;
  VECTOR PositionCandidates[MAX_CHAIN_SIZE][NUMBER_TRIAL_ORIENTATIONS];
  VECTOR* TrialPositionFirstBead;
  double WeightBeadTrials[MAX_CHAIN_SIZE][NUMBER_TRIAL_ORIENTATIONS];
  double WeightBeadTotal[NUMBER_TRIAL_ORIENTATIONS] = {0};
  double WeightChain;
  double BondLength;
  double VolumeMeters = SimulationBox.volume*Cube(ANGSTRON);
  double AcceptanceCriteria;
  double DeltaPotentialLongRangeCorrection;
  bool* IsCavity;
  int SelectedTrialOrientation;
  int countTrials;
  int ChosenMolecule = GetRandomIntegerInterval(0, OldConfiguration.NumberMolecules - 1);
  int	CountCavity=0;

  DeletionStepsAttempted++;

  CopyConfiguration(OldConfiguration, NewConfiguration);

  for(int j=0; j<ChainSize; j++){
    WeightBeadTotal[j] = 0.0;
    PartialPotential = GetPartialExternalPotential(*NewConfiguration, ChosenMolecule, j);
    if(j == 0){
      WeightBeadTotal[j] = !PartialPotential.overlap ? exp(-Beta*PartialPotential.potential) : 0.0;
      WeightChain = WeightBeadTotal[j];
    }else{
      WeightBeadTotal[j] = !PartialPotential.overlap ? exp(-Beta*PartialPotential.potential) : 0.0;
      for(int k=1; k<NUMBER_TRIAL_ORIENTATIONS; k++){
        NewConfiguration->Molecules[ChosenMolecule].Atoms[j].Position = SampleBeadPosition(NewConfiguration->Molecules[ChosenMolecule], j);
        PartialPotential = GetPartialExternalPotential(*NewConfiguration, ChosenMolecule, j);
        WeightBeadTotal[j] += !PartialPotential.overlap ? exp(-Beta*PartialPotential.potential) : 0.0;
      }
      NewConfiguration->Molecules[ChosenMolecule].Atoms[j].Position = OldConfiguration.Molecules[ChosenMolecule].Atoms[j].Position;
      WeightChain *= WeightBeadTotal[j]/NUMBER_TRIAL_ORIENTATIONS;
    }
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
    for(int i=0; i<OldConfiguration.NumberMolecules; i++){
      if(i!=ChosenMolecule){
        NewConfiguration->Molecules[CountNew].Size = OldConfiguration.Molecules[i].Size;
        NewConfiguration->Molecules[CountNew].MolarMass = OldConfiguration.Molecules[i].MolarMass;
        if(NewConfiguration->Molecules[CountNew].Atoms == NULL){
          NewConfiguration->Molecules[CountNew].Atoms = (ATOM*) calloc(NewConfiguration->Molecules[CountNew].Size, sizeof(ATOM));
        }else{
          NewConfiguration->Molecules[CountNew].Atoms = (ATOM*) realloc(
            NewConfiguration->Molecules[CountNew].Atoms, 
            NewConfiguration->Molecules[CountNew].Size*sizeof(ATOM)
          );
        }

        for (int j = 0; j < NewConfiguration->Molecules[CountNew].Size; j++)
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
void GetVolumeChange(CONFIGURATION OldConfiguration, CONFIGURATION* NewConfiguration) {
  BOX OldBox;
  BOX NewBox;
  double PotentialOld, PotentialNew;
  double ScalingFactor;
  double DeltaH;
  double AcceptanceCriteria;
  VECTOR CenterOfMassDisplacement;
  VECTOR CenterOfMassNew;

  NewBox = OldBox = SimulationBox;

  PotentialOld = GetPotentialNonbonded(OldConfiguration, ReferencePotential)
                  + GetPotentialBonded(OldConfiguration)
                  + GetPotentialLongRangeCorrection(OldConfiguration);

  GetCenterOfMassAllMolecules(&OldConfiguration);

  double randomInterval = GetRandomDoubleInterval(-1, 1);
  double lnVolumeChange = randomInterval * MaxLnVolumeChange;
  NewBox.volume = SimulationBox.volume * exp(lnVolumeChange);

  ScalingFactor = cbrt(NewBox.volume / OldBox.volume);
  double scaledSize = SimulationBox.xSize * ScalingFactor;

  NewBox.xSize = NewBox.ySize = NewBox.zSize = scaledSize;
  NewBox.xMin = NewBox.xMax = SimulationBox.xMin * ScalingFactor;
  NewBox.yMin = NewBox.yMax = SimulationBox.yMin * ScalingFactor;
  NewBox.zMin = NewBox.zMax = SimulationBox.zMin * ScalingFactor;

  SimulationBox = NewBox;

  *NewConfiguration = OldConfiguration;

  for (int i = 0; i < NewConfiguration->NumberMolecules; i++) {
    CenterOfMassNew = MultiplyVectorScalar(OldConfiguration.Molecules[i].CenterOfMass, ScalingFactor);
    CenterOfMassDisplacement = VectorSubtraction(CenterOfMassNew, OldConfiguration.Molecules[i].CenterOfMass);

    MOLECULE* molecule = &NewConfiguration->Molecules[i];
    for (int j = 0; j < molecule->Size; j++) {
      molecule->Atoms[j].Position = VectorSum(molecule->Atoms[j].Position, CenterOfMassDisplacement);
    }
  }

  PotentialNew = GetPotentialNonbonded(*NewConfiguration, ReferencePotential)
                  + GetPotentialBonded(*NewConfiguration)
                  + GetPotentialLongRangeCorrection(*NewConfiguration);

  DeltaH = ((PotentialNew - PotentialOld)
            + SimulationPressure * (NewBox.volume - OldBox.volume) / Cube(METER_TO_ANGSTRON)
            - ((NewConfiguration->NumberMolecules + 1) / Beta) * log(NewBox.volume / OldBox.volume));

  AcceptanceCriteria = exp(-Beta * DeltaH);

  MovementAccepted = Metropolis(AcceptanceCriteria);

  if (MovementAccepted) {
    VolumeStepsAccepted++;
  } else {
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