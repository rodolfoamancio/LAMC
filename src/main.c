/* ***************************************************************************************************************
 * This program was developed to perform a configurational-bias monte carlo simulation for a pure component
 * alkene system in the NVT ensemble. The chain is treated in terms of pseudo-atoms (CH4, CH3, CH2) with TraPPE 
 * force-field parameters.
 * 
 * Developer: Rodolfo J Amancio (rodolfojamancio@gmail.com)
 * University of Campinas, School of Chemical Engineering, Campinas/SP -  Brazil
 * **************************************************************************************************************/

// ---------------------------------------------------------------------------------------------------------------
// Libraries loading
// ---------------------------------------------------------------------------------------------------------------

#include<stdio.h>
#include<stdbool.h>
#include<math.h>
#include"constants.h"
#include"vectors.h"
#include"molecule.h"
#include"potential.h"
#include"monte_carlo.h"
#include"properties.h"
#include"structural.h"
#include"files_io.h"
#include"simulation_setup.h"

double RunSimulation(char* InputFileName){
  FILE *PropertiesDataFile, *ConfigurationsFile, *ProfileFile;
  double PressureExcess, PressureTotal, PressureLongRangeCorrection, PressureIdealGas;
  double PotentialNonbonded, PotentialBonded, PotentialWalls, PotentialLongRangeCorrection, PotentialTotal, PotentialPerturbed;
  double WeightGhostMolecule;
  double DensityMolar, DensityMass;
  double WeightIdealChainSample, WeightIdealChainSquared, STDWeightIdealChain;  
  double CellLength;
  int CountMovements = 0;
  int ProductionFlag = 0;
  int OptimizeAcceptanceEverySteps;
  CONFIGURATION OldConfiguration, NewConfiguration;

  // Variables initialization
  StartTime = time(&StartTime);
  CountProductionStates=0;
  MaxTranslationDistance = 0.5;
  MaxLnVolumeChange = 0.1;
  DisplacementStepsAttempted=InsertionStepsAttempted=DeletionStepsAttempted=VolumeStepsAttempted=0;
  DisplacementStepsAccepted=InsertionStepsAccepted=DeletionStepsAccepted=VolumeStepsAccepted=0;
  VolumeAcceptanceCalculated=DisplacementAcceptanceCalculated=0;

  // Placeholder parameters before reading input files
  ReferencePotential = MIE;
  PerturbationPotential = MIE;

  // Input file data
  ReadInputFile(InputFileName);

  Beta = 1/(BOLTZMANN_CONSTANT*Temperature);
  SimulationBox.volume = SimulationBox.xSize*SimulationBox.ySize*SimulationBox.zSize;
  Ensemble = EquilibraionEnsemble;
  OptimizeAcceptanceEverySteps = (int) InitialNumberMolecules;
  NumberTotalSteps = (int) NumberTotalCycles*InitialNumberMolecules;
  NumberEquilibrationSteps = (int) NumberEquilibrationCycles*InitialNumberMolecules;
  StepsCalculateProperties = (int) CyclesCalculateProperties*InitialNumberMolecules;
  StepsRecordConfiguration = (int) CyclesRecordConfiguration*InitialNumberMolecules;
  StepsGetProfiles = (int) CyclesGetProfiles*InitialNumberMolecules;

  printf(
    "-----------------------------------------------------Inputs-------------------------------------------\n"
    "Reference potential:.......................: %s\n"
    "Perturbation potential:....................: %s\n"
    "Equilibration ensemble:....................: %s \n"
    "Simulation ensemble:.......................: %s \n"
    "Temperature:...............................: %f K\n"
    "Number total cycles:.......................: %d \n"
    "Number equilibration cycles:...............: %d \n"
    "Number total steps:........................: %d \n"
    "Number equilibration steps:................: %d \n"
    "Acceptance ratio:..........................: %f \n"
    "Cycles to calculate properties:............: %d \n"
    "Cycles to record configuration:............: %d \n"
    "Cycles to get profiles:....................: %d \n"
    "Steps to calculate properties:.............: %d \n"
    "Steps to record configuration:.............: %d \n"
    "Steps to get profiles:.....................: %d \n"
    "Simulation box:............................: %s \n"
    "Displacement probability:..................: %f \n"
    "Inserion probability:......................: %f \n"
    "Deletion probability:......................: %f \n"
    "Number molecules:..........................: %d \n"
    "Chain size:................................: %d \n"
    "Fugacity:..................................: %f MPa\n"
    "Pressure:..................................: %f MPa\n"
    "X size.....................................: %f \n"
    "Y size.....................................: %f \n"
    "Z size.....................................: %f \n"
    "------------------------------------------------------------------------------------------------------\n\n",
    GetPotentialTypeLabel(ReferencePotential), GetPotentialTypeLabel(PerturbationPotential), GetEnsembleLabel(EquilibraionEnsemble),
    GetEnsembleLabel(SimulationEnsemble), Temperature, NumberTotalCycles, NumberEquilibrationCycles, NumberTotalSteps,
    NumberEquilibrationSteps, DisplacementAcceptanceTarget, CyclesCalculateProperties, CyclesRecordConfiguration, 
    CyclesGetProfiles, StepsCalculateProperties, StepsRecordConfiguration, StepsGetProfiles, SimulationBox.ClosedBox ? "Closed" : "Open",
    DisplacementAttemptProbability, InsertionAttemptProbability, DeletionAttemptProbability, InitialNumberMolecules, ChainSize,
    SimulationFugacity/1E6, SimulationPressure/1E6, SimulationBox.xSize, SimulationBox.ySize, SimulationBox.zSize
  );

  InitializeConfiguration(&OldConfiguration);
  InitializeConfiguration(&NewConfiguration);

  printf("--------------------------------------------Ideal chain calculations----------------------------------\n");
  WeightIdealChain=WeightIdealChainSquared=0.0;
  for(int i=0; i<NUMBER_SAMPLES_IDEAL_CHAIN; i++){
    WeightIdealChainSample = GetRosenbluthWeightIdealChain(OldConfiguration.Molecules[0]);
    WeightIdealChain += WeightIdealChainSample;
    WeightIdealChainSquared += Squared(WeightIdealChainSample);
  }
  WeightIdealChain /= NUMBER_SAMPLES_IDEAL_CHAIN;
  WeightIdealChainSquared /= NUMBER_SAMPLES_IDEAL_CHAIN;
  STDWeightIdealChain = sqrt(WeightIdealChainSquared - WeightIdealChain);

  printf(
    "Weight ideal chain:........................: %f\n"
    "STD Weight ideal chain:....................: %f\n"
    "------------------------------------------------------------------------------------------------------\n\n", 
    WeightIdealChain, STDWeightIdealChain
  );

  printf("\n----------------------------------------Generating initial configuration------------------------------\n\n");
  CellLength = GenerateInitialConfiguration(&OldConfiguration);
  printf("Initial configuration generated\n");
  if(ReferencePotential == HARD_SPHERE && CellLength <= OldConfiguration.Molecules[0].Atoms[0].Sigma){
    printf("Cell length (%.2f) is to small for Hard-sphere simulation. Reduce system density!. Simulation aborted!", CellLength);
    return -1;
  }
  // Initialize output files
  PropertiesDataFile = InitializePropertiesDataFile(InputFileName);
  ConfigurationsFile = InitializeConfigurationFile(InputFileName);
  
  if(SimulationBox.ClosedBox) ProfileFile = InitializeProfilesFile(InputFileName);
  
  RecordConfiguration(ConfigurationsFile, OldConfiguration);

  // Initialize Statistics variables

  AveragePressure=AveragePressureExcess=AveragePressureIdealGas=AveragePressureLongRangeCorrection=0;
  AveragePotential=AveragePotentialBonded=AveragePotentialNonbonded=AveragePotentialLongRangeCorrection=AveragePotentialWalls=0;
  AverageNumberMolecules=AverageDensityMass=AverageDensityMolar=0;

  DensityMass = GetDensityMass(OldConfiguration);
  DensityMolar = GetDensityMolar(OldConfiguration);

  printf("\n-------------------------------------------Simulation starting----------------------------------------\n\n");
  
  // --------------------------------------------------------------------------------------------------------------
  // Monte carlo loop
  // --------------------------------------------------------------------------------------------------------------

  CopyConfiguration(OldConfiguration, &NewConfiguration);
  
  for(int i=1; i<=NumberTotalSteps + 1; i++){
    if(i>=NumberEquilibrationSteps/2) Ensemble = SimulationEnsemble;
    if(i>NumberEquilibrationSteps) ProductionFlag = 1;

    GetMonteCarloMove(OldConfiguration, &NewConfiguration);
    
    if(MovementAccepted) CopyConfiguration(NewConfiguration, &OldConfiguration);
    
    if(StepsGetProfiles>0){
      if((i%StepsGetProfiles==0 && SimulationBox.ClosedBox)){
        GetCenterOfMassAllMolecules(&OldConfiguration);
        CalculateOrderParameter(&OldConfiguration);
        RecordProfileData(ProfileFile, OldConfiguration);
      }
    }

    
    if(i%StepsCalculateProperties==0){
      if(Ensemble!=NVT){
        DensityMass = GetDensityMass(OldConfiguration);
        DensityMolar = GetDensityMolar(OldConfiguration);
      }
      WeightGhostMolecule = (Ensemble != muVT) ? GetRosenbluthWeightGhostMolecule(OldConfiguration) : 0.0;
      PressureIdealGas = GetPressureIdealGas(OldConfiguration.NumberMolecules, SimulationBox.volume);
      if(!SimulationBox.ClosedBox){
        if(Ensemble!=NPT){
          PressureLongRangeCorrection = GetPressureLongRangeCorrection(OldConfiguration);
          PressureExcess = GetPressureExcess(OldConfiguration);
          PressureTotal = PressureIdealGas + PressureLongRangeCorrection + PressureExcess;
        }else{
          PressureLongRangeCorrection = -1;
          PressureExcess = -1;
          PressureTotal = SimulationPressure;
        }
      }
      
      PotentialNonbonded = GetPotentialNonbonded(OldConfiguration, ReferencePotential);
      PotentialBonded = GetPotentialBonded(OldConfiguration);
      PotentialLongRangeCorrection = GetPotentialLongRangeCorrection(OldConfiguration, ReferencePotential);
      PotentialWalls = GetTotalPotentialExternal(OldConfiguration);
      
      PotentialPerturbed = (
        GetPotentialNonbonded(OldConfiguration, PerturbationPotential)
        + GetPotentialLongRangeCorrection(OldConfiguration, PerturbationPotential)
      );
      PotentialTotal = PotentialNonbonded + PotentialBonded + PotentialLongRangeCorrection + PotentialWalls;

      fprintf(
        PropertiesDataFile,
        "%d %s %d %f %f %f %f %f %f %d %f %f %e %e %e %e %e %e %e %e %e %e %e %e %e\n",
        i, GetEnsembleLabel(Ensemble), ProductionFlag, DisplacementAcceptanceCalculated,
        MaxTranslationDistance, Temperature, SimulationBox.xSize, SimulationBox.ySize,
        SimulationBox.zSize, OldConfiguration.NumberMolecules, DensityMass, DensityMolar,
        PotentialTotal, PotentialBonded, PotentialNonbonded, PotentialLongRangeCorrection, 
        PotentialWalls, PotentialPerturbed, WeightIdealChain, STDWeightIdealChain, 
        WeightGhostMolecule, PressureTotal, PressureExcess, PressureIdealGas, PressureLongRangeCorrection
      );

      if(ProductionFlag==1){

        CountProductionStates++;

        AveragePressure += PressureTotal;
        AveragePressureExcess += PressureExcess;
        AveragePressureIdealGas += PressureIdealGas;
        AveragePressureLongRangeCorrection += PressureLongRangeCorrection;

        AveragePotential += PotentialTotal;
        AveragePotentialBonded += PotentialBonded;
        AveragePotentialNonbonded += PotentialNonbonded;
        AveragePotentialLongRangeCorrection += PotentialLongRangeCorrection;
        AveragePotentialWalls += PotentialWalls;

        AverageWeightGhostMolecule += WeightGhostMolecule;

        AverageNumberMolecules += OldConfiguration.NumberMolecules;
        AverageDensityMass += DensityMass;
        AverageDensityMolar += DensityMolar;

      }
      
    }
        
    // Record configuration
    if(i % StepsRecordConfiguration == 0) RecordConfiguration(ConfigurationsFile, OldConfiguration);
    
    // adjust maximum displacement according to the acceptance ratio
    DisplacementAcceptanceCalculated = (double) DisplacementStepsAccepted/DisplacementStepsAttempted;
    if(i%OptimizeAcceptanceEverySteps==0 && i<NumberEquilibrationSteps){
      MaxTranslationDistance *= LimMinMax(DisplacementAcceptanceCalculated/DisplacementAcceptanceTarget, 0.1, 1.9);
      MaxTranslationDistance = LimMinMax(MaxTranslationDistance, 0.01, 20.0);
    }

    // adjust maximum volume change according to the acceptance ratio
    VolumeAcceptanceCalculated = (double) VolumeStepsAccepted/VolumeStepsAttempted;
    if(Ensemble==NPT && VolumeStepsAttempted%10==0 && i<NumberEquilibrationSteps/2){
      MaxLnVolumeChange *= 1 + LimMinMax((VolumeAcceptanceCalculated - VolumeAcceptanceTarget), -0.05, 0.05);
      MaxLnVolumeChange = LimMinMax(MaxLnVolumeChange, 0.01, 0.7);
    }
  }
  
  RecordSimulationLog(InputFileName);
  RecordConfiguration(ConfigurationsFile, OldConfiguration);
  if(SimulationBox.ClosedBox) fclose(ProfileFile);
  fclose(ConfigurationsFile);
  fclose(PropertiesDataFile);

  printf("\n-----------------------------------------Simulation has finished--------------------------------------\n\n");
  return 0;
}

int main(int argc, char *argv[]){
  double value;
  value = RunSimulation(argv[1]);

  return value;
}