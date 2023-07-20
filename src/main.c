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

void RunSimulation(char* InputFileName){
    FILE *PropertiesDataFile, *ConfigurationsFile, *DensityProfileFile, *OrientationProfileFile, *ProfileFile;
	double PressureExcess, PressureTotal, PressureLongRangeCorrection, PressureIdealGas, CalculatedPressurePengR;
	double PotentialNonbonded, PotentialBonded, PotentialWalls, PotentialLongRangeCorrection, PotentialTotal, PotentialPerturbed;
	double WeightGhostMolecule, ChemicalPotentialLongRangeCorrection, PhiCalculated, FugacityCalculated;
	double DensityNumeric, DensityMolar, DensityMass;
	double TemperatureDecrement;
    double DisplacementAcceptanceRatio;
	double *WeightIdealChainSamples, STDWeightIdealChain;	
	int CountMovements = 0, CountBlocks = 0, CountStatesInBlock = 0;
	int ProductionFlag = 0;
	int DeltaTime = 0;
	int statesPerBlock;
	int i, j, k;
	int OptimizeAcceptanceEverySteps, OptimizeAcceptanceEveryCycles;
	int StepsShiftedTemperature;
	CONFIGURATION OldConfiguration, NewConfiguration;
	PROFILE DensityProfile, OrientationProfile;

	StartTime = time(&StartTime);
	CountProductionStates=0;
	ProbabilityCavityArray = (double*) calloc(MAX_NUMBER_MOLECULES, sizeof(double));
	ProbabilityCavitySummed = (double*) calloc(MAX_NUMBER_MOLECULES, sizeof(double));
	NumberCavityScan = (int*) calloc(MAX_NUMBER_MOLECULES, sizeof(int));
	
	MaxTranslationDistance = 0.5;
	MaxLnVolumeChange = 0.1;
	DisplacementStepsAttempted=InsertionStepsAttempted=DeletionStepsAttempted=VolumeStepsAttempted=0;
	DisplacementStepsAccepted=InsertionStepsAccepted=DeletionStepsAccepted=VolumeStepsAccepted=0;
	VolumeAcceptanceCalculated=DisplacementAcceptanceCalculated=0;

	ReferencePotential = LENNARD_JONES;
	PerturbationPotential = LENNARD_JONES;

	ReadInputFile(InputFileName);
	SimulationBox.volume = SimulationBox.xSize*SimulationBox.ySize*SimulationBox.zSize;
	Ensemble = EquilibraionEnsemble;

	OptimizeAcceptanceEveryCycles = 1;
	OptimizeAcceptanceEverySteps = OptimizeAcceptanceEveryCycles*InitialNumberMolecules;

	NumberTotalSteps = NumberTotalCycles*InitialNumberMolecules;
	NumberTotalSteps = 10*(NumberTotalSteps/10);

	NumberEquilibrationSteps = NumberEquilibrationCycles*InitialNumberMolecules;
	NumberEquilibrationSteps = 10*(NumberEquilibrationSteps/10);

	StepsCalculateProperties = CyclesCalculateProperties*InitialNumberMolecules;
	StepsCalculateProperties = 10*(StepsCalculateProperties/10);

	StepsRecordConfiguration = CyclesRecordConfiguration*InitialNumberMolecules;
	StepsRecordConfiguration = 10*(StepsRecordConfiguration/10);

	StepsGetProfiles = CyclesGetProfiles*InitialNumberMolecules;
	StepsGetProfiles = 10*(StepsGetProfiles/10);
	printf("Reference = %d, perturbation = %d\n", ReferencePotential, PerturbationPotential);
	printf("-----------------------------------------------------Inputs-------------------------------------------\n");
	printf("Equilibration ensemble:....................: %s\n", GetEnsembleLabel(EquilibraionEnsemble));
	printf("Simulation ensemble:.......................: %s\n", GetEnsembleLabel(SimulationEnsemble));
	printf("Temperature:...............................: %f K\n", SimulationTemperature);
	printf("Number total cycles:.......................: %d\n", NumberTotalCycles);
	printf("Number equilibration cycles:...............: %d\n", NumberEquilibrationCycles);
	printf("Number total steps:........................: %d\n", NumberTotalSteps);
	printf("Number equilibration steps:................: %d\n", NumberEquilibrationSteps);
	printf("Acceptance ratio:..........................: %f\n", DisplacementAcceptanceTarget);
	printf("Cycles to calculate properties:............: %d\n", CyclesCalculateProperties);
	printf("Cycles to record configuration:............: %d\n", CyclesRecordConfiguration);
	printf("Cycles to get profiles:....................: %d\n", CyclesGetProfiles);
	printf("Steps to calculate properties:.............: %d\n", StepsCalculateProperties);
	printf("Steps to record configuration:.............: %d\n", StepsRecordConfiguration);
	printf("Steps to get profiles:.....................: %d\n", StepsGetProfiles);
	printf("Simulation box:............................: %s\n", SimulationBox.ClosedBox ? "Closed" : "Open");
	printf("Displacement probability:..................: %f\n", DisplacementAttemptProbability);
	printf("Inserion probability:......................: %f\n", InsertionAttemptProbability);
	printf("Deletion probability:......................: %f\n", DeletionAttemptProbability);
	printf("Number molecules:..........................: %d\n", InitialNumberMolecules);
	printf("Chain size:................................: %d\n", ChainSize);
	printf("Fugacity:..................................: %f MPa\n", SimulationFugacity/1E6);
	printf("Pressure:..................................: %f MPa\n", SimulationPressure/1E6);
	printf("Reference potential........................: %s\n", GetPotentialTypeLabel(ReferencePotential));
	printf("Perturbation potential.....................: %s\n", GetPotentialTypeLabel(PerturbationPotential));
	printf("X size.....................................: %f\n", SimulationBox.xSize);
	printf("Y size.....................................: %f\n", SimulationBox.ySize);
	printf("Z size.....................................: %f\n", SimulationBox.zSize);
	printf("------------------------------------------------------------------------------------------------------\n\n");
	
	Temperature = SimulationTemperature;
	Beta = 1/(BOLTZMANN_CONSTANT*Temperature);

	InitializeConfiguration(&OldConfiguration);
	InitializeConfiguration(&NewConfiguration);

	printf("--------------------------------------------Ideal chain calculations----------------------------------\n");
	WeightIdealChainSamples = (double*) calloc(NUMBER_SAMPLES_IDEAL_CHAIN, sizeof(double));
	WeightIdealChain = 0.0;
	for(int i=0; i<NUMBER_SAMPLES_IDEAL_CHAIN; i++){
		WeightIdealChainSamples[i] = GetRosenbluthWeightIdealChain(OldConfiguration.Molecules[0]);
		WeightIdealChain += WeightIdealChainSamples[i];
	}
	WeightIdealChain /= NUMBER_SAMPLES_IDEAL_CHAIN;
	STDWeightIdealChain = 0.0;
	for(int i=0; i<NUMBER_SAMPLES_IDEAL_CHAIN; i++){
		STDWeightIdealChain += Squared(WeightIdealChainSamples[i] - WeightIdealChain);
	}
	STDWeightIdealChain /= (NUMBER_SAMPLES_IDEAL_CHAIN - 1);


	
	printf("Weight ideal chain:........................: %f\n", WeightIdealChain);
	printf("------------------------------------------------------------------------------------------------------\n\n");

	printf("\n");
	printf("----------------------------------------Generating initial configuration------------------------------\n");
	printf("\n");
	GenerateInitialConfiguration(&OldConfiguration);
	
	printf("Initial configuration generated\n");

	// Initialize output files
	PropertiesDataFile = InitializePropertiesDataFile(InputFileName);
	ConfigurationsFile = InitializeTrajectoryFile(InputFileName);
	
	if(SimulationBox.ClosedBox && ChainSize > 1)
		ProfileFile = InitializeProfilesFile(InputFileName);
	
	RecordConfiguration(ConfigurationsFile, OldConfiguration);

	PotentialNonbonded = GetPotentialNonbonded(OldConfiguration, ReferencePotential);
	PotentialBonded = GetPotentialBonded(OldConfiguration);
	PotentialLongRangeCorrection = GetPotentialLongRangeCorrection(OldConfiguration);

	if(!SimulationBox.ClosedBox){
		PotentialWalls = 0;		
		PressureIdealGas = GetPressureIdealGas(OldConfiguration.NumberMolecules, SimulationBox.volume);
		PressureExcess = GetPressureExcess(OldConfiguration);
		PressureLongRangeCorrection = GetPressureLongRangeCorrection(OldConfiguration);
		PressureTotal = PressureIdealGas + PressureLongRangeCorrection + PressureExcess;
	}else{
		PotentialWalls = GetTotalPotentialExternal(OldConfiguration);

		PressureIdealGas = 0;
		PressureExcess = 0;
		PressureLongRangeCorrection = 0;		

		if(ChainSize>1){
			DensityProfile = ConfigureProfile(NumberBinsForProfiles);
			OrientationProfile = ConfigureProfile(NumberBinsForProfiles);
			DetermineDensityProfile(OldConfiguration, &DensityProfile);
			DetermineOrientationProfile(OldConfiguration, &OrientationProfile);
		}
	}

	DensityMass = GetDensityMass(OldConfiguration);
	DensityMolar = GetDensityMolar(OldConfiguration);

	PressureTotal = PressureIdealGas + PressureLongRangeCorrection + PressureExcess;
	PotentialTotal = PotentialNonbonded + PotentialBonded + PotentialLongRangeCorrection + PotentialWalls;

	AveragePressure=AveragePressureExcess=AveragePressureIdealGas=AveragePressureLongRangeCorrection=0;
	AveragePotential=AveragePotentialBonded=AveragePotentialNonbonded=AveragePotentialLongRangeCorrection=AveragePotentialWalls=0;
	AverageNumberMolecules=AverageDensityMass=AverageDensityMolar=0;

	printf("\n");
	printf("-------------------------------------------Simulation starting----------------------------------------\n");
	printf("\n");
	
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
			if((i%StepsGetProfiles==0 && SimulationBox.ClosedBox && ChainSize > 1)){
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

			if(Ensemble!=muVT){
				WeightGhostMolecule = GetRosenbluthWeightGhostMolecule(OldConfiguration);
			}else{
				WeightGhostMolecule = 0.0;
			}

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
			PotentialLongRangeCorrection = GetPotentialLongRangeCorrection(OldConfiguration);
			PotentialWalls = GetTotalPotentialExternal(OldConfiguration);
			
			PotentialPerturbed = GetPotentialNonbonded(OldConfiguration, PerturbationPotential);
			PotentialTotal = PotentialNonbonded + PotentialBonded + PotentialLongRangeCorrection + PotentialWalls;

			fprintf(PropertiesDataFile, "%d ", i);
			fprintf(PropertiesDataFile, "%s ", GetEnsembleLabel(Ensemble));
			fprintf(PropertiesDataFile, "%d ", ProductionFlag);
			fprintf(PropertiesDataFile, "%f ", DisplacementAcceptanceCalculated);
			fprintf(PropertiesDataFile, "%f ", MaxTranslationDistance);
			fprintf(PropertiesDataFile, "%f ", Temperature);
			fprintf(PropertiesDataFile, "%f ", SimulationBox.xSize);
			fprintf(PropertiesDataFile, "%f ", SimulationBox.ySize);
			fprintf(PropertiesDataFile, "%f ", SimulationBox.zSize);
			fprintf(PropertiesDataFile, "%d ", OldConfiguration.NumberMolecules);
			fprintf(PropertiesDataFile, "%f ", DensityMass);
			fprintf(PropertiesDataFile, "%f ", DensityMolar);
			fprintf(PropertiesDataFile, "%e ", PotentialTotal);
			fprintf(PropertiesDataFile, "%e ", PotentialBonded);
			fprintf(PropertiesDataFile, "%e ", PotentialNonbonded);
			fprintf(PropertiesDataFile, "%e ", PotentialLongRangeCorrection);
			fprintf(PropertiesDataFile, "%e ", PotentialWalls);
			fprintf(PropertiesDataFile, "%e ", PotentialPerturbed);
			fprintf(PropertiesDataFile, "%e ", WeightIdealChain);
			fprintf(PropertiesDataFile, "%e ", STDWeightIdealChain);
			fprintf(PropertiesDataFile, "%e ", WeightGhostMolecule);
			fprintf(PropertiesDataFile, "%e ", PressureTotal);
			fprintf(PropertiesDataFile, "%e ", PressureExcess);
			fprintf(PropertiesDataFile, "%e ", PressureIdealGas);
			fprintf(PropertiesDataFile, "%e ", PressureLongRangeCorrection);
		
			fprintf(PropertiesDataFile, "\n");
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
          DisplacementAcceptanceRatio = DisplacementAcceptanceCalculated/DisplacementAcceptanceTarget;
      		DisplacementAcceptanceRatio = LimMinMax(DisplacementAcceptanceRatio, 0.3, 1.7);
			MaxTranslationDistance *= DisplacementAcceptanceRatio;
			MaxTranslationDistance = LimMinMax(MaxTranslationDistance, 0.05, 10.0);
		}

        // adjust maximum volume change according to the acceptance ratio
		VolumeAcceptanceCalculated = (double) VolumeStepsAccepted/VolumeStepsAttempted;
		if(Ensemble==NPT && VolumeStepsAttempted%10==0 && i<NumberEquilibrationSteps/2){
			if(VolumeAcceptanceCalculated - VolumeAcceptanceTarget > 0.05){
				MaxLnVolumeChange *= 1.05;
			}else if(VolumeAcceptanceCalculated - VolumeAcceptanceTarget < -0.05){
				MaxLnVolumeChange *= 0.95;
			}
			MaxLnVolumeChange = LimMinMax(MaxLnVolumeChange, 0.01, 0.7);
		}
	}
	
	RecordSimulationLog(InputFileName);
	RecordConfiguration(ConfigurationsFile, OldConfiguration);
	if(SimulationBox.ClosedBox) fclose(ProfileFile);
	fclose(ConfigurationsFile);
	if(SimulationBox.ClosedBox && ChainSize > 1)
		fclose(PropertiesDataFile);

	printf("\n");
	printf("-----------------------------------------Simulation has finished--------------------------------------\n");
	printf("\n");
}

int main(int argc, char *argv[]){
	
	RunSimulation(argv[1]);

	return 0;
}