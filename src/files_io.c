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

#include"files_io.h"

// ---------------------------------------------------------------------------------------------------------------
// File parsing functions
// ---------------------------------------------------------------------------------------------------------------

/* ***************************************************************************************************************
 * Name       | GetEnsembleLabel                                                           
 * --------------------------------------------------------------------------------------------------------------- 
 * Function   | Returns the label corresponding to the given ensemble.                             
 * Parameters | InputEnsemble: An enumeration value representing an ensemble (NVT, NPT, or muVT).                           
 * **************************************************************************************************************/
char* GetEnsembleLabel(enum ensemble InputEnsemble){
  if(InputEnsemble==NVT){
   return "NVT";
  }else if(InputEnsemble==NPT){
   return "NPT";
  }else if(InputEnsemble==muVT){
   return "muVT";
  }
}


/* ***************************************************************************************************************
 * Name       | GetPotentialTypeLabel                                                           
 * --------------------------------------------------------------------------------------------------------------- 
 * Function   | Returns the label corresponding to the given potential type.                             
 * Parameters | Potential: An enumeration value representing a potential type (LENNARD_JONES or HARD_SPHERE).                           
 * **************************************************************************************************************/
char* GetPotentialTypeLabel(enum PotentialType Potential){
  if(Potential==LENNARD_JONES){
   return "Lennard-Jones";
  }else if(Potential==HARD_SPHERE){
   return "Hard-sphere";
  }
}

/* ***************************************************************************************************************
 * Name       | RemoveSubstring                                                           
 * --------------------------------------------------------------------------------------------------------------- 
 * Function   | Removes the given substring from the given string.                             
 * Parameters | string: A pointer to the string from which the substring needs to be removed.                           
 *            | sub: A pointer to the substring that needs to be removed from the string.  
 * **************************************************************************************************************/
void RemoveSubstring (char *string, char *sub) {
  char *match;
  int len = strlen(sub);
  while ((match = strstr(string, sub))) {
   *match = '\0';
   strcat(string, match+len);
  }
}

/* ***************************************************************************************************************
 * Name       | ReadInitialConfiguration                                                            
 * --------------------------------------------------------------------------------------------------------------- 
 * Function   | Reads the initial configuration from a file and stores it in a CONFIGURATION structure.      
 * Parameters | InitialConfigurationFilePath: The path of the file containing the initial configuration.               
 *            | Configuration: A pointer to a CONFIGURATION structure where the data will be stored.    
 * **************************************************************************************************************/   
void ReadInitialConfiguration(char InitialConfigurationFilePath[], CONFIGURATION* Configuration){
  FILE* InitialConfiguraitonFile = fopen(InitialConfigurationFilePath, "r");
  CONFIGURATION InitialConfiguraiton;
  VECTOR AtomPosition;
  char  str[100], aux1[100], aux2[100], aux3[100];
  double x, y, z;
  int  NumberOfParticles, CCounter = 0, MoleculeIndex, ParticleIndex;

  fscanf(InitialConfiguraitonFile, "%d", &NumberOfParticles);
  fscanf(InitialConfiguraitonFile, "%s %s %s", aux1, aux2, aux3);
  for(int i = 0; i < NumberOfParticles; i++){
   fscanf(InitialConfiguraitonFile, "%s %lf %lf %lf", aux1, &AtomPosition.x, &AtomPosition.y, &AtomPosition.z);
   if(strcmp(aux1, "C")==0){
    MoleculeIndex = CCounter/ChainSize;
    ParticleIndex = CCounter % ChainSize;
    Configuration->Molecules[MoleculeIndex].Atoms[ParticleIndex].Position = AtomPosition;
    CCounter++;
   }
  }

  if(CCounter == InitialNumberMolecules*ChainSize){
   printf("Initial configuration file loaded successfully. Don't panic!\n");
   InitialConfiguraiton.NumberMolecules = InitialNumberMolecules;
 }else{
  printf("Something is wrong in the initial configuration file!\n");
  printf("Number of carbon atoms read %d, expected %d\n", CCounter, InitialNumberMolecules*ChainSize);
 }
 fclose(InitialConfiguraitonFile);
}


/* ***************************************************************************************************************
 * Name       | ReadInputFile                                                            
 * --------------------------------------------------------------------------------------------------------------- 
 * Function   | Reads the input file and sets the corresponding variables.      
 * Parameters | inputsFilePath: The path to the input file.               
 * **************************************************************************************************************/
void ReadInputFile(char inputsFilePath[]){
  FILE* InputFile = fopen(inputsFilePath, "r");
  char  InputName[100], InputData[100];
  int  int_aux;
  double double_aux;

  while(!feof(InputFile)){
   fscanf(InputFile, "%s %s", InputName, InputData);
   
   if(strcmp(InputName, "NUBMER_OF_MOLECULES")==0){
    sscanf(InputData, "%d", &int_aux);
    InitialNumberMolecules = int_aux;

   }else if(strcmp(InputName, "CHAIN_SIZE")==0){
    sscanf(InputData, "%d", &int_aux);
    ChainSize = int_aux;

   }else if(strcmp(InputName, "TEMPERATURE")==0){
    sscanf(InputData, "%lf", &double_aux);
    SimulationTemperature = double_aux;
     
   }else if(strcmp(InputName, "BOX_X_SIZE")==0){
    sscanf(InputData, "%lf", &double_aux);
    SimulationBox.xSize = double_aux;
    SimulationBox.xMin = -double_aux/2;
    SimulationBox.xMax = double_aux/2;
     
   }else if(strcmp(InputName, "BOX_Y_SIZE")==0){
    sscanf(InputData, "%lf", &double_aux);
    SimulationBox.ySize = double_aux;
    SimulationBox.yMin = -double_aux/2;
    SimulationBox.yMax = double_aux/2;

   }else if(strcmp(InputName, "BOX_Z_SIZE")==0){
    sscanf(InputData, "%lf", &double_aux);
    SimulationBox.zSize = double_aux;
    SimulationBox.zMin = -double_aux/2;
    SimulationBox.zMax = double_aux/2;

   }else if(strcmp(InputName, "NUMBER_TOTAL_STEPS")==0){
    sscanf(InputData, "%lf", &double_aux);
    int_aux = (int) double_aux;
    NumberTotalSteps = int_aux;
     
   }else if(strcmp(InputName, "NUMBER_EQUILIBRATION_STEPS")==0){
    sscanf(InputData, "%lf", &double_aux);
    int_aux = (int) double_aux;
    NumberEquilibrationSteps = int_aux;

   }else if(strcmp(InputName, "ACCEPTANCE_RATIO")==0){
    sscanf(InputData, "%lf", &double_aux);
    DisplacementAcceptanceTarget = double_aux;
    VolumeAcceptanceTarget = double_aux;

   }else if(strcmp(InputName, "CLOSED_BOX")==0){
    sscanf(InputData, "%d", &int_aux);
    SimulationBox.ClosedBox = int_aux == 1 ? true : false;

   }else if(strcmp(InputName, "NUMBER_OF_BINS_FOR_PROFILE")==0){
    sscanf(InputData, "%d", &int_aux);
    NumberBinsForProfiles = int_aux;

   }else if(strcmp(InputName, "DISPLACEMENT_PROBABILITY")==0){
    sscanf(InputData, "%lf", &double_aux);
    DisplacementAttemptProbability = double_aux;

   }else if(strcmp(InputName, "INSERTION_PROBABILITY")==0){
    sscanf(InputData, "%lf", &double_aux);
    InsertionAttemptProbability = double_aux;

   }else if(strcmp(InputName, "DELETION_PROBABILITY")==0){
    sscanf(InputData, "%lf", &double_aux);
    DeletionAttemptProbability = double_aux;

   }else if(strcmp(InputName, "FUGACITY")==0){
    sscanf(InputData, "%lf", &double_aux);
    SimulationFugacity = double_aux;
     
   }else if(strcmp(InputName, "SIMULATION_ENSEMBLE")==0){
    if(strcmp(InputData, "NVT")==0){
     SimulationEnsemble = NVT;
    }else if(strcmp(InputData, "NPT")==0){
     SimulationEnsemble = NPT;
    }else if(strcmp(InputData, "muVT")==0){
     SimulationEnsemble = muVT;
    }

   }else if(strcmp(InputName, "EQUILIBRATION_ENSEMBLE")==0){
    if(strcmp(InputData, "NVT")==0){
     EquilibraionEnsemble = NVT;
    }else if(strcmp(InputData, "NPT")==0){
     EquilibraionEnsemble = NPT;
    }else if(strcmp(InputData, "muVT")==0){
     EquilibraionEnsemble = muVT;
    }

   }else if(strcmp(InputName, "REFERENCE_POTENTIAL")==0){
    if(strcmp(InputData, "LENNARD_JONES")==0){
     ReferencePotential = LENNARD_JONES;
    }else if(strcmp(InputData, "HARD_SPHERE")==0){
     ReferencePotential = HARD_SPHERE;
    }

   }else if(strcmp(InputName, "PERTURBED_POTENTIAL")==0){
    if(strcmp(InputData, "LENNARD_JONES")==0){
     PerturbationPotential = LENNARD_JONES;
    }else if(strcmp(InputData, "HARD_SPHERE")==0){
     PerturbationPotential = HARD_SPHERE;
    }

   }else if(strcmp(InputName, "NUMBER_EQUILIBRATION_CYCLES")==0){
    sscanf(InputData, "%lf", &double_aux);
    int_aux = (int) double_aux;
    NumberEquilibrationCycles = int_aux;
     
   }else if(strcmp(InputName, "NUMBER_TOTAL_CYCLES")==0){
    sscanf(InputData, "%lf", &double_aux);
    int_aux = (int) double_aux;
    NumberTotalCycles = int_aux;
     
   }else if(strcmp(InputName, "CYCLES_TO_CALCULATE_PROPERTIES")==0){
    sscanf(InputData, "%lf", &double_aux);
    int_aux = (int) double_aux;
    CyclesCalculateProperties = int_aux;

   }else if(strcmp(InputName, "CYCLES_TO_RECORD_CONFIGURAITON")==0){
    sscanf(InputData, "%lf", &double_aux);
    int_aux = (int) double_aux;
    CyclesRecordConfiguration = int_aux;

   }else if(strcmp(InputName, "CYCLES_TO_GET_PROFILES")==0){
    sscanf(InputData, "%lf", &double_aux);
    int_aux = (int) double_aux;
    CyclesGetProfiles = int_aux;

   }
  }
  SimulationBox.volume = SimulationBox.xSize*SimulationBox.ySize*SimulationBox.zSize;
}

/* ***************************************************************************************************************
 * Name       | RecordSimulationLog                                                            
 * --------------------------------------------------------------------------------------------------------------- 
 * Function   | Records the final results of the full simulation
 * Parameters | BaseName name for the output name (based on the .inp input) with .out.log
 * **************************************************************************************************************/
void RecordSimulationLog(char BaseName[]){
  char           LogFileName[100];
  char           EnsembleLabel[10];
  FILE*          LogFile;
  int i, j;

  // Configures log file name
  strcpy(LogFileName, BaseName);
  RemoveSubstring(LogFileName, ".inp");
  strcat(LogFileName, ".out.log");

  // Open file
  LogFile = fopen(LogFileName, "w");

  // Print basic header
  fprintf(
    LogFile, 
    "================================================================================================================\n"
    "This program performs a configurational-bias monte carlo simulation for pure linear alkanes. The simulations can\n"
    "be performed either in the canonical or grand-canonical ensembles, it is also possible to run the simuation with\n"
    "open or closed boxes. For the latter case, the external walls are modelled using a Steele-10-5-4 potential which\n"
    "represent the adsorbent walls. For such cases, it is calculated the density and orienation profiles following \n"
    "second Legendre polynomial.\n\n"
    "Author: Rodolfo J Amancio (rodolfojamancio@gmal.com)\n"
    "University of Campinas, School of chemical engineering, Campinas/SP - Brazil\n"
    "================================================================================================================\n\n"
  );

  // Print basic input data for the simulation
  fprintf(
    LogFile, 
    "================================================================================================================\n"
    "                     Simulation setup:\n"
    "----------------------------------------------------------------------------------------------------------------\n"
    "Chain size:.....................................: %d\n"
    "Temperature:....................................: %.2f K\n"
    "Number of total steps:..........................: %d\n"
    "Number of equilibration steps:..................: %d\n"
    "Acceptance ratio target:........................: %f\n"
    "Steps to calculate properties...................: %d\n"
    "Steps to record configuration...................: %d\n"
    "Simulation ensemble:............................: %s\n", 
    ChainSize, Temperature, NumberTotalSteps,
    NumberEquilibrationSteps, DisplacementAcceptanceTarget, 
    StepsCalculateProperties, StepsRecordConfiguration,
    GetEnsembleLabel(SimulationEnsemble)
  );

  if (SimulationEnsemble == NVT) {
      fprintf(
        LogFile, 
        "Number of molecules:............................: %d\n"
        "x size:.........................................: %.2f A\n"
        "y size:.........................................: %.2f A\n"
        "z size:.........................................: %.2f A\n"
        "Box volume:.....................................: %.2f A^3\n",
        InitialNumberMolecules, SimulationBox.xSize, SimulationBox.ySize, 
        SimulationBox.zSize, SimulationBox.volume
      );
  } else if (SimulationEnsemble == NPT) {
      fprintf(
        LogFile, 
        "Number of molecules:............................: %d\n"
        "Pressure:.......................................: %f MPa\n",
        InitialNumberMolecules, SimulationPressure / 1E6
      );
  } else if (SimulationEnsemble == muVT) {
      fprintf(
        LogFile, 
        "Fugacity:.......................................: %f MPa\n"
        "x size:.........................................: %.2f A\n"
        "y size:.........................................: %.2f A\n"
        "z size:.........................................: %.2f A\n"
        "Box volume:.....................................: %.2f A^3\n",
        SimulationFugacity / 1E6, SimulationBox.xSize, SimulationBox.ySize, 
        SimulationBox.zSize, SimulationBox.volume
      );
  }

  fprintf(LogFile, "----------------------------------------------------------------------------------------------------------------\n");

  if (SimulationBox.ClosedBox) {
      fprintf(
        LogFile, 
        "----------------------------------------------------------------------------------------------------------------\n"
        "                    Simulation in a closed box\n"
        "Number of bins for profiles.....................: %d\n",
        NumberBinsForProfiles
      );
  }

  fprintf(
    LogFile, 
    "================================================================================================================\n\n"
    "================================================================================================================\n"
    "                     Simulation summary:\n"
    "Acceptance ratio calculated:....................: %f\n"
    "Deletions attempted:............................: %d\n"
    "Deletions accepted:.............................: %d\n"
    "Insertions attempted:...........................: %d\n"
    "Insertions accepted:............................: %d\n"
    "Displacements attempted:........................: %d\n"
    "Displacements accepted:.........................: %d\n"
    "Volume changes attempted:.......................: %d\n"
    "Volume changes accepted:........................: %d\n"
    "----------------------------------------------------------------------------------------------------------------\n"
    "Potential:......................................: %e J\n"
    "Potential bonded:...............................: %e J\n"
    "Potential nonbonded:............................: %e J\n"
    "Potential walls:................................: %e J\n"
    "Potential long range correction:................: %e J\n"
    "----------------------------------------------------------------------------------------------------------------\n"
    "----------------------------------------------------------------------------------------------------------------\n"
    "Average density mass:...........................: %f kg/m³\n"
    "Average density molar:..........................: %f mol/m³\n"
    "----------------------------------------------------------------------------------------------------------------\n",
    DisplacementAcceptanceCalculated, DeletionStepsAttempted, DeletionStepsAccepted, InsertionStepsAttempted,
    InsertionStepsAccepted, DisplacementStepsAttempted, DisplacementStepsAccepted, VolumeStepsAttempted,
    VolumeStepsAccepted, (AveragePotential / CountProductionStates), (AveragePotentialBonded / CountProductionStates),
    (AveragePotentialNonbonded / CountProductionStates), (AveragePotentialWalls / CountProductionStates),
    (AveragePotentialLongRangeCorrection / CountProductionStates), (AverageDensityMass / CountProductionStates),
    (AverageDensityMolar / CountProductionStates)
  );


  if (!SimulationBox.ClosedBox) {
    if (Ensemble != NPT) {
        float pressure = 1E-6 * (AveragePressure / CountProductionStates);
        fprintf(LogFile, "----------------------------------------------------------------------------------------------------------------\n");
        fprintf(LogFile, "Pressure:.......................................: %f MPa\n", pressure);
        fprintf(LogFile, "Pressure ideal gas:.............................: %f MPa\n", 1E-6 * (AveragePressureIdealGas / CountProductionStates));
        fprintf(LogFile, "Pressure excess:................................: %f MPa\n", 1E-6 * (AveragePressureExcess / CountProductionStates));
        fprintf(LogFile, "Pressure long range correction:.................: %f MPa\n", 1E-6 * (AveragePressureLongRangeCorrection / CountProductionStates));
        fprintf(LogFile, "----------------------------------------------------------------------------------------------------------------\n");
    }
    if (Ensemble != muVT) {
        float weightGhostMolecule = AverageWeightGhostMolecule / CountProductionStates;
        float fugacity = 1E-6 * (AveragePressureIdealGas / weightGhostMolecule);
        float fugacityCoefficient = fugacity / (1E-6 * (AveragePressure / CountProductionStates));

        fprintf(LogFile, "----------------------------------------------------------------------------------------------------------------\n");
        fprintf(LogFile, "Weight ghost molecule:..........................: %f \n", weightGhostMolecule);
        fprintf(LogFile, "Fugacity:.......................................: %f Mpa\n", fugacity);
        fprintf(LogFile, "Fugacity coefficient:...........................: %f \n", fugacityCoefficient);
        fprintf(LogFile, "----------------------------------------------------------------------------------------------------------------\n");
    }
  }

  time(&EndTime);
  int ElapsedTimeSeconds = EndTime - StartTime;
  int Days = ElapsedTimeSeconds / (3600 * 24);
  int DaysRest = ElapsedTimeSeconds % (3600 * 24);
  int Hours = DaysRest / 3600;
  int HoursRest = DaysRest % 3600;
  int Minutes = HoursRest / 60;
  int MinutsRest = HoursRest % 60;
  int Seconds = MinutsRest % 60;

  fprintf(LogFile, "\n----------------------------------------------------------------------------------------------------------------\n");
  fprintf(LogFile, "Simulations has finished. Elapsed time: %d days, %d hours, %d minutes and %d seconds.\n", Days, Hours, Minutes, Seconds);
  fprintf(LogFile, "----------------------------------------------------------------------------------------------------------------\n");

  fclose(LogFile);
}

/* ***************************************************************************************************************
 * Name       | InitializePropertiesDataFile                                                            
 * --------------------------------------------------------------------------------------------------------------- 
 * Function   | Creates the header and return the properties file
 * Parameters | BaseName name for the output name (based on the .inp input) _properties_data.out.dat
 * Returns    | FILE*: Pointer to the created properties file
 * **************************************************************************************************************/
FILE* InitializePropertiesDataFile(char BaseName[]) {
    FILE* PropertiesDataFile;

    // Generates base name
    char PropertyDataFileName[100];
    strcpy(PropertyDataFileName, BaseName);
    RemoveSubstring(PropertyDataFileName, ".inp");
    strcat(PropertyDataFileName, "_properties_data.out.dat");

    PropertiesDataFile = fopen(PropertyDataFileName, "w");

    // Header for the pressure-potential data file
    const char* headers[] = {
      "i", "ensemble", "production", "displacement_acceptance", "max_displacement", 
      "temperature", "x_size", "y_size", "z_size", "number_molecules", "density_mass", 
      "density_mol", "potential_total", "potential_bonded", "potential_nonbonded", 
      "potential_lrc", "potential_walls", "potential_perturbed", "weight_ideal_chain", 
      "std_weight_ideal_chain", "weight_ghost_molecule", "pressure_total", 
      "pressure_excess", "pressure_ideal", "pressure_lrc"
    };

    int numHeaders = sizeof(headers) / sizeof(headers[0]);

    for (int i = 0; i < numHeaders; i++) {
        fprintf(PropertiesDataFile, "%s ", headers[i]);
    }

    fprintf(PropertiesDataFile, "\n");

    return PropertiesDataFile;
}

/* ***************************************************************************************************************
 * Name       | InitializeProfilesFile                                                            
 * --------------------------------------------------------------------------------------------------------------- 
 * Function   | Creates the header and returns the properties file.
 * Parameters | BaseName: Name for the output file (based on the .inp input)
 * Returns    | FILE*: Pointer to the created properties file
 * **************************************************************************************************************/
FILE* InitializeProfilesFile(char BaseName[]){
  FILE* ProfileFile;
  char ProfileFileName[100];
  strcpy(ProfileFileName, BaseName);
  RemoveSubstring(ProfileFileName, ".inp");
  strcat(ProfileFileName, "_raw_profile.out.dat");
  ProfileFile = fopen(ProfileFileName, "w");
  return ProfileFile;
}

/* ***************************************************************************************************************
 * Name       | RecordProfileData                                                            
 * --------------------------------------------------------------------------------------------------------------- 
 * Function   | Stores profile data to ProfileFile
 * Parameters | 
 *              - FILE* ProfileFile: file pointer to record data
 *              - CONFIGURATION Configuration current configuration with profile data
 * **************************************************************************************************************/
void RecordProfileData(FILE* ProfileFile, CONFIGURATION Configuration){
  fprintf(
   ProfileFile, 
   "%f %f %f %d %f\n", 
   SimulationBox.xSize, 
   SimulationBox.ySize, 
   SimulationBox.zSize, 
   NumberBinsForProfiles, 
   Configuration.Molecules[0].MolarMass
  );
  for(int i = 0; i < Configuration.NumberMolecules; i++){
   fprintf(ProfileFile, "%7.2f %6.3f\n", Configuration.Molecules[i].CenterOfMass.z, Configuration.Molecules[i].OrderParameter);
  }
  fprintf(ProfileFile, "\n");
}


/***************************************************************
 * Name       | RecordInitialConfiguraiton
 * -------------------------------------------------------------
 * Function   | Records the initial configuration to a file.
 * Parameters |
 *              - BaseName: Name for the output file (based on the .inp input)
 *              - Configuration: Structure holding the configuration details
 *              - InitialConfigurationFile: Pointer to the output file
 ****************************************************************/
void RecordInitialConfiguraiton(char BaseName[], CONFIGURATION Configuration, FILE* InitialConfigurationFile) {
  // Generates base name
  char ConfigurationsFileName[100];
  strcpy(ConfigurationsFileName, BaseName);
  RemoveSubstring(ConfigurationsFileName, ".inp");
  strcat(ConfigurationsFileName, "_initial_configuraion.out.xyz");
  
  fprintf(InitialConfigurationFile, "%d\n\n", Configuration.NumberMolecules * Configuration.Molecules[0].Size);
  
  char buffer[256];  // Buffer to store formatted output
  
  for (int i = 0; i < Configuration.NumberMolecules; i++) {
    for (int j = 0; j < ChainSize; j++) {
      sprintf(buffer, "C  %.4f  %.4f  %.4f\n", Configuration.Molecules[i].Atoms[j].Position.x,
              Configuration.Molecules[i].Atoms[j].Position.y, Configuration.Molecules[i].Atoms[j].Position.z);
      fputs(buffer, InitialConfigurationFile);
    }
  }
}


/***************************************************************
 * Name       | InitializeTrajectoryFile
 * -------------------------------------------------------------
 * Function   | Initializes the trajectory file for writing.
 * Parameters |
 *              - BaseName: Name for the output file (based on the .inp input)
 * Returns    | Pointer to the trajectory file
 ****************************************************************/
FILE* InitializeTrajectoryFile(char BaseName[]) {
  // Generates base name
  char ConfigurationsFileName[100];
  strcpy(ConfigurationsFileName, BaseName);
  RemoveSubstring(ConfigurationsFileName, ".inp");
  strcat(ConfigurationsFileName, "_configurations.out.xyz");

  FILE* ConfigurationsFile = fopen(ConfigurationsFileName, "w");
  if (ConfigurationsFile == NULL) {
    perror("Error opening file");
    return NULL;  // Return NULL on error
  } 
  return ConfigurationsFile;
}

/***************************************************************
 * Name       | RecordConfiguration
 * -------------------------------------------------------------
 * Function   | Records the given configuration to the specified file.
 * Parameters |
 *              - IntermediateConfigurationFile: File pointer for writing
 *              - Configuration: The configuration to be recorded
 * Returns    | None
 ****************************************************************/
void RecordConfiguration(FILE* IntermediateConfigurationFile, CONFIGURATION Configuration) {
  int numAtoms = Configuration.NumberMolecules * ChainSize;

  fprintf(IntermediateConfigurationFile, "%d\n\n", numAtoms);

  for (int i = 0; i < Configuration.NumberMolecules; i++) {
    for (int j = 0; j < Configuration.Molecules[i].Size; j++) {
      char atomString[50];
      sprintf(
        atomString, 
        "C  %.4f  %.4f  %.4f\n", 
        Configuration.Molecules[i].Atoms[j].Position.x, 
        Configuration.Molecules[i].Atoms[j].Position.y, 
        Configuration.Molecules[i].Atoms[j].Position.z
      );
      fprintf(IntermediateConfigurationFile, "%s", atomString);
    }
  }
}
