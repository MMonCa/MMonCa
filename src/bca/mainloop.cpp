 /*
 * Original author: jesus.hernandez.mangas@tel.uva.es
 *
 * Copyright 2014 University of Valladolid, Valladolid, Spain.
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 * 
 *    http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 *
 * Adapted to MMonCa: I. Martin-Bragado. IMDEA Materials Institute.
 *
 */ 
///////////////////////////////////////////////////////////////////////////
//  
// MAINLOOP.CPP
//
///////////////////////////////////////////////////////////////////////////
#include "mainloop.h"

using std::cerr;
using std::endl;
using std::cout;
using std::ifstream;
using std::ofstream;
using std::cin;

namespace BCA {

extern Hist     HstPROJ[N_ATOMS];    // Statistics 1D
extern Histo3D  HstPROJ3D[N_ATOMS];  // Statistics 3D

// Assign the default value to the variables /////////////////////////////
void MainLoop::Default()
 {
   InteractionRadius = 2.7155;  //    (A)
   XTalSize          = 0.5;     //    (lattice units)
   ThresholdEnergy   = 10;      // To follow cascade (eV)
   GhiLimit          =  0.1;
   SimultaneousDistance = 0.5;  //  SIMULTANEOUS COLLISION (A)
   SimultaneousDistance2 = 0.25;// (A)
   POS(0,0,0);
   DivType           = 1;       // 0 : Uniform distribution
                                // 1 : Cosine distribution
   Energy_Spread     = 0;       // 0 : No energy distribution
                                // 1 : Uniform distribution
                                // 2 : Gaussian distribution
   Energy_Sigma      = 0.0;     // Sigma value for energy distribution
   Energy_Percentage = 0.0;
   NumberOfImplants  = 10;
   Verbose = 0;                 // 0 = Do not show extra information
                                // 1 = Show extra information
   Randomize         = 1;       // Say if must be random or not
   Seed1         = 1234;
   Seed2         = 5678;

   ABC_X (1,0,0);               // Eje <100> = (1,0,0) DIR=(1,0,0)
                                //     <110> = (1,1,0) DIR=(1,1,0)
   ABC_Y (0,0,0);               // Wafer flat direction 
                                //  (default value dependant on ABC_X)
   Tha_Cut = 0.0;               // Cut tilt angle
   Phi_Cut = 0.0;               // Cut rotation angle
   WinA          = 5.431;       // Window Width (A) (Axe Z) Horizontal
   WinB          = 5.431;       // Window Width (A) (Axe Y) Vertical
   IP            = 1e10;        // At first time randomize
   MaxItPerImplant   = 100000L;
   MaxDepth          = 1e10;    // Profundidad m�xima despu�s de la cual
                                //  descarto seguir el ion. Lo coloco ah�.
   MaxFlyPerImplant  = 50L;     // Modificado para que cuente el numero
                                //  max. de no colisiones SEGUIDAS en
                                //  una implantacion
   LevelOfSimulation = 1;       // 0 = All levels. No restrictions
                                // 1 = Only follows ION trajectory
   Therm = 1;                   // Select or not thermal displacements
                                //  ( =0)    No thermal displacements
                                //  (!=0) Thermal displacements allowed
   PearsonIV   = 0;             // Print Pearson IV fitted coefficients
                                //   (=0)  No print
                                //   (!=0) Print
   pE = 0.0; //0.10;            // Energy percentage for surface rare event
   pX = 0.0; //0.05;            // Depth percentage for surface rare event

   strcpy(HSTFile,"Histo1D. "); // Last space is required
   strcpy(sFrenado,"Our ");     // Last space is required
   strcpy(sIonization,"ZBL ");  // Last space is required
   Frenado = 2;     // Variable global definida en stopping.h
   Ionizacion = 1;  // Global variable defined in stopping.h

   rs0 = 1.8;   // (a.u)

   for(int i=0;i<100;i++)
    for(int j=0;j<100;j++)
     SpecificPotential[i][j]=0; // By default ZBL screening function

   PAUSE=0; // No pause at the end of calculation

   modif=0;
   // Damage accumulation  variables
   RecombinationFactor   = 0.06; // Default value for Boron
   KinchinPeaseConstant  = 0.8;
   AmorphizationDensity  = 0.1*4.99e22; // (defects/cm^3) Default value for Si
   CutOffEnergy          = 0.15; // (eV/A^3) equivale a 24.0 eV/celda en Si
   DamageDisplacementEnergy = 15.0; // (eV) Default value for silicon
   strcpy(PreviousDamageFile,""); // No previous damage
   FullDoseOutput = 1;  // (=0) Implanted dose is used in temporary outputs
                        // (=1) Full dose is used in temporary outputs
   DoseSplitting =0;     
   int i_;
   for(i_=0; i_<10 ; i_++ ) DoseDistribution[i_]=0.0;                   
   for(i_=0; i_<N_ATOMS; i_++ ) PROJS[i_]=0;
   ImplantationArea = 100.0; // 100 (A^2) = 1 nm^2
   WEffective       =  1.0; // Adimensional
   ascii=0;

   DoseRate = 1e13; /* cm^-2 s^-1 */
   
   max_r=0.0;

 // 20100309 LuisJou MOD by J.H
   Extraccion = 0;
 //
    

  }

 // Read and set the parameters ///////////////////////////////////////////

void MainLoop::SetParameters( int argc, char **argv )
 {
   SimulationDefinition SIM;
  
   Input IN( "input_file", 1 );
   if( argc > 1 ) IN = Input( argv[1], 1 ); // Read from.. and case sensitive ON
   else 
    {
     cerr<< "# Error, I do not known the name of input file.\n"; 
     exit(READ_ERROR);
    };

   IN.Add( Variable("Projectiles",    INTMATRIX, &SIM.PROJS            ) );
   IN.Add( Variable("ENERGY",              REAL, &SIM.ENERGY           ) );
   IN.Add( Variable("DIR",               VECTOR, &SIM.DIR              ) );
   IN.Add( Variable("Miller",           INTEGER, &SIM.Miller           ) );
   IN.Add( Variable("Tha",                 REAL, &SIM.Tha              ) );
   IN.Add( Variable("Phi",                 REAL, &SIM.Phi              ) );
   IN.Add( Variable("Divergence",          REAL, &SIM.Divergency       ) );
   IN.Add( Variable("NumberOfImplants", INTEGER, &SIM.NumberOfImplants ) );
   IN.Add( Variable("Dose",                REAL, &SIM.Dose             ) );
   IN.Add( Variable("Temperature",         REAL, &SIM.Temperature      ) );
   IN.Add( Variable("RS0",               MATRIX, &SIM.RS0              ) );
   IN.Add( Variable("NextSimulation",SIMULATIONDEFINITION, &SIMULACIONES, &SIM));
   
   IN.Add( Variable("EnergySpread",     INTEGER, &Energy_Spread        ) );
   IN.Add( Variable("EnergySigma",         REAL, &Energy_Sigma         ) );
   IN.Add( Variable("EnergyPercentage",    REAL, &Energy_Percentage    ) );
   IN.Add( Variable("POS",           VECTOR, &POS          ) );
   IN.Add( Variable("DivType",              INTEGER, &DivType              ) );
   IN.Add( Variable("InteractionRadius",       REAL, &InteractionRadius    ) );
   IN.Add( Variable("ThresholdEnergy",         REAL, &ThresholdEnergy      ) );
   IN.Add( Variable("GhiLimit",            REAL, &GhiLimit             ) );
   IN.Add( Variable("SimultaneousDistance",    REAL, &SimultaneousDistance ) );
   IN.Add( Variable("SD2",    REAL, &SimultaneousDistance2 ) );
   IN.Add( Variable("Verbose",             BOOLEAN, &Verbose         ) );

   IN.Add( Variable("Seed1",                INTEGER, &Seed1        ) );
   IN.Add( Variable("Seed2",                INTEGER, &Seed2        ) );
   IN.Add( Variable("Randomize",            INTEGER, &Randomize            ) );
   IN.Add( Variable("ABC",                   VECTOR, &ABC_X                ) );
   IN.Add( Variable("FLAT",                  VECTOR, &ABC_Y                ) );
   IN.Add( Variable("Cut_Tha",             REAL, &Tha_Cut          ) );
   IN.Add( Variable("Cut_Phi",             REAL, &Phi_Cut          ) );
   IN.Add( Variable("WinA",            REAL, &WinA         ) );
   IN.Add( Variable("WinB",            REAL, &WinB         ) );
   IN.Add( Variable("Therm",                BOOLEAN, &Therm        ) );
   IN.Add( Variable("MaxIterationsPerImplant", LONG, &MaxItPerImplant      ) );
   IN.Add( Variable("MaxDepth"               , REAL, &MaxDepth             ) );
   IN.Add( Variable("LevelOfSimulation",    INTEGER, &LevelOfSimulation    ) );
   IN.Add( Variable("Atom",      ATOMDEFINITION, &Atoms        ) );
   IN.Add( Variable("SpecificPotential",      ARRAY, &SpecificPotential    ) );

   IN.Add( Variable("Amorphous",        INTEGER, &LAYER.Amorphous));
   IN.Add( Variable("LatticeParameter",  VECTOR, &LAYER.LatticeParameter));
   IN.Add( Variable("Angles",            VECTOR, &LAYER.Angles));
   IN.Add( Variable("XTal",      XTALDEFINITION, &LAYER.XTal));
   IN.Add( Variable("XTal2",    XTALDEFINITION2, &LAYER.XTal));
   IN.Add( Variable("XMin",                REAL, &LAYER.XMin));
   IN.Add( Variable("XMax",                REAL, &LAYER.XMax));
   IN.Add( Variable("EDTFile[",            TEXT_, &LAYER.EDTFile));
   IN.Add( Variable("EDTCreate",        INTEGER, &LAYER.EDTCreate));
   IN.Add( Variable("NextLayer",LAYERDEFINITION, &Bulk, &LAYER));

   IN.Add( Variable("pE",                  REAL, &pE               ) );
   IN.Add( Variable("pX",                  REAL, &pX               ) );

   IN.Add( Variable("PearsonIV",        BOOLEAN, &PearsonIV        ) );
   IN.Add( Variable("FullDoseOutput",   BOOLEAN, &FullDoseOutput   ) );
   IN.Add( Variable("HSTFile[",           TEXT_, &HSTFile          ) );
   IN.Add( Variable("Stopping[",          TEXT_, &sFrenado         ) );
   IN.Add( Variable("Ionization[",        TEXT_, &sIonization      ) );

   IN.Add( Variable("PAUSE",                INTEGER, &PAUSE                ) );

   IN.Add( Variable("RecombinationFactor",     REAL, &RecombinationFactor  ) );
   IN.Add( Variable("KinchinPeaseConstant",    REAL, &KinchinPeaseConstant ) );
   IN.Add( Variable("AmorphizationDensity",    REAL, &AmorphizationDensity ) );
   IN.Add( Variable("CutOffEnergy",            REAL, &CutOffEnergy         ) );
   IN.Add( Variable("DisplacementEnergy",      REAL, &DamageDisplacementEnergy  ) );
   IN.Add( Variable("PreviousDamageFile[",    TEXT_, &PreviousDamageFile   ) );

   IN.Add( Variable("DoseSplitting",        BOOLEAN, &SIM.DoseSplitting    ) );
   IN.Add( Variable("DoseDistribution",     MATRIX2, &SIM.DoseDistribution ) );

   IN.Add( Variable("ImplantationArea",        REAL, &ImplantationArea     ) );
   IN.Add( Variable("ASCII",                BOOLEAN, &ascii                ) );

   IN.Add( Variable("DoseRate",    REAL, &DoseRate     ) );
 // 20100309 LuisJou MOD by J.H
   IN.Add( Variable("Extract",    INTEGER, &Extraccion     ) );
 //
   

#ifdef DAMAGE1D
#ifdef MODEL2
   Amorfizacion.CteA = 1.0;
   Amorfizacion.CteB = 100e-23;
   IN.Add( Variable("CteA",        REAL, &Amorfizacion.CteA     ) );
   IN.Add( Variable("CteB",        REAL, &Amorfizacion.CteB     ) );
   IN.Add( Variable("CteC",        REAL, &Amorfizacion.CteC     ) );
   IN.Add( Variable("CteD",        REAL, &Amorfizacion.CteD     ) );
   IN.Add( Variable("CteK",        REAL, &Amorfizacion.CteK     ) );

   IN.Add( Variable("AnchoCajas",  REAL, &Amorfizacion.ancho  ) );

#endif   
#endif
    
   IN.ReadParameters();
   if( Atoms.N == 0 )   // Atoms by default
   {
     AtomP = Atoms.Add( AtomDefinition( "B" , 11.000, 5, 519.0 ));
             Atoms.Add( AtomDefinition( "Si", 28.086,14, 519.0 ));
   };
   if(Atoms.N>N_ATOMS)
   {
    cerr << "Maximum number of atoms reached. Please recompile.\n"; 
    exit(MAX_ATOMS_ERROR);
   }
   if(ABC_Y==NullVector)
    {
     if(ABC_X.X==1 && ABC_X.Y==0 && ABC_X.Z==0) ABC_Y(0.0,1.0,0.0);
     else if(ABC_X.X==1 && ABC_X.Y==1 && ABC_X.Z==0) ABC_Y(0.0,0.0,1.0);
     else
      {
       cerr << "# Error: You must define de Y axe with FLAT token\n";
       exit(FLAT_NOT_DEFINED_ERROR);
      }
    }
   Temperature=SIM.Temperature;
   if( Therm ) Thermal.Generate( Atoms, Temperature );
   else        Thermal.Generate( Atoms, -1 );
   // Care! It uses the last Temperature value of the input file
   // The value shown in  output  file is not correct

// Dose Splitting preparation 
  List<SimulationDefinition> SIMULACIONES_2;
  LinkObject<SimulationDefinition>* Sim;
  SimulationDefinition Sim2;
  Sim = SIMULACIONES.PointerToFirstData();
  do
  {
   if(Sim->Data.DoseSplitting)
    {
     Sim2=Sim->Data;
     double TotalDose=Sim->Data.Dose, InitialAutomaticDose = 1e14;
     int i_=0;
     while(TotalDose >0 && i_<10)
      {
       if(Sim2.DoseDistribution[i_]==0)
        {
         Sim2.DoseDistribution[i_] = InitialAutomaticDose;
         TotalDose = TotalDose - Sim2.DoseDistribution[i_];
         InitialAutomaticDose = InitialAutomaticDose*10.0;
         if(InitialAutomaticDose>TotalDose) InitialAutomaticDose = TotalDose;
        }
       else
        {
         TotalDose = TotalDose - Sim2.DoseDistribution[i_];
         InitialAutomaticDose = TotalDose;
        }
       if(TotalDose < 0 )
          {
           cerr << "# Error: Dose distribution greater than specified Dose\n";
           exit(DOSE_DISTRIBUTION_ERROR);
          }
       Sim2.Dose = Sim2.DoseDistribution[i_];
       if(TotalDose==0) Sim2.DoseSplitting=0;
       i_++;
       SIMULACIONES_2.Add(Sim2);
      } 
    }
   else
    {
     SIMULACIONES_2.Add(Sim->Data);
    } 
   Sim = SIMULACIONES.PointerToNextData();
  } while(Sim!=NULL); 
    
 
  SIMULACIONES.DeleteAll();
  Sim = SIMULACIONES_2.PointerToFirstData();
  do
  {
   SIMULACIONES.Add(Sim->Data);
   Sim = SIMULACIONES_2.PointerToNextData();
  } while(Sim!=NULL);                                              

                                                       
   char ORDEN[300];
   unsigned j=0,k=0;
   for(unsigned i=0;i<strlen(argv[1]);i++)
    if(argv[1][i]=='/') j=i+1;
   for(unsigned i=j;i<strlen(argv[1]);i++) ORDEN[k++]=argv[1][i];
   ORDEN[k]=0;

   sprintf(DIRECTORY,"D%s", ORDEN );
#ifdef UNIX
   sprintf(ORDEN,"mkdir -p %s", DIRECTORY);
#else
    sprintf(ORDEN,"mkdir %s", DIRECTORY);
#endif

   if(system(ORDEN) == -1)
	   cerr << "Error issuing command " << ORDEN << endl;

   sprintf(ORDEN,"%s/Parameters.out", DIRECTORY );
   IN.Set(ORDEN,1);         // Write parameters used to 'output'
   IN.WriteParameters();
   
   HSTFile[ strlen(HSTFile)-1 ] ='\0';
}

 // Init a kind of things /////////////////////////////////////////////////

void MainLoop::Init( int argc, char **argv, int Initialize, int Pipe )
{
   // Select Electronic Stopping Function
   int i = 0;
   sFrenado[ strlen(sFrenado)-1 ] = '\0';
   while( strcmp( sFrenado, ESA[i].Name )!=0 )
    {
     if( strcmp( ESA[i].Name,"?" ) == 0 )
      {
       cerr << "# The stopping you can select are:\n";
       do
        {
         cerr << "# " << ESA[i].Name << endl;
         cerr << "#  " << ESA[i].Comment << endl;
         i--;
        } while( i>-1 );
       exit(ELECTRONIC_STOPPING_ERROR);
      };
     i++;
    };
   Frenado = i;

   if(Verbose)
   {
     cout << " Using " << ESA[i].Name << " for stopping. ";
     cout << ESA[i].Comment << _OK_ << endl;
   }
   // Select Ionization Function
   i = 0;
   sIonization[ strlen(sIonization)-1 ] = '\0';
   while( strcmp( sIonization, IS[i].Name )!=0 )
    {
     if( strcmp( IS[i].Name,"?" ) == 0 )
      {
       cerr << "# The ionization you can select are:\n";
       do
        {
         cerr << "# " << IS[i].Name << endl;
         cerr << "#  " << IS[i].Comment << endl;
         i--;
        } while( i>-1 );
       exit(IONIZATION_ERROR);
      };
     i++;
    };
   Ionizacion = i;

   if(Verbose)
   {
     cout << " Using " << IS[i].Name << " as ionization. ";
     cout << IS[i].Comment << _OK_ << endl;
   }

   GhiLimit = GhiLimit + DELTA1;

   XTalSize = 1.5;          // lattice units FIXED VALUE

   if( SimultaneousDistance < 0.0 ) SimultaneousDistance = 0.0;
   if( SimultaneousDistance2 < 0.0 ) SimultaneousDistance2 = 0.0;

   RND.Init(Seed1,Seed2);
   OR.Orienta( ABC_X, ABC_Y, Tha_Cut, Phi_Cut, Verbose );
   if( Initialize ) Solve.FirstInit( );
   SUPERLAYER.Init(XTalSize,InteractionRadius,GhiLimit,SimultaneousDistance);
   SUPERLAYER.MakeLayers( Bulk, Initialize, Verbose );
   Crystallite = SUPERLAYER.Layers.Search(1)->Data;
   LAYER       = Bulk.Search(1)->Data;
   LayerNumber = 0;
   Targets.Init(
          Crystallite.LatticeParameter,
          InteractionRadius,
          GhiLimit,
          ThresholdEnergy,
          SimultaneousDistance2, Temperature
        );
	
//   mean_cross_section = 1e-9;  // cm^2  20040521
   mean_cross_section = 0;  // cm^2  20040521	Esta es la buena 20040524
#ifdef DAMAGE1D  
   Amorfizacion.Mean_CrossSection = mean_cross_section; // Model W
   Amorfizacion.TotalGeneratedFrenkelPairs = 0.0;
   TMD = 0.0;
   SigmaMartinMean = 0.0;
#endif

 }
 
 // Init Simulation ///////////////////////////////////////////////////////
void MainLoop::InitSimulation( SimulationDefinition &SIM_ACTUAL, int ImplantNumber)
 {
    int i_;

    for(i_=0; i_<N_ATOMS; i_++ ) PROJS[i_]=SIM_ACTUAL.PROJS[i_];
    for(i_=0; i_<10; i_++ ) DoseDistribution[i_]=SIM_ACTUAL.DoseDistribution[i_];
    DoseSplitting    = SIM_ACTUAL.DoseSplitting;
    AtomP            = SIM_ACTUAL.AtomP;
    ENERGY           = SIM_ACTUAL.ENERGY;  
    Miller           = SIM_ACTUAL.Miller;
    DIROld           = SIM_ACTUAL.DIR;
    Tha              = SIM_ACTUAL.Tha;
    Phi              = SIM_ACTUAL.Phi;
    Divergency       = SIM_ACTUAL.Divergency;
    NumberOfImplants = SIM_ACTUAL.NumberOfImplants;
    Dose             = SIM_ACTUAL.Dose;
    Temperature      = SIM_ACTUAL.Temperature;
    rs0              = SIM_ACTUAL.RS0[0];
    for(int k=0;k<MAXLAYER;k++) RS0_LAYERS[k]=SIM_ACTUAL.RS0[k];

    if( !Miller )
    {
     DIROld = cos( Tha*M_PI/180 )                    *OR.EjeX +
              sin( Tha*M_PI/180 )*sin( Phi*M_PI/180 )*OR.EjeY +
              sin( Tha*M_PI/180 )*cos( Phi*M_PI/180 )*OR.EjeZ ;
    };

    Targets.DIROld = DIROld; // 20040521

    if( Therm ) Thermal.Generate( Atoms, Temperature );
    else        Thermal.Generate( Atoms, -1 );
   
    // Initialize rare-event algorithm
    if(PROJS[0]==0) PROJS[0]=AtomP;
    for(i_=0; i_<N_ATOMS; i_++)
     if(PROJS[i_]!=0) N_atoms=i_+1;
     else break;
          
    for(i_=0;i_<N_ATOMS;i_++)
     if(PROJS[i_]!=0)
     {
      HstPROJ  [PROJS[i_]].Initialize(MaxDepth, Dose );
      HstPROJ3D[PROJS[i_]].Initialize();
     }
    else break;

    pE = 0;   // Avoid activation of surface rare event

    TotalIterations   = 0;
    Interst           = 0;
    NTarg             = 0;
    ROT               = 0;
    TotalMeanDistance = 0;
    NSimMean          = 0.0;
    KK                = 1;

    WinAEff = WinA;
    WinBEff = WinB;

    // Solve WEffective
    WEffective = ImplantationArea * Dose /( NumberOfImplants * 1e16);

    #ifdef DAMAGE1D
     // Initialize Amorfizacion object variables
     Amorfizacion.Init( Dose, ImplantationArea, LevelOfSimulation,
               RecombinationFactor, KinchinPeaseConstant,
               AmorphizationDensity*1e-24,                // (defects/A^3)
               CutOffEnergy,                              // 
               DamageDisplacementEnergy,
               ImplantNumber-1, PreviousDamageFile, DoseRate, WEffective );
    #endif
    #ifdef DAMAGE3D  
     // Initialize Amorfizacion3D object variables
     Amorfizacion3D.Init3D( Dose, ImplantationArea, LevelOfSimulation,
               RecombinationFactor, KinchinPeaseConstant,
               AmorphizationDensity*1e-24,                // (defects/A^3)
               CutOffEnergy,                              // 
               DamageDisplacementEnergy,
               WinA, WinB,
               ImplantNumber-1, PreviousDamageFile );
     #ifndef DAMAGE3D_FASE1
      // Only ####### DAMAGE3D #######
      WinAEff = WinBEff = sqrt( ImplantationArea );
      if(WinA<WinAEff) WinAEff = WinA;
      if(WinB<WinBEff) WinBEff = WinB;   
     #endif   
    #endif


    if(Verbose)
    {
      cout << "\tDose             = " << Dose << endl;
      cout << "\tN                = " << NumberOfImplants << endl;
      cout << "\tImplantationArea = " << ImplantationArea << endl;
      cout << "\tWEffective       = " << WEffective << endl;
      cout << "\tNEffective       = " << NumberOfImplants * WEffective << endl;
    }
    cout << endl << _BLUE_;
    cout << "    Ion      Depth (nm)    Stat weigth      Path (nm)    Cascade time (sec.)\n";
    cout << " -----------------------------------------------------------------------------\n";
    cout << _NORMAL_;
    cout.setf(std::ios::fixed);
    cout.precision(4);
 } 

 // Do it /////////////////////////////////////////////////////////////////

void MainLoop::DoIt()
 {
  Interstitials.DeleteAll();
  Targets.Anteriores.DeleteAll();             

  DIR = DIROld;
  RND.MakeDivergency( DIR, DivType, Divergency );
  if ( Randomize )
   {
    Vector TTT;
//    POS = TTT(0 ,(RND.Generate()-0.5)*WinBEff,(RND.Generate()-0.5)*WinAEff);
    POS = TTT(0 ,(RND.Generate())*WinBEff,(RND.Generate())*WinAEff);
    POS = POS -  DIR;
   }
  Vector T, T2;
  T = T2 = POS % Crystallite.LatticeParameter;
  T.Integer2( Targets.TranslationParameter,Crystallite.LatticeParameter );
  T = T2 - T;

  double Energy, rnd, rnd2, s;
  if(Energy_Percentage!=0.0) Energy_Sigma = ENERGY*Energy_Percentage/100.0;
  switch( Energy_Spread)
  {
   case 1:  // Uniform distribution -------------
        rnd = RND.Generate();
        Energy = (rnd-0.5)*Energy_Sigma+ENERGY;
    break;
   case 2:  // Gaussian distribution ------------
        do {
        rnd  = 2*RND.Generate()-1;
        rnd2 = 2*RND.Generate()-1;
            s = rnd*rnd+rnd2*rnd2;
            } while( s>=1);
        Energy = ENERGY + Energy_Sigma*rnd*sqrt(-2.0*log(s)/s);
    break;
   default: // No distribution ------------------
        Energy = ENERGY;
    break;
  }
  
  // Simultaneous Multiple Implant =======================
  if(PROJS[0]==0)
   {
    CurrentProjectile = Projectile( Atom( AtomP, 0, T ), POS, GhiLimit,
                                    DIR, Energy, PROJECTILE );
    CurrentProjectile.Weigth = WEffective;                                    
    CurrentProjectile.InitialTime = 0.0; 
    CurrentProjectile.FinalTime   = 0.0;
    
    Projectiles.InsertAt( 1, CurrentProjectile );
    PROJS[0]=AtomP;
   }
  else
  {
   int i_;
   TotalMass=0.0;
   for(i_=0; i_<N_ATOMS; i_++)
    if(PROJS[i_]!=0)
     {
      // Calculate total mass involved to divide the energy
      N_atoms=i_+1;
      TotalMass+=Atoms.Search(PROJS[i_])->Data.W;
     }
    else break;
     
   for(i_=0; i_<N_atoms; i_++)
    { 
     AtomP = PROJS[i_];
     double TheEnergy = Energy * Atoms.Search(AtomP)->Data.W / TotalMass;
     CurrentProjectile = Projectile( Atom( AtomP, 0, T ), POS, GhiLimit,
                                    DIR, TheEnergy, PROJECTILE );
     CurrentProjectile.InitialTime = 0.0;				    
     CurrentProjectile.FinalTime   = 0.0;
     CurrentProjectile.Weigth = WEffective;
     Projectiles.InsertAt( 1, CurrentProjectile );
    }
  }
  
  if( (LAYER.Amorphous==1) )            // POLICRYSTALINE
   {
    double RotationMatrix[3][3];
    RND.MakeRotationMatrix( RotationMatrix );
    Vector R = CurrentProjectile.R;
    R = Crystallite.Nearest( ByDistance, R ).R;
    CurrentProjectile.R = CurrentProjectile.R - R;
    Projectiles.PointerToFirstData()->Data = CurrentProjectile;
    Crystallite.Rotate( RotationMatrix  );
   }

  // MAIN WHILE LOOP ///////////////////////////////////////////////////////////
  
  modif=0;
  while( Projectiles.N != 0 )
  {
   CurrentProjectile = Projectiles.PointerToFirstData()->Data;
   AtomP = CurrentProjectile.Index;

   if( CurrentProjectile.Status!=PROJECTILE && LevelOfSimulation==1 )
       {
        Projectiles.DeleteAll();
        #ifdef DAMAGE1D
          Crystallite.UnRotate();
        #endif
        #ifdef DAMAGE3D
          Crystallite.UnRotate();
        #endif
        break;                    // Exit from current implantation
       }

#ifdef DAMAGE1D
        double numero_aleatorio = RND.Generate();
        double posicionx        = OR.ProyectaX(CurrentProjectile.Abs);
        int numero_caja       = int (Amorfizacion.DefectDensity.GetIndex(posicionx));
    if(numero_caja>=0)
    {
         double nnnnn =  Amorfizacion.DefectDensity.M[numero_caja]
                        /Amorfizacion.Na;

         if( nnnnn > numero_aleatorio && IP<Crystallite.MeanAtomicRadius  ) 
         {   
          double RotationMatrix[3][3];
          RND.MakeRotationMatrix( RotationMatrix );
          if (NTarg!=0)  
           Crystallite.Near = Crystallite.Nearest( ByDistance, CurrentProjectile.R );
          CurrentProjectile.R = CurrentProjectile.R - Crystallite.Near.R;
          Projectiles.PointerToFirstData()->Data = CurrentProjectile;
          Crystallite.Rotate( RotationMatrix );
          ROT++;
          modif=1;
         }
        } 
#endif //DAMAGE1D
#ifdef DAMAGE3D
        double numero_aleatorio = RND.Generate();
        Vector posicion       = OR.Proyecta(CurrentProjectile.Abs);
        Vector numero_caja    = Amorfizacion3D.DefectDensity.GetIndex(posicion);
    if(numero_caja.X>=0)
    {
         double nnnnn =  
           Amorfizacion3D.DefectDensity.M[(int)numero_caja.X][(int)numero_caja.Y][(int)numero_caja.Z]
             /Amorfizacion3D.Na;

         if( nnnnn > numero_aleatorio && IP<Crystallite.MeanAtomicRadius  ) 
         {   
          double RotationMatrix[3][3];
          RND.MakeRotationMatrix( RotationMatrix );
          if (NTarg!=0)  
           Crystallite.Near = Crystallite.Nearest( ByDistance, CurrentProjectile.R );
          CurrentProjectile.R = CurrentProjectile.R - Crystallite.Near.R;
          Projectiles.PointerToFirstData()->Data = CurrentProjectile;
          Crystallite.Rotate( RotationMatrix );
          ROT++;
          modif=1;
         }
        } 
#endif //DAMAGE3D

     //
     // Multilayer section. Here we test the current coordinates of the
     //  projectile. If it is outside the current layer we update the
     //  Crystallite variable with the apropiate one content in SUPERLAYER
     //  variable. Also we update the Target's LatticeParameter variable.
     //
     Vector V = OR.Proyecta( Projectiles.PointerToFirstData()->Data.Abs );

     if( V.X > MaxDepth )
      {
       Projectiles.PointerToFirstData()->Data.Status = ION_NOT_FOLLOWED;
       Interstitials.Add( Projectiles.PointerToFirstData()->Data );

       Projectiles.DeleteAt( 1 );
      }

     if( V.X > SUPERLAYER.FrontSurface &&
         V.X < SUPERLAYER.BackSurface )
     { // First IF  // Inside the bulk
      if( (V.X < LAYER.XMin) || (V.X > LAYER.XMax) )
       {  // Second IF
      // The projectile is outside the current layer. We must search
      //    in the LayerDefinition list for another one layer
          LinkObject<LayerDefinition>* ItemLayer;
          LinkObject<Crystal>*   ItemCrystal;
          ItemLayer   =              Bulk.PointerToFirstData();
          LayerNumber = 0;
      ItemCrystal = SUPERLAYER.Layers.PointerToFirstData();
      do
        {
         if (V.X > ItemLayer->Data.XMin)
          if(V.X < ItemLayer->Data.XMax)
           {
         // The projectile is inside the ItemLayer
         // Change the layer
         LAYER       = ItemLayer->Data;
         Crystallite = ItemCrystal->Data;
         Targets.LatticeParameter     = Crystallite.LatticeParameter;
         Targets.TranslationParameter = Crystallite.LatticeParameter *
                        (double )TRANSLATION_PARAMETER;
         break;
           };
         ItemLayer   =      Bulk.PointerToNextData();
         ItemCrystal = SUPERLAYER.Layers.PointerToNextData();
         LayerNumber++;
        }
      while( ItemLayer != NULL );
      rs0 = RS0_LAYERS[LayerNumber];
      if( ItemLayer == NULL )
       {
          cerr << "MAIN.CPP: No layers found at " << V << endl;
          exit( NO_LAYERS_ERROR );
       };
       };       // End of second IF
     };         // End of first  IF

     if( (LAYER.Amorphous==2) && (IP<Crystallite.MeanAtomicRadius) )
      {
    double RotationMatrix[3][3];
    RND.MakeRotationMatrix( RotationMatrix );
    Vector R = CurrentProjectile.R;
    if (NTarg!=0) R = Crystallite.Nearest( ByDistance, R ).R;
        else          R = Crystallite.Near.R;
    CurrentProjectile.R = CurrentProjectile.R - R;
    Projectiles.PointerToFirstData()->Data = CurrentProjectile;
    Crystallite.Rotate( RotationMatrix );
    ROT++;
      }

      NTarg = Crystallite.SearchTargets(
                        Targets,
                        CurrentProjectile,
                        IP,
                        RND,
                        SUPERLAYER.FrontSurface,
                        SUPERLAYER.BackSurface,
                        Lista2
                      );
      Projectiles.PointerToFirstData()->Data = CurrentProjectile;

      int NS = 0;
      if( NTarg!=0 && Targets.PointerToFirstData()->Data.Index!=0 )
        NS = Targets.MakeCollision( Interstitials, Projectiles, Solve,
                           &Crystallite.DENS,&LAYER );
      else
        Targets.IncrementPosition( Interstitials, Projectiles,
                   SUPERLAYER.FrontSurface,
                   SUPERLAYER.BackSurface,
                   &Crystallite.DENS   );

      if(NS==-1)
        Targets.IncrementPosition( Interstitials, Projectiles,
                   SUPERLAYER.FrontSurface,
                   SUPERLAYER.BackSurface,
                   &Crystallite.DENS   );

      Targets.DeleteAll();  // Very important to recover memory

      if(    CurrentProjectile.Iterations > MaxItPerImplant
      || CurrentProjectile.Fly        > MaxFlyPerImplant )
       {
        switch( CurrentProjectile.Status )
         {
          case PROJECTILE :
            Projectiles.PointerToFirstData()->Data.Status = ION_NOT_FOLLOWED;
            break;
          default:
            Projectiles.PointerToFirstData()->Data.Status = NOT_FOLLOWED;
         };
        Interstitials.Add( Projectiles.PointerToFirstData()->Data );
        Projectiles.DeleteAt( 1 );
       };

#ifdef DAMAGE1D
   if(modif)
      {
        Crystallite.UnRotate();
        modif=0;
      }
#endif //DAMAGE1D
#ifdef DAMAGE3D
   if(modif)
      {
        Crystallite.UnRotate();
        modif=0;
      }
#endif //DAMAGE3D

   } // END OF WHILE LOOP //////////////////////////////////////////////////////////////////7

 }

 // Final stage of DoIt(). To divide the work in the paralell version
long unsigned int MainLoop::DoItRareEvent()
 {
  double Depth, Weigth;
#ifdef DAMAGE3D
  double min_InitPosY;
  double max_InitPosY;
  double min_InitPosZ;
  double max_InitPosZ;
  double yy,zz;
  Vector ThePOS;
  int jj;  
#endif  
  // 20100309 LuisJou MOD by J.H
  FILE *O_I,*O_V, *O_P, *O_PI, *O_BS, *O_MMONCA;
  double JX, JY, JZ, WW, EE;
  int II;
  if(Extraccion)
  {
   O_I  = fopen("Interstitials.txt","a+");    
   O_V  = fopen("Vacancies.txt","a+");    
   O_P  = fopen("Projectiles.txt","a+");    
   O_PI = fopen("Projectiles_InitPos.txt","a+");    
   O_BS = fopen("Backscattered.txt","a+");
   O_MMONCA = fopen("Cascades.txt", "a+");    
  }
  //

  Interst       += Interstitials.N;
  // Busco contar todas las iteraciones y distancias recorridas por los iones
  //  reales o virtuales
  if(Interstitials.N!=0)
  {
    LinkObject<Projectile>* III;

    double CascadeTime = 0.0;
    double Sum_r_by_t  = 0.0; // 20040521
    double Sum_t       = 0.0;
  
    III = Interstitials.PointerToFirstData();
    do
     {
      
      if(CascadeTime<III->Data.FinalTime) CascadeTime =III->Data.FinalTime;

      // 20040521
      Sum_r_by_t += ( III->Data.FinalTime - III->Data.InitialTime ) * III->Data.mean_r ;
      Sum_t      += ( III->Data.FinalTime - III->Data.InitialTime );
      if(max_r<III->Data.max_r) max_r = III->Data.max_r;

      TMD += III->Data.LocalPath;

//      printf("  mean_r = %g nm, time = %g sec.\n",  III->Data.mean_r/WEffective*0.1, 
//				( III->Data.FinalTime - III->Data.InitialTime ) );
      
//      printf("%s \t %12g %12g\n", AtomTypeStr[III->Data.Status],
//                    III->Data.InitialTime,III->Data.FinalTime);
//      printf("Mean R_i = %g\n", III->Data.mean_r/ WEffective );      
      
      switch( III->Data.Status )
       {
         case IMPLANTED_ION  :
         case BACK_SCATTERED :
         case ION_NOT_FOLLOWED :
         case OUTSIDE_ION    :
            TotalIterations   += III->Data.Iterations;
            TotalMeanDistance += III->Data.LocalPath;

            Depth = OR.ProyectaX( III->Data.Abs );
            Weigth = III->Data.Weigth;
            AtomP = III->Data.Index;
            HstPROJ  [AtomP].Add( Depth, Weigth, Depth );

 #ifdef DAMAGE3D
            min_InitPosY = -WinA/2-III->Data.InitPos.Y;
            max_InitPosY =  WinA/2-III->Data.InitPos.Y;
            min_InitPosZ = -WinB/2-III->Data.InitPos.Z;
            max_InitPosZ =  WinB/2-III->Data.InitPos.Z;
            
            ThePOS = OR.Proyecta( III->Data.Abs );
            Vector Despl;
            //1�.- Barrido por todas las zonas
            Despl.Y = -WinAEff*N_SECTOR/2;
            for(jj=0;jj<N_SECTOR;jj++)
            {
             Despl.Z = -WinBEff*N_SECTOR/2;
             for(int kk=0;kk<N_SECTOR;kk++)
             {
              //2�.- Validaci�n de la zona
              if(Despl.Y>=min_InitPosY && Despl.Y<=max_InitPosY && 
                 Despl.Z>=min_InitPosZ && Despl.Z<=max_InitPosZ)    
              {       
               //3�.- Actualizaci�n de los proyectiles por esa zona
               HstPROJ3D[AtomP].Add( ThePOS+Despl, Weigth, III->Data.LocalPath );
              }
              Despl.Z += WinAEff;
             }
             Despl.Y += WinBEff; 
            } 
 #else  // DAMAGE1D or NODAMAGE          
            HstPROJ3D[AtomP].Add( OR.Proyecta( III->Data.Abs ), Weigth, 
                                                   III->Data.LocalPath );
 #endif // DAMAGE3D

 if(KK%1==0)
  {
            cout << _BLUE_ << _BOLD_;
            cout.width(7);
            cout <<  KK <<_NORMAL_;     // Ion number
            cout.width(15);
            cout << Depth*0.1;      // Depth (nm)
            cout.width(15);
            cout.precision(10);
            cout << Weigth;         // Statistical weigth
            cout.precision(4);
            cout.width(15);
            cout << III->Data.LocalPath*0.1;// Path (nm)
            cout << " ";
        
	    printf(" %15g ", CascadeTime);  // Time (sec.)
        
	    if( III->Data.Status != IMPLANTED_ION )
             {
              cout << _RED_;
              cout << AtomTypeStr[ III->Data.Status  ];
              cout << _NORMAL_;
             }
           cout << endl;
		   		   
  }
    //
    // You can add as much output parameters as you want in this place
    //
#ifdef DAMAGE1D
#ifdef ENERGIA_DEPOSITADA
        Amorfizacion.Suma(III->Data.Abs.X,III->Data.Abs.Y,III->Data.Abs.Z,
              III->Data.Energy*III->Data.Weigth,"R");
#endif
#endif
#ifdef DAMAGE3D
#ifdef ENERGIA_DEPOSITADA
        Amorfizacion3D.Suma(III->Data.Abs.X,III->Data.Abs.Y,III->Data.Abs.Z,
              III->Data.Energy*III->Data.Weigth,"R");
#endif
#endif

        break;
    default:
#ifdef DAMAGE1D    
#ifdef ENERGIA_DEPOSITADA
            Depth = OR.ProyectaX( III->Data.Abs );
            Weigth = III->Data.Weigth;

        Amorfizacion.Suma(III->Data.Abs.X,III->Data.Abs.Y,III->Data.Abs.Z,
              III->Data.Energy*III->Data.Weigth,"R");
#endif
#endif
#ifdef DAMAGE3D    
#ifdef ENERGIA_DEPOSITADA
            Depth = OR.ProyectaX( III->Data.Abs );
            Weigth = III->Data.Weigth;

        Amorfizacion3D.Suma(III->Data.Abs.X,III->Data.Abs.Y,III->Data.Abs.Z,
              III->Data.Energy*III->Data.Weigth,"R");
#endif
#endif
        break;
       }; //Fin del switch

    // 20100309 LuisJou MOD by J.H
	if(Extraccion)
    {
     JX = OR.ProyectaX( III->Data.Abs );
     JY = OR.ProyectaY( III->Data.Abs );
	 JZ = OR.ProyectaZ( III->Data.Abs );
	 WW = III->Data.Weigth;
	 II = III->Data.Index;
	 EE = III->Data.Energy;
   	 if(III->Data.Status == INTERSTITIAL ) 
	 {
       	   if(Extraccion & 0x01) 
		    fprintf(O_I,"%+15lf %+15lf %+15lf %15lf %15lf %3d\n", JX, JY, JZ, WW, EE, II );
	   if(Extraccion & 0x20)
                fprintf(O_MMONCA, "I %f %f %f\n", JX/10, JY/10, JZ/10);

  	   JX = OR.ProyectaX( III->Data.InitPos );
	   JY = OR.ProyectaY( III->Data.InitPos );
	   JZ = OR.ProyectaZ( III->Data.InitPos );			  
	   if(Extraccion & 0x02) 
		    fprintf(O_V,"%+15lf %+15lf %+15lf %15lf %3d\n", JX, JY, JZ, WW, II );			  
	if(Extraccion & 0x20)
		fprintf(O_MMONCA, "V %f %f %f\n", JX/10, JY/10, JZ/10);	 
} 
	if(III->Data.Status == IMPLANTED_ION )
	 {
             if(Extraccion & 0x04)
	        fprintf(O_P,"%+15lf %+15lf %+15lf %15lf %15lf %3d\n", JX, JY, JZ, WW, EE, II );
  	     if(Extraccion & 0x20)
                fprintf(O_MMONCA, "B %f %f %f\n", JX/10, JY/10, JZ/10);

	   JX = OR.ProyectaX( III->Data.InitPos );
	   JY = OR.ProyectaY( III->Data.InitPos );
	   JZ = OR.ProyectaZ( III->Data.InitPos );			  
           if(Extraccion & 0x08)
	        fprintf(O_PI,"%+15lf %+15lf %+15lf %15lf %3d\n", 	JX, JY, JZ, WW, II );			   
	 }  
	if((III->Data.Status == BACK_SCATTERED)||(III->Data.Status == SPUTTERED) )
	 {
           if(Extraccion & 0x10)
            fprintf(O_BS,"%+15lf %+15lf %+15lf %+15lf %+15lf %+15lf %15lf %15lf %3d\n", 
				JX, JY, JZ, III->Data.Dir.X, III->Data.Dir.Y, III->Data.Dir.Z, WW, EE, II );
	 }  
	}
    //	   
	   
     }
    while( III = Interstitials.PointerToNextData(), III != NULL );
    
	// 20100309 LuisJou MOD by J.H
	if(Extraccion)
	{
	 fclose(O_I);
	 fclose(O_P);
	 fclose(O_V);
	 fclose(O_PI);
	 fclose(O_MMONCA);
	 fclose(O_BS);
	} 
	//
	
    
    // 20040521
    // Solve the final mean_r
    double mean_r;
    mean_r = ( Sum_r_by_t / Sum_t ) / WEffective; // In Angstroms
    
    //printf ( "Mean_R  = %g nm\n ", mean_r * 0.1);
    
//    fprintf(stderr, "%d %g", KK,  M_PI * mean_r * mean_r * 1e-16);			 
    // Solve the mean_cross_section
    mean_cross_section = (    
			      M_PI * mean_r * mean_r * 1e-16 // In cm^2
			    + (KK-1)*mean_cross_section
			  )/KK;
#ifdef DAMAGE1D			  
   Amorfizacion.Mean_CrossSection = mean_cross_section; // Model W			  

//    fprintf(stderr, " %g\n", mean_cross_section );			 
//    fprintf(stderr, "Max_r = %g, Sigma_max = %g \n", max_r, M_PI* max_r*max_r* 1e-16);
   
//  printf("TotalGeneratedFrenkelPairs = %g\n", Amorfizacion.TotalGeneratedFrenkelPairs);    

    SigmaMartinMean = (Amorfizacion.TotalGeneratedFrenkelPairs/(5e22*TMD*1e-8)
     			    + (KK-1)*SigmaMartinMean  )/KK;
			  
    Amorfizacion.TotalGeneratedFrenkelPairs = 0.0;
#endif			  	
    TMD = 0.0;
  };

#ifdef DAMAGE1D
 #ifdef MODEL2
  FILE *OutputFile;
  char PP[300];
  sprintf(PP,"%s/DDBC_1D.dat",DIRECTORY);
  OutputFile=fopen(PP,"a");
//  OutputFile=fopen("DDBC_1D.dat","a");
  fprintf(OutputFile,"%lf\n", Amorfizacion.DeltaDamageByCascade);
  Amorfizacion.DeltaDamageByCascade = 0;
  fclose(OutputFile);
 #else
 Amorfizacion.UpdateDefectDensity();
 #endif
#endif  // DAMAGE1D

#ifdef DAMAGE3D
 #ifdef MODEL2
  FILE *OutputFile;
  OutputFile=fopen("DDBC.dat","a");
  fprintf(OutputFile,"%lf\n", Amorfizacion3D.DeltaDamageByCascade);
  Amorfizacion3D.DeltaDamageByCascade = 0;
  fclose(OutputFile);
 #else
  Amorfizacion3D.UpdateDefectDensity3D(Interstitials.PointerToFirstData()->Data.InitPos);
 #endif
#endif // DAMAGE3D

  KK = KK + 1;
  return (KK);
 }

 // Done //////////////////////////////////////////////////////////////////

void MainLoop::Done(int KK_Sim)
 {
  TotalMeanDistance = TotalMeanDistance / TotalIterations; 
  NSimMean = NSimMean / TotalIterations;
  double EffDens =NSimMean/(TotalMeanDistance*M_PI
          *InteractionRadius*InteractionRadius)*1e24;
  if(Verbose)
  {
    cout.setf(std::ios::scientific);
    cout << endl
         << "Mean Distance     = " << TotalMeanDistance << " (A)" << endl;
    cout << "Effective density = " << EffDens
         << " at/cm^3 (aprox.)" << endl;
    cout << "NSimMean          = " << NSimMean << " at/col" << endl;
    cout << "NColMean          = " << (double )TotalIterations/NumberOfImplants
         << " col/cascade" << endl;
    cout << "RotationMean      = " << (double )ROT/TotalIterations
         << " rot/total_iterations" << endl << endl;
    cout << "Simultaneous col. = " << PEPE
         << " vs " << TotalIterations << " total collisions ("
         << (double )PEPE/TotalIterations*100.0 << " %)" << endl;
  }

//  fprintf(stderr, " Final (N = %d) cross section %g\n", KK, mean_cross_section );			 
  
  /*
  FILE * OO;
  OO = fopen("CrossSection_vs_Energy","a+");
  
  fprintf(OO,"%g %g\n",ENERGY, mean_cross_section );
  
  fclose(OO);
 */ 
  
#ifdef DAMAGE1D



  Amorfizacion.WriteAmorphization();
#endif
#ifdef DAMAGE3D
  Amorfizacion3D.WriteAmorphization3D(0);
#endif

  //
  // Rare event. Show histogram section.
  //
  char sTemp[300];
  sprintf(sTemp,"%s%02d",HSTFile,KK_Sim); // Number simulation
  for(int i_=0; i_<N_atoms; i_++)
    {
     AtomP=PROJS[i_];
     HstPROJ[AtomP].Do();
     if(PearsonIV) HstPROJ[AtomP].SolvePearsonIV(AtomP,HSTFile);
     HstPROJ[AtomP].SaveHisto(AtomP, 0, sTemp,
                                 FullDoseOutput,NumberOfImplants,
                                 NumberOfImplants);
     if(Verbose) HstPROJ[AtomP].ShowFinalInfo(AtomP);

     HstPROJ3D[AtomP].Solve3D(AtomP,Dose,WinA*WinB, ascii);
     if(ascii) HstPROJ3D[AtomP].Show(AtomP,Dose);
    };
      
  STPO.Clear();         // Clear current stopping tables and recover memory
   
  if(PAUSE) cin >> PAUSE;

#ifdef MODEL2

  printf(" Next parameter is valid only if RareEvent is deactivated\n"); 
  printf("\nSigma_0 = %g cm^2\n\n", SigmaMartinMean);

#endif

 }

void MainLoop::JoinDoseSplitting()
{
  VarHst A[N_ATOMS],B[N_ATOMS];
  ifstream ifsA,ifsB;
  ofstream ofsC;
  char cc, Nombre[300];
  int nIons, nIV,i, KK_Sim, KK, i_Atom, theAtom;
  double nIR;
  double depth,conc,conc1,conc2, depth1[N_ATOMS],depth2[N_ATOMS];
  double outside_conc;

  LinkObject<SimulationDefinition>* Sim;
  KK_Sim=1;
  KK = 1;
  Sim = SIMULACIONES.PointerToFirstData();
  do
  {
   // I read the first result
   for(i_Atom=0; i_Atom<N_ATOMS; i_Atom++) 
    if(Sim->Data.PROJS[i_Atom]!=0)
     {
      theAtom=Sim->Data.PROJS[i_Atom];
      sprintf(Nombre,"%s/.%s%02d_%02d",DIRECTORY,HSTFile,KK_Sim,theAtom);  
//      cout << "Reading ... " << Nombre << endl;
      ifsA.open(Nombre);
      ifsA >> cc >> nIons >> nIR >> nIV >> outside_conc;    
      i=0;
      while(!ifsA.eof()&&i<N_BIN)
      {
       ifsA >> depth >> conc >> conc1 >> conc2;
       A[theAtom].M[i++]=conc;
      }
//      cout << "Profile 1 (Prof= "<< depth << " ) N="<< i-1  <<endl;        
      depth1[theAtom]=depth*N_BINS/(N_BIN-1);
      assert(depth1[theAtom]>0);
      A[theAtom].Set(0.0,depth1[theAtom]);
      ifsA.close();
     }
     
   if(Sim->Data.DoseSplitting!=0)
    {
     do
     {
      KK_Sim++;
      Sim = SIMULACIONES.PointerToNextData();

      for(i_Atom=0; i_Atom<N_ATOMS; i_Atom++) 
       if(Sim->Data.PROJS[i_Atom]!=0)
        {
         theAtom=Sim->Data.PROJS[i_Atom];
         sprintf(Nombre,"%s/.%s%02d_%02d",DIRECTORY,HSTFile,KK_Sim,theAtom);   
         //cout << "Reading ... " << Nombre << endl;
         ifsB.open(Nombre);
         ifsB >> cc >> nIons >> nIR >> nIV >> outside_conc;
         i=0;
         while(!ifsB.eof()&&i<=N_BIN)
         {
          ifsB >> depth >> conc >> conc1 >> conc2;
          B[theAtom].M[i++]=conc; 
         }
//         cout << "Perfil 2 (Prof= "<< depth << " ) N="<< i-1  <<endl;       
         depth2[theAtom]=depth*N_BINS/(N_BIN-1);
         assert(depth2[theAtom]>0);
         B[theAtom].Set(0.0,depth2[theAtom]);
         ifsB.close();  

         if(depth1[theAtom]>depth2[theAtom]) 
          {
            depth2[theAtom]=depth1[theAtom];
            B[theAtom].ReScale(0.0,depth1[theAtom]);
          }
         else                               
          { 
            depth1[theAtom]=depth2[theAtom];
            A[theAtom].ReScale(0.0,depth2[theAtom]);
          }

         // Suma
         for(i=0;i<N_BINS;i++) A[theAtom].M[i] += B[theAtom].M[i];
        }
      if(!Sim) 
       {
        cerr << "# Error en MAINLOOP.CPP: JoinSplitting()" <<endl;
        exit(NO_SIM_ERROR);
       }
     }while(Sim->Data.DoseSplitting!=0 );
    }
   else KK_Sim++;

   // Write
   for(i_Atom=0; i_Atom<N_ATOMS; i_Atom++) 
    if(Sim->Data.PROJS[i_Atom]!=0)
     {
      theAtom=Sim->Data.PROJS[i_Atom];
      sprintf(Nombre,"%s/%s%02d_%02d",DIRECTORY,HSTFile,KK,theAtom);   
//      cout << "Writting ... " << Nombre << endl;
      int N_Bin;
      for(N_Bin=N_BINS;N_Bin>0;N_Bin--)
       if(A[theAtom].M[N_Bin]!=0) break;
      N_Bin++; 
      
      ofsC.open(Nombre);
      ofsC << cc << " "<< nIons << " " << nIR << " " << nIV 
                << " " << outside_conc << endl;
      for(i=0;i<N_Bin;i++) 
        ofsC << i*A[theAtom].mAncho << " " << A[theAtom].M[i] << endl;
      ofsC << "# The two first numbers after the # means:\n";
      ofsC << "#  Trajectories, double ions = NumberOfImplants*EffectiveWeigth, Virtual ions\n";
      ofsC << "# Next numbers are:\n";
      ofsC << "#  Depth(nm), Conc. (at/cm^3)\n";      
      ofsC.close(); 
     }
   KK++; 
   Sim = SIMULACIONES.PointerToNextData();    
  } while(Sim!=NULL); 
    
}
}
