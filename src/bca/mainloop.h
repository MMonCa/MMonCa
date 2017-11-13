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
// MAINLOOP.H
//
///////////////////////////////////////////////////////////////////////////
#ifndef _INCLUDE_MAINLOOP_H_
#define _INCLUDE_MAINLOOP_H_

#include "list.h"
#include "layer.h"
#include "simul.h"
#include "solver2.h"
#include "collisio.h"
#include "histo.h"
#include "read.h"
namespace BCA {
extern List<LayerDefinition> Bulk;  // Define all layers
extern List<SimulationDefinition> SIMULACIONES; // Define all simulations
extern StoppingStruct ESA[ ];
extern IonizationStruct IS[ ];
extern double rs0;
extern double RS0_LAYERS[MAXLAYER];
extern int Frenado;
extern int Ionizacion;
extern char AtomTypeStr[][20] ;
extern int LayerNumber;
extern int PROJS[N_ATOMS];
///////////////////////////////////////////////////////////////////////////

class MainLoop  
{
 public:
    double InteractionRadius,  XTalSize, ThresholdEnergy,
         GhiLimit, SimultaneousDistance, SimultaneousDistance2, Tha,
         Phi, Divergency, ENERGY, Dose,
         WinA, WinB, IP, Temperature,MaxDepth, pE, pX,
         Tha_Cut, Phi_Cut, Energy_Sigma, Energy_Percentage,
             RecombinationFactor, KinchinPeaseConstant, AmorphizationDensity,
             CutOffEnergy, DamageDisplacementEnergy,
             ImplantationArea, WEffective, WinAEff, WinBEff;
    double DoseRate;
    double mean_r, mean_cross_section, max_r;
    double DoseDistribution[10];
      Vector POS, DIR, DIROld, ABC_X, ABC_Y;
     int N_atoms;
     double TotalMass;
     int DivType, Miller, Randomize, Seed1, Seed2, Verbose,
         LevelOfSimulation, AtomP, Therm, NTarg, PearsonIV,
         modif, Energy_Spread, FullDoseOutput, DoseSplitting, ascii;
     char sFrenado[80], sIonization[80], HSTFile[256], PreviousDamageFile[300];
 unsigned
        long NumberOfImplants, MaxItPerImplant, MaxFlyPerImplant, KK, ROT,
             TotalIterations, Interst;
 long double TotalMeanDistance;
 double TMD, SigmaMartinMean;
 
 // 20100309 LuisJou MOD by J.H
	int Extraccion;
 //
      Random RND;
      Solver Solve;
      Target Targets;
  LayerClass SUPERLAYER;
     Crystal Crystallite;
  Projectile CurrentProjectile;
 LayerDefinition LAYER;
 List<Projectile> Interstitials, Projectiles;
 List<Atom> Lista2;
 std::ofstream OS;
     int PAUSE;

 // Assign the default value to the variables
 void Default();

 // Read and set the parameters
 void SetParameters( int argc, char **argv );

 // Init a kind of things
 void Init( int argc, char **argv, int Initialize=1, int Pipe=1 );
 
 // Init Simulation
 void InitSimulation( SimulationDefinition &SIM_ACTUAL, int ImplantNumber=0 );

 // Do it
 void DoIt();

 // Last stage of DoIt()
 long unsigned int DoItRareEvent();

 // Done
 void Done(int KK_Sim=1);
 
 // Join the splitted doses
 void JoinDoseSplitting();
};
}
// End of MainLoop
#endif  // _INCLUDE_MAINLOOP_H_
