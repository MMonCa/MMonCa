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
// COLLISIO.H
//
///////////////////////////////////////////////////////////////////////////
#ifndef _INCLUDE_COLLISIO_H_
#define _INCLUDE_COLLISIO_H_

#include <iomanip>
#include "defs.h"
#include "inel.h"
#include "atom.h"
#include "list.h"
#include "solver2.h"
#include "edt3d.h"
#include "layer.h"

#ifdef DAMAGE1D
#include "amorf.h"
namespace BCA {
extern clase_Amorfizacion Amorfizacion; }
#endif

#ifdef DAMAGE3D
#include "amorf3d.h"
namespace BCA {
extern clase_Amorfizacion3D Amorfizacion3D; }
#endif

namespace BCA {

extern long double NSimMean;
extern long PEPE;

///////////////////////////////////////////////////////////////////////////

class Body 
{
 public:
  Vector R; // Position
  Vector V; // Velocity
  Vector A; // Acceleration
  Vector A0;    // Last Acceleration
  double   M; // Mass
  int    Z; // Atomic number
};

///////////////////////////////////////////////////////////////////////////

class Target : public List<Projectile>
{
 private:
  // Variables ------------------------------------------------------------
  double InteractionRadius, GhiLimit, ThresholdEnergy, SimultaneousDistance;
  LinkObject<Projectile> *Item;
  Body T[BODY_MAX]; // Targets
  double Temperature;

 public:
  int  nT;      // Number of targets
  Vector Ref_Crystal; 
  List<Projectile> Results;
  Body PBody;       // Incident projectile
  List<Vector> Anteriores;  // Abs & freezed coordinates of previous targets
  Projectile NewP, P;
  Vector TranslationParameter,LatticeParameter;
  Vector DIROld;

//-------------------------------------------------------------------------
 void Init( Vector LP=NullVector, double IR=2.7155, double GL=0.1, double TE=10, 
             double SD=0.25, double Temp=300 );

//-------------------------------------------------------------------------
 int MakeCollision( List<Projectile> &Interstitials,
                    List<Projectile> &Projectiles, 
                    Solver &Solve,      
            DensityClass *DENSITY, LayerDefinition *LAYER );

//-------------------------------------------------------------------------
 int Phase1( Solver &Solve, double &RatioEP, double &RatioET, double &Einel,
            Vector &VFinal, Vector &AdvanceV );

//-------------------------------------------------------------------------
 double ElectronicLosses( DensityClass *DENSITY, double NewEnergy, Vector Distance );

//-------------------------------------------------------------------------
 double INELASTIC_LOSSES( DensityClass *DENSITY, double NewEnergy,
                                Vector Distance, Vector DirFinal  );

//-------------------------------------------------------------------------
 double ScaleEnergies( double RatioEnergyProjectiles,
            double RatioEnergyTargets,
            double InelasticEnergy );

//-------------------------------------------------------------------------
 void UpdateInterstitials( List<Projectile> &Interstitials,
              List<Projectile> &Projectiles,
              double Coef_Energy,
              LayerDefinition *LAYER  );

//-------------------------------------------------------------------------
 void UpdateProjectile( List<Projectile> &Interstitials,
               List<Projectile> &Projectiles, 
               Vector Advance, Vector FinalDirection );

//-------------------------------------------------------------------------
 void IncrementPosition( List<Projectile> &Interstitials,
                        List<Projectile> &Projectiles, 
                    double FrontSurface, double BackSurface, DensityClass *DENSITY );

}; // End of Target class
}
///////////////////////////////////////////////////////////////////////////
#endif  // _INCLUDE_COLLISIO_H_
