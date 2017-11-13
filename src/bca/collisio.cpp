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
#include "collisio.h"
#include "mainloop.h"

using std::cerr;

namespace BCA {

EDT1D DENS;

extern char AtomTypeStr[][20];


extern MainLoop MAINLOOP; // 20040521


 // Init procedure ////////////////////////////////////////////////////////

void Target::Init( Vector LP, double IR, double GL, double TE, double SD, double Temp )
 {
  if(LP==NullVector) LatticeParameter(5.431,5.431,5.431);
  else               LatticeParameter     = LP;
  TranslationParameter = (double )TRANSLATION_PARAMETER*LatticeParameter;
  InteractionRadius    = IR;
  GhiLimit         = GL + DELTA1;
  ThresholdEnergy      = TE;
  SimultaneousDistance = SD;
  Temperature = Temp;
 }

//=========================================================================

// Make a collision with each target of the list //////////////////////////
 
int Target::MakeCollision( List<Projectile> &Interstitials,
                    List<Projectile> &Projectiles, 
                    Solver &Solve,      
            DensityClass *DENSITY, LayerDefinition *LAYER  )
{
 double RatioET, RatioEP, Einel;
 Vector Advance, VFinal;
 int NSim;
 P = Projectiles.PointerToFirstData()->Data;

 NSim = Phase1( Solve, RatioEP, RatioET, Einel, VFinal, Advance ); 

 if( NSim == 0 ) return -1;
 NSimMean = NSimMean + NSim;
 Einel = INELASTIC_LOSSES( DENSITY, RatioEP*P.Energy, Advance, VFinal);
 double Coef_Energy = ScaleEnergies(  RatioEP, RatioET, Einel ); 
 UpdateInterstitials( Interstitials, Projectiles, Coef_Energy, LAYER );
 UpdateProjectile( Interstitials, Projectiles, Advance, VFinal );
 Results.DeleteAll();
#ifdef DAMAGE1D
#ifdef ENERGIA_DEPOSITADA
 Amorfizacion.Suma(NewP.Abs.X,NewP.Abs.Y,NewP.Abs.Z,
        (Coef_Energy*Einel)*P.Weigth,"Q");
#endif
#endif
#ifdef DAMAGE3D
#ifdef ENERGIA_DEPOSITADA
 Amorfizacion3D.Suma(NewP.Abs.X,NewP.Abs.Y,NewP.Abs.Z,
        (Coef_Energy*Einel)*P.Weigth,"Q");
#endif
#endif
 return NSim;
}

// PHASE 1. Solve the position and final direction of proj. ///////////////

int Target::Phase1( Solver &Solve, double &RatioEP, double &RatioET, double &Einel,
            Vector &VFinal, Vector &AdvanceV )
{
 LinkObject<Projectile>* TheTarget, * Item2; 
 LinkObject<Vector>    * ItemV;
 Projectile NewT, T; 
 double MinimunAdvance, QElect = 0;
 int i1, i2, NSim = 0;
 
 // Simultaneous collision correction. Avoid recollide with the same targets.
 Vector Ref = P.Abs - P.R;

 if(Anteriores.N!=0) 
 {
  Item = PointerToFirstData();
  for( i1=0; i1<N; i1++ )
  {
   ItemV = Anteriores.PointerToFirstData();
   for( i2=0; i2<Anteriores.N; i2++)
   {
    if( ItemV->Data == (Item->Data.InitPos + Ref) )
     {
       PEPE ++;
       Delete( Item );   
     }
    ItemV = Anteriores.PointerToNextData();
   }
   Item = PointerToNextData();
  }
 }
 if( N==0 ) return 0; //////////// There are no targets ///////////

 // Main LOOP for simultaneous collisions -------------------------------
 Einel =0;
 MinimunAdvance  = 1e10;
 TheTarget = PointerToFirstData();
 do
 {
  T = NewT = TheTarget->Data;
  Vector T_GetR = NewT.R;
  NewT.InitPos  = Ref + T_GetR;
  NewP = P;
  Solve.All( NewP, NewT, MinimunAdvance, Einel, QElect );
  NewT.Abs  = Ref + NewT.R;
  T.Einel   = NewP.Energy;
  TheTarget->Data = T;

  Results.Add( NewT );

  NSim ++;
  TheTarget = PointerToNextData();
 }
 while ( NSim < N );

 // SECOND Pass. -----------------------------------------------------------
 // Repeat the 'simultaneity' test if apropiate

 Anteriores.DeleteAll();              
 if( NSim>=1 )
  {
   MinimunAdvance = MinimunAdvance + SimultaneousDistance;
   MinimunAdvance = MinimunAdvance + fabs( MinimunAdvance*DELTA1 );
   double Advance = -1e10;
   RatioET =  0.0;
   Einel   =  0.0;
   VFinal(0.0,0.0,0.0);
   Item  = Results.PointerToFirstData();
   Item2 =         PointerToFirstData();
   do
    {
      if( Item->Data.Ghi>=MinimunAdvance )
       {
        Item->Data.Index = 0 ;
        NSim--;
       }
      else 
       {
        Anteriores.Add( Item2->Data.InitPos + Ref );
        if( Advance < Item->Data.Ghi ) Advance = Item->Data.Ghi;
        RatioET = RatioET + Item->Data.Energy;   
        VFinal  = VFinal  + Item->Data.Dir;
        Einel   = Einel   + Item->Data.Einel; 
       }   
      Item  = Results.PointerToNextData();
      Item2 =         PointerToNextData();          
    }
   while ( Item!=NULL );

   VFinal  = Precision( P.Dir - VFinal, DELTA1, 0.0 );
   RatioEP = VFinal * VFinal;
   VFinal  = Unitary( VFinal );
   AdvanceV = Advance * P.Dir;
   return NSim;
  }
 else // NO TARGET ATOMS IN THE CURRENT GROUP OF SITES.
  {
   cerr << "# Error: Target::Phase1, No targets atoms\n";
   exit(NOATOMS_ERROR);
  };
}

//-------------------------------------------------------------------------
double Target::ElectronicLosses( DensityClass *DENSITY, double NewEnergy, Vector Distance )
{
 int Z2 = 0; // Para el calculo del frenado en vuelo libre
 AtomDefinition AD = Atoms.Search( P.Index )->Data;
 return DENSITY->LossesNew( AD.Z, AD.W, P.Abs,
          Unitary( Distance ), P.Energy, NewEnergy, (double )Distance, Z2 );
}

//-------------------------------------------------------------------------
double Target::INELASTIC_LOSSES( DensityClass *DENSITY, double NewEnergy,
                                Vector Distance, Vector DirFinal  )
{
 int Z2 = 0;
 double Inel = 0;  
 double rDistance=(double ) Distance;
 Distance = Unitary(Distance);
 AtomDefinition AD = Atoms.Search( P.Index )->Data;

 // LOCAL Interaction (Projectile's Electrons - Target's Electrons) =========
 Item = Results.PointerToFirstData();
 do{
    if( Item->Data.Index>0 )
    {
      AtomDefinition AD2 = Atoms.Search( Item->Data.Index )->Data;
      if(Z2==0)Z2=AD2.Z;
      Inel += InelLosses( AD.Z, AD.W, P.Abs,
                  Distance, P.Energy, NewEnergy, rDistance, 
                  AD2.Z, Item->Data.InitPos, Item->Data.Abs, Item->Data.Dir, 
                  DirFinal );
    };
    Item = Results.PointerToNextData();
 }while(Item!=NULL);

 // NON LOCAL Interaction (Projectile's Nucleus - Target's Electrons) =======
 Inel += DENSITY->LossesNew( AD.Z, AD.W, P.Abs, Distance, P.Energy, NewEnergy, 
                                    rDistance, Z2 );
 return Inel;
}

//-------------------------------------------------------------------------

double Target::ScaleEnergies( double RatioEnergyProjectiles,
            double RatioEnergyTargets,
            double InelasticEnergy )
{
 double EP0, Coef, Coef2, Tope, EnergyTargets, EnergyProjectile;

 EP0   = P.Energy;
 Coef  = 1.0 / 
        ( RatioEnergyProjectiles + RatioEnergyTargets + InelasticEnergy/EP0 );
 EnergyTargets    = Coef * RatioEnergyTargets * EP0;
 InelasticEnergy  = Coef * InelasticEnergy; 
 Tope = EnergyTargets + InelasticEnergy;
 if( EP0 < Tope )
  {
   Coef2            = EP0/Tope;
   //EnergyTargets    = Coef2 * EnergyTargets;
   //InelasticEnergy  = Coef2 * InelasticEnergy;
   Coef             = Coef2 * Coef;
   EnergyProjectile = 0.0;   
  }
 else
  { 
   EnergyProjectile = EP0 - Tope;
  };
 NewP.Energy = EnergyProjectile;
 return Coef;
}

//-------------------------------------------------------------------------

void Target::UpdateInterstitials( List<Projectile> &Interstitials,
              List<Projectile> &Projectiles,
              double Coef_Energy,
              LayerDefinition *LAYER  )
{
 Vector Translation;

 Item = Results.PointerToFirstData();
 do
 {
  if( Item->Data.Index != 0 ) // If NOT OUTSIDE
   {
    Item->Data.Energy = Coef_Energy * Item->Data.Energy*P.Energy;
    Item->Data.Dir    = Unitary(Item->Data.Dir);

    double BindingEnergy;
    if( Item->Data.Index == Item->Data.Index1)
     BindingEnergy = LAYER->XTal.Search( Item->Data.LS )
                                     ->Data.BindingEnergy;
    else
     BindingEnergy = LAYER->XTal.Search( Item->Data.LS )
                                     ->Data.BindingEnergy2;

    if ( Item->Data.Energy > BindingEnergy + ThresholdEnergy )
     {
      // Update the relative coordinates. Move into the proper suboctant
      Translation = Item->Data.R;
      Translation.Integer2( TranslationParameter, LatticeParameter );

      Item->Data.R = Item->Data.R -Translation;
      Item->Data.Energy  -= BindingEnergy;
      Item->Data.Ghi      = GhiLimit;
      Item->Data.Status   = SECONDARY_PROJECTILE;
      
      // The statistical weight of the secondary recoils
      //   is the same as the primary projectile 20000602
      Item->Data.Weigth = P.Weigth;   
      
      Item->Data.InitialTime = P.FinalTime;
      Item->Data.FinalTime   = P.FinalTime;
      Item->Data.mean_r      = 0.0;	// 20040521
      Item->Data.max_r       = 0.0;

#ifdef DAMAGE1D
      #ifdef ENERGIA_DEPOSITADA
      Amorfizacion.Suma(
          OR.ProyectaX(Item->Data.Abs),
          OR.ProyectaY(Item->Data.Abs),
          OR.ProyectaZ(Item->Data.Abs),
        (Item->Data.Energy+BindingEnergy)*P.Weigth,"B", P.Weigth);
      #endif
      #ifdef MODEL2
      Amorfizacion.SeqDamageAccumulation( 
          OR.ProyectaX(NewP.Abs), 
        (Item->Data.Energy+BindingEnergy)*P.Weigth,'B', P.Weigth, Temperature );
      #else  
      Amorfizacion.DamageAccumulation( 
          OR.ProyectaX(NewP.Abs), 
        (Item->Data.Energy+BindingEnergy)*P.Weigth,'B', P.Weigth);
      #endif	
#endif

#ifdef DAMAGE3D
      #ifdef ENERGIA_DEPOSITADA
      Amorfizacion3D.Suma(
          OR.ProyectaX(Item->Data.Abs),
          OR.ProyectaY(Item->Data.Abs),
          OR.ProyectaZ(Item->Data.Abs),
        (Item->Data.Energy+BindingEnergy)*P.Weigth,"B", P.Weigth);
      #endif
 #ifdef MODEL2
      Amorfizacion3D.SeqDamageAccumulation3D(OR.Proyecta(NewP.Abs),
        (Item->Data.Energy+BindingEnergy)*P.Weigth,'B', P.Weigth, Temperature);
 #else
      Amorfizacion3D.DamageAccumulation3D( 
          OR.Proyecta(NewP.Abs),
        (Item->Data.Energy+BindingEnergy)*P.Weigth,'B', P.Weigth);
 #endif
#endif


      Projectiles.Add( Item->Data );
     }
    else
    {
#ifdef DAMAGE1D
      #ifdef ENERGIA_DEPOSITADA
      Amorfizacion.Suma(Item->Data.Abs.X,Item->Data.Abs.Y,Item->Data.Abs.Z,
        Item->Data.Energy*P.Weigth,"I");
      #endif 
      #ifdef MODEL2
      Amorfizacion.SeqDamageAccumulation( OR.ProyectaX(NewP.Abs),
            Item->Data.Energy*P.Weigth,'I',P.Weigth, Temperature );
      #else
      Amorfizacion.DamageAccumulation( OR.ProyectaX(NewP.Abs),
            Item->Data.Energy*P.Weigth,'I',P.Weigth );
      #endif
#endif

#ifdef DAMAGE3D
      #ifdef ENERGIA_DEPOSITADA
      Amorfizacion3D.Suma(Item->Data.Abs.X,Item->Data.Abs.Y,Item->Data.Abs.Z,
        Item->Data.Energy*P.Weigth,"I");
      #endif
 #ifdef MODEL2
      Amorfizacion3D.SeqDamageAccumulation3D(OR.Proyecta(NewP.Abs),
        (Item->Data.Energy+BindingEnergy)*P.Weigth,'B', P.Weigth, Temperature);
 #else
      Amorfizacion3D.DamageAccumulation3D( OR.Proyecta(NewP.Abs),
          Item->Data.Energy*P.Weigth,'I',P.Weigth );
 #endif
#endif
    }; 

   };
  Item = Results.PointerToNextData();
 } while( Item != NULL );

}

//-------------------------------------------------------------------------

void Target::UpdateProjectile( List<Projectile> &Interstitials,
               List<Projectile> &Projectiles, 
               Vector Advance, Vector FinalDirection )
{
 double XIM, GhiMin, GhiMax;
 Vector Translation;

 // FIND MINIMUM FORWARD SEARCH DISTANCE ON NEW PATH
 // Solve XIM parameter. Correction (back-scattered collisions) for the
 //   GhiMax parameter.
 
 GhiMin  = P.Ghi;
 if ( Advance*P.Dir >= GhiMin ) XIM = 0.0;
 else   XIM = ( GhiMin - Advance*P.Dir )*( FinalDirection *P.Dir );
 if (XIM<0.0) XIM = 0.0;

 NewP.Dir = FinalDirection;
 NewP.R   = P.R + Advance;
 NewP.Abs = P.Abs + Advance;
 Translation = NewP.R;
 Translation.Integer2( TranslationParameter, LatticeParameter );
 NewP.R  = NewP.R - Translation;

 // CALCULATE MINIMUM SEARCH DISTANCE FOR THE NEXT COLLISION

 GhiMax = -2.0; 
 Item = Results.PointerToFirstData();
 do
 {
  if( Item->Data.Index > 0 ) // If NOT OUTSIDE
  {
    // The new Ghi parameter
    double Ghi = NewP.Dir * ( Item->Data.InitPos - NewP.Abs );
    if (Ghi>GhiMax) GhiMax  = Ghi;
  };
  Item = Results.PointerToNextData();
 } while (Item!=NULL);

 if ( GhiMax < GhiLimit ) GhiMax = GhiLimit;
 GhiMax = GhiMax + XIM;
 NewP.Ghi = GhiMax;
 NewP.Fly = 0;
 NewP.Iterations++;
 NewP.LocalPath += (double ) Advance;
 
// if(isnan(P.Energy)) {P.Energy=0;NewP.Energy=0;}; // REVISAR <------------------

 double delta_time;
 // 7.21958e-15 = 1e-10/sqrt(2*1.602e-19/1.67e-27 ), conversion to S.I.
 AtomDefinition AD = Atoms.Search( P.Index )->Data;
 if(P.Energy > 0) delta_time = (double ) Advance / sqrt( P.Energy/ AD.W ) * 7.21958e-15;
 else             delta_time = 0.0;
 
 NewP.FinalTime += delta_time ;
 
 // 20040521
 Vector Rel           = (NewP.Abs - NewP.InitPos);
 Vector Perpendicular = Rel  -  (Rel * DIROld) * DIROld;
 double  rrr = (double ) Perpendicular;
 
 if( NewP.FinalTime !=NewP.InitialTime )
 NewP.mean_r = 1/(NewP.FinalTime-NewP.InitialTime)
               *(  rrr * NewP.Weigth * delta_time  
	         + (NewP.FinalTime-NewP.InitialTime - delta_time)*NewP.mean_r
		);
 if( NewP.max_r<rrr) NewP.max_r = rrr;

// if(isnan(NewP.mean_r)) printf( " TFinal = %g  TInicial = %g  rrr = %g  W = %g  deltaT = %g \n", 
//           NewP.FinalTime, NewP.InitialTime, rrr, NewP.Weigth, delta_time );
// cerr << NewP.FinalTime<< " " << NewP.mean_r/0.20 <<  " " << rrr << " " << delta_time << endl;
 
 
 // Add current projectile to interstitials list --------------------------

 if ( NewP.Energy <= ThresholdEnergy )
   {
     switch ( NewP.Status )
     {
      case PROJECTILE   : NewP.Status = IMPLANTED_ION;
              break;
      case SECONDARY_PROJECTILE : NewP.Status = INTERSTITIAL;
              break;
      default       : cerr << "# Error: Target::Makecollision, Wrong status " 
                           << AtomTypeStr[NewP.Status] << " " 
                           << NewP.Status << "\n";
              break;
     };
     Interstitials.Add( NewP );     // ADD interstitial
     Projectiles.DeleteAt( 1 );     // Delete it as projectile
   }

 // ELSE Add current projectile to projectiles list -----------------------

 else Projectiles.PointerToFirstData()->Data = NewP;
}

///////////////////////////////////////////////////////////////////////////

void Target::IncrementPosition( List<Projectile> &Interstitials,
                        List<Projectile> &Projectiles, 
                    double FrontSurface, double BackSurface, DensityClass *DENSITY )
{
 // Variables
 Projectile  P;
 double OtherGhi, XPosition;
 Vector Forward, Translation;

 P = Projectiles.PointerToFirstData()->Data;
 OtherGhi = P.Ghi;
 if(OtherGhi<GhiLimit) OtherGhi = GhiLimit;
 P.Ghi   = OtherGhi;
 Forward = P.Dir*OtherGhi;
 P.Abs   = P.Abs + Forward;
 P.R     = P.R   + Forward ;
 P.Fly++;
 P.LocalPath += (double ) Forward;
 Translation =  P.R;
 Translation.Integer2( TranslationParameter, LatticeParameter );
 P.R = P.R - Translation;

 
// if(isnan(P.Energy)) P.Energy=0; // REVISAR <-----------------
 
 double delta_time;
 AtomDefinition AD = Atoms.Search( P.Index )->Data;
 if(P.Energy > 0)  delta_time = (double ) Forward / sqrt( P.Energy/ AD.W ) * 7.21958e-15;
 else              delta_time = 0.0;
 P.FinalTime += delta_time;
 
 // 20040521
 Vector Rel           = (NewP.Abs - NewP.InitPos);
 Vector Perpendicular = Rel  -  (Rel * DIROld) * DIROld;
 double  rrr = (double ) Perpendicular;

 if( NewP.FinalTime !=NewP.InitialTime )
 NewP.mean_r = 1/(NewP.FinalTime-NewP.InitialTime)
               *(  rrr * NewP.Weigth * delta_time  
	         + (NewP.FinalTime-NewP.InitialTime - delta_time)*NewP.mean_r
		);
 if( NewP.max_r<rrr) NewP.max_r = rrr;


 // Test if the projectile is in the Xtal. FRONT surface ------------------ 
 XPosition = OR.ProyectaX( P.Abs );
 Projectiles.PointerToFirstData()->Data = P;
 
                            // Solo si ya entro en
                            // el cristal
 
 if( (-XPosition > FrontSurface + InteractionRadius) && (P.S!=0) )
 {
    switch( P.Status )
    {
      case PROJECTILE    : P.Status = BACK_SCATTERED;
               break;
      case SECONDARY_PROJECTILE : P.Status = SPUTTERED;
               break; 
      default        : cerr << "# Error: Target::IncrementPosition1, Wrong status "
                            << AtomTypeStr[NewP.Status] << " " 
                            << NewP.Status << "\n";
    };  
   Interstitials.Add( P );
   Projectiles.DeleteAt( 1 );
 }
 else  // Test if the projectile is in the Xtal. BACK surface -------------------
 if( XPosition > BackSurface + InteractionRadius )
 {
    switch( P.Status )
    {
      case PROJECTILE : P.Status = OUTSIDE_ION;
               break;
      case SECONDARY_PROJECTILE    : P.Status = OUTSIDE;
               break;
      default        : cerr << "Error: Target::IncrementPosition2, Wrong status "
                            << AtomTypeStr[NewP.Status] << " " 
                            << NewP.Status << "\n";
    };
   Interstitials.Add( P );
   Projectiles.DeleteAt( 1 );
 }
 else // The projectile is inside the Xtal
 { 
   // We stop electronically the ion
   double Einel =  ElectronicLosses(DENSITY,P.Energy,Forward);
   P.Energy = P.Energy - Einel;

   if ( P.Energy <= ThresholdEnergy )
   {
     switch ( P.Status )
     {
      case SECONDARY_PROJECTILE   : P.Status = INTERSTITIAL;
              break;
      case PROJECTILE : P.Status = IMPLANTED_ION;
                  break;
      default       : cerr << "# Error: Target::IncrementPosition, Wrong status "
                            << AtomTypeStr[NewP.Status] << " " 
                            << NewP.Status << "\n";
              break;
     };
     Interstitials.Add( P );     // ADD interstitial
     Projectiles.DeleteAt( 1 );     // Delete like projectile
   }
   else    
   Projectiles.PointerToFirstData()->Data = P; 
      
#ifdef DAMAGE1D
#ifdef ENERGIA_DEPOSITADA
 Amorfizacion.Suma(P.Abs.X,P.Abs.Y,P.Abs.Z,Einel*P.Weigth,"F");
#endif
#endif
#ifdef DAMAGE3D
#ifdef ENERGIA_DEPOSITADA
 Amorfizacion3D.Suma(P.Abs.X,P.Abs.Y,P.Abs.Z,Einel*P.Weigth,"F");
#endif
#endif
 }
}
}
///////////////////////////////////////////////////////////////////////////
