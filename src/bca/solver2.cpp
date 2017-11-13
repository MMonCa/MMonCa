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
// SOLVER2.CPP
//
///////////////////////////////////////////////////////////////////////////
#include "solver2.h"
namespace BCA {
///////////////////////////////////////////////////////////////////////////
//
// Init function
//
void Solver::FirstInit( )
{
 Table.ReadTable( );
}

/////////////////////////////////////////////////////////////////////////
//
//  Solve all calculations

void Solver::All( Projectile &P,Projectile &T,double &AdvanceMinimum,
        double &Einel, double &QElect )
{
 Q=0;
 AtomP   = P.Index;
 AtomT   = T.Index;
 RT  = T.R;             // Target Position  (A)
 RP  = P.R;             // Projectile Position  (A) 
 VP  = P.Dir;
 EP  = P.Energy;            // Initial Energy   (eV)
 DeltaR  = RT - RP;
 Ghi     = VP * DeltaR;             // Forward impact param (A)
 p   = sqrt( DeltaR*DeltaR - Ghi*Ghi ); // Impact parameter (A)

 double x1, x2, cos2, e1, e2, Advance, sin2, ratioTE, qElect;
 Vector NewVT;
 if( Ghi<=0 ) Ghi = 1e-5;

// Solve  --------------------------------------------------------- 
 Table.InterpolateAll_( AtomP, AtomT, EP, p, x1, x2, cos2, e1,e2, Q, qElect );

// Post-solver treatment -----------------------------------------
 if ( x2  < DELTA1 )   x2 = DELTA1;
 if (cos2 <   -1.0 ) cos2 = -1.0;
 if (cos2 >    1.0 ) cos2 =  1.0;
 sin2 = sqrt( 1-cos2*cos2 );
 Advance = Ghi + x1;
 if( Advance < AdvanceMinimum ) AdvanceMinimum = Advance;
 if ( p > 1e-7 )
  {
   double C2 = sin2 / p;
   NewVT = (cos2-Ghi*C2)*VP + C2*DeltaR;
  }
 else NewVT = VP;       // Frontal collision
 ratioTE = e2/EP;
 NewVT = (double ) sqrt( ratioTE*Atoms.Search( AtomT )->Data.W
               /Atoms.Search( AtomP )->Data.W) * NewVT;
 NewVT = Precision( NewVT, DELTA1, 0.0 );
 T.Energy   = ratioTE;          // ratioTE => T.Energy
 T.Ghi      = Advance;          // Store temporary in this field
 T.S        = p;
 T.R        = RT + x2 * VP;
 T.Dir      = NewVT;
 P.Energy   = e1;
 P.R        = RP + Advance * VP;

 P.S = p; 

 T.Einel    = Q ;
 Einel += Q;
 QElect+=qElect;
}
/////////////////////////////////////////////////////////////////////////
}
