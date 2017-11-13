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
//////////////////////////////////////////////////////////////////////////
//
// THERMAL.H
//
///////////////////////////////////////////////////////////////////////////

#include "thermal.h"

// Generate probability table & Amplitude displacements
namespace BCA {
void ThermClass::Generate( List<AtomDefinition> &Atoms, double T )
{
   double CVAMP[10]= { 2.7777777777777778E-02, - 2.7777777777777778E-04,
             4.7241118669690098E-06, - 9.1857730746619636E-08,
             1.8978869988970999E-09, - 4.0647616451442255E-11,
             8.9216910204564526E-13, - 1.9939295860721076E-14,
             4.5189800296199182E-16, - 1.0356517612182470E-17  };

   double SCALE = 12.063464; // (A)

   // Construct the probability table                   OK.
   double QZ,QA,QB,QX, QQ;
   int  I, J;
   LinkObject<AtomDefinition> *ATOM;


   QZ = 1.0 / ( 2*THERMAL_NTBL+1 );
   QA = 0.5;
   Table[0] = 0.0;
   for( J = 1; J < THERMAL_NTBL; J++ )
   {
    QA = QA - QZ;
    QB = - log( QA );
    QX = sqrt( 2*QB );
    QB = ( QX * (QX * 0.010328 + 0.802853) + 2.515517 ) /
     ( QX * (QX * (QX * 0.001308 + 0.189269) + 1.432788) + 1.0 );

    Table[J] = QX - QB;
   };

   // Evaluate the RMS Displacement Amplitudes              OK.

   assert( Atoms.NumberOfElements() > 0 ); // For debugging

   ATOM = Atoms.PointerToFirstData();

   for( J=1;
    J<=Atoms.NumberOfElements();
    J++, ATOM=Atoms.PointerToNextData() )

    if( ATOM->Data.TDebye > 0.0 )
     {
      if( T<=0 )  // Negative or zero ambient temperature
       {
    if( T==0 ) QX = 0.25;
    else       QX = 0.00;
       }
      else    // Positive ambient temperature
       {
    QZ = ATOM->Data.TDebye / T;
    QA = QZ * QZ;

    if( QZ <= 2.16 )    // Convergent expansion of Debye Function
     {
       QB = 1.0 / QZ;
       QX = QB;
       for( I=0; I<10; I++)
        {
         QB = QB * QA;
         QX = QX + QB * CVAMP[I];
        };
     }
    else            // Asymptotic expansion of Debye Function
     {
       QB = QZ;
       QX = 0.25 + M_PI * M_PI / (6.0 * QA );
       QA = exp( - QZ );
       QQ = QA;
       for( I=0; I<10; I++)
        {
         QX = QX - ( (1.0 + 1.0 / QB) / QB ) * QQ;
         QB = QB + QZ;
         if( QB > 160.0 ) break;
         QQ = QA * QQ;
        };
     };
       };

      // Update Temperature Amplitude                   OK.

      QQ = SCALE * sqrt( QX/ (ATOM->Data.W*ATOM->Data.TDebye) );

      ATOM->Data.TAmp = QQ;

      TAmp[J]=QQ;   // For fast access
     };
}

///////////////////////////////////////////////////////////////////////////

ThermClass Thermal;
}
