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
///////////////////////////////////////////////////////////////////////////f
//
//  INEL.H
//  Inelastic stopping due to electron-electron interaction
//
///////////////////////////////////////////////////////////////////////////
#ifndef _INCLUDE_INEL_H_
#define _INCLUDE_INEL_H_

#include "defs.h"
#include "ssf.h"
#include <math.h>

#define NTABLE 1000
#define DLIM     15.0   // (A)
#define MAX_ATOM 60

//extern char TABLEDIRDEFAULT[300];
namespace BCA {
class InelT
{
 public:
     double * TABLA[MAX_ATOM][MAX_ATOM]; // Pointer table to stopping tables
     SSF Ssf;
     double au;

 InelT()
  {
   int i,j;
   for(i=0; i<MAX_ATOM; i++)
    for( j=0; j<MAX_ATOM; j++)
     TABLA[i][j]=NULL;
  }

 inline double II( double R  )
  {
   double X = Ssf.ScreeningOnly( R*au, X );
   return X*X/R; 
  }

 inline double I( double R )
  {
   double a = 0;
   double x = R;
   double w = 0.02;

   if(R>=30) return 1e-10;
   for( x = R; x < 30; x+=w )
    a = a + II(x)*w;
   return a; 
  }

  double FFF(  double R, double Z1, double Z2 );

};

// ------------------------------------------------------------------------
// Integration in a trajectory
//

#define NPASOS_ 6
#define Vcritica 0.7

double InelLosses( double Z1, double M1, Vector PAbs, Vector Dir,
           double ENERGY, double E1, double AllowedDistance,
           double Z2, Vector TAbs, Vector TAbs1, Vector TDir1,
           Vector DirFinal );
}
#endif    // _INCLUDE_INEL_H_
