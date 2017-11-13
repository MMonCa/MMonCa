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
// COL3D.H
//
///////////////////////////////////////////////////////////////////////////
#ifndef _INCLUDE_COL3D_H_
#define _INCLUDE_COL3D_H_

#include "defs.h"
#include "vector.h"
#include "edt3d.h"
#include "inel.h"
#include "indexes.h"
#include "table3d.h"
#include "ssf.h"

#include <time.h>
namespace BCA {
///////////////////////////////////////////////////////////////////////////
// Class Trajectorie Calculus 
class TC
{
 private:
  double Z1, Z2, M1, M2, S_3D, MD[3][3], MI[3][3], au, h,vx1,vy1,vx2,vy2,
         vz1,vz2,rx1,ry1,rx2,ry2,rz1,rz2,ax1,ay1,ax2,ay2,az1,az2,ax10,ay10,
         az10,ax20,ay20,az20,R,hh,h2h,KL;
  double q;      // Non local inelastic stopping // Perdidas frenado inelastico       
  double qElect; // Local inelastic stopping //Perdidas frenado electronico
  EDT    EDTable;
  
 public:
  float  *TABLE;
  std::ofstream OUT3D;
  int    Out3D;
  SSF    Screen;

void Init( double Z1_= 5.0, double Z2_= 14.0, double M1_= 11.0,
       double M2_= 28.086, const std::string &fn = "", int size = N_ENERGY,
       int Initialize = 1, int SaveIt = 1 ); // (YES)  NO = 0 YES = Other
       
void MakeTrajectorie( int PAINT = 1 );

void SetParamPSystem( double ENERGY, double Xx, double Yy, double Zz );

void MakeMatrix( double sinTHETA, double cosTHETA,
             double sinPHI, double cosPHI,
               double sinALPHA, double cosALPHA );

double GetResultPSystem( double &x1, double &x2, double &cos1, double &cos2,
             double &E1, double &E2, double &Q, double &QElect );

void SolveAll(
    double ENERGY, double S, double Ghi, double THETA, double PHI, double ALPHA,
    double &x1, double &x2, double &cos2, double &E1, double &E2, double &Q,
    double &QElect );

void MAKEOUT( int i1=0, int i2=0, int i3=0, int i4=0, int i5=0, int Print = 0);

void SAVE_ASCII( char* pfnt, int size=N_ENERGY );

void SAVE_BIN( char *pfnt, int size=N_ENERGY );

void DOIT_PVM( int i_ENERGY, int i_PHI, int i_E = 0, int i_P = 0 );

void DOIT( int i_ENERGY );

void DOITTEST( int i_ENERGY = 0, int i_PHI = 0, int i_THETA  = 0,
           int i_ALPHA  = 0, int i_S   = 0, int Print = 0 );
};
}
#endif  // _INCLUDE_COL3D_H_
