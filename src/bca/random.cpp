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
// RANDOM.CPP
//
///////////////////////////////////////////////////////////////////////////

//
//  Implement the algorithm originally appeared in "Toward a Universal
// Random Number Generator" by George Marsaglia and Arif Zaman.
// Florida State University Report: FSU-SCRI-87-50 (1987)
//
//  It was later modified by F. James and published in " A Review of
// Pseudorandom Number Generators"
//
//  It passes ALL of the tests for random number generators and has a
// period of 2^144, is completely portable (gives bit identical results
// on all machines with at least 24-bit mantissas in the floating point
// representation).
//
//  The algorithm is a combination of a Fibonacci sequence (with lags of 97
// and 33, and operation "substraction plus one, modulo one") and an
// "arithmetic sequence" (using substraction).
//

#include "vector.h"
#include "random.h"

using std::cerr;
using std::endl;

namespace BCA {
// Initialize Random Number Generator
void Random::Init( int Seed1 = 1802, int Seed2 = 9373 )
{
  if ((Seed1<0)||(Seed1>31328)||(Seed2<0)||(Seed2>30081))
  {
    cerr << "\nThe first SEED must have a value beetwen 0 and 31328"
         << "\nThe second SEED beetwen 0 and 30081\n";
  };
  unsigned long int S1 = Seed1;
  unsigned long int S2 = Seed2;
 
  int i=((S1/177)%177)+2;
  int j=( S1 %177)+2;
  int k=((S2/169)%178)+1;
  int l=  S2 %169;
  int a,b;

  for (a=0;a<97;a++)
  {
    unsigned long int sum = 0;
    unsigned long int t   = 16777216L;

    for (b=0;b<24;b++)
    {
      unsigned long int m=(((i*j)%179)*k)%179;
      i=j;
      j=k;
      k=m;
      l=(53*l+1)%169;
      t >>=1;
      if (((l*m)%64)>=32) sum += t;
    };
    _u[a] = sum;
  };
  
  _i = 96;
  _j = 32;
  _carry = 362436L;
}

// Test for Generator
int Random::Test()
{
 long int Values[6] = { 6533892L,14220222L, 7275067L,
            6172232L, 8354498L,10633180L };
 int i;
 Random RND(1802,9373);     // Initialize

 for( i = 0; i < 20000; i++ ) RND.Generate();
 for( i = 0; i < 6; i++ )
  {
   if (4096.0*4096.0*RND.Generate()-Values[i]==0)
    cerr << "TEST " << i << " OK" << endl;
   else return -1; 
  } 
 return 0;
}

// Generate a random rotation matrix
int Random::MakeRotationMatrix( double M[3][3] )
{
 double QA,QB,QC,QD,QE,QF,QG,QH, QQ;

 do
 {
  QB = Generate();
  QD = QB*QB;
  QC = Generate();
  QE = QC*QC;
  QF = QD+QE;
 } while( QF>1.0 );
 QD = (QD-QE)/QF;
 QC = 2*QB*QC/QF;
 if( Generate() > 0.5 ) QC = -QC;

 QA = Generate();
 if( Generate() > 0.5 ) QA = -QA;
 QB = sqrt( 1-QA*QA );

 do
 {
  QE = Generate();
  QG = QE*QE;
  QF = Generate();
  QH = QF*QF;
  QQ = QH+QG;
 } while( QQ>1.0 );
 QE = 2*QE*QF/QQ;
 QF = (QG-QH)/QQ;
 if( Generate() > 0.5 ) QE = -QE;

 M[0][0] = QA*QD*QF-QC*QE;
 M[1][0] = QA*QC*QF+QD*QE;
 M[2][0] =-QB*QF;

 M[0][1] =-QC*QF-QA*QD*QE;
 M[1][1] = QD*QF-QA*QC*QE;
 M[2][1] = QB*QE;

 M[0][2] = QB*QD;
 M[1][2] = QB*QC;
 M[2][2] = QA;

 return 0;
}

// Generate a random divergency

int Random::MakeDivergency( Vector &DIR, int DivType, double Divergency )
{

 int L;
 double Beam[3], PCos[3], DVY[3];
 double DivCos, DivSin, RND1, Rnd1, RND2, Rnd2, R1, R2, R3, Sum;

//19990917

if(DivType==1 || DivType==2) Divergency=Divergency*0.707106781;

////

 Beam[0] = DIR.Z;
 Beam[1] = DIR.Y;
 Beam[2] = DIR.X;

 if( 1.0 - fabs( Beam[2] ) < 1e-6 ) L=0;
 else                   L=2;

 DivCos = cos( Divergency*M_PI/180 );
 DivSin = 1.0 - DivCos;

 do
 {
  RND1 = Generate();
  Rnd1 = RND1*RND1;

  RND2 = Generate();
  Rnd2 = RND2*RND2;

  Sum  = Rnd1 + Rnd2;
 }
 while( Sum > 1.0 );

 R1 = ( Rnd1 - Rnd2 )/Sum;
 R2 = ( 2*RND1*RND2 )/Sum;

 if( Generate() > 0.5 ) R2 = -R2;

 switch (DivType)
 {
  case 0: // Isotropically distributed
      R3 = Generate();
      break;
  case 1: // Non uniform cosine distribution to azimuthal angle
      do { R3 = M_PI*( 2*Generate() - 1 ); }
      while( cos( R3 ) + 1 < 2*Generate() );
      R3 =cos( R3 );
      break;
  default: // No distribution
      R3 = 0.0;
 };

 R3 = DivCos + DivSin * R3;

 PCos[L] = 0;

 if( L < 2 )
  {
   PCos[1]= Beam[2];
   PCos[2]=-Beam[1];
  }
 else
  {
   PCos[0]= Beam[1];
   PCos[1]=-Beam[0];
  };

 DVY[0] = Beam[L]*Beam[0];
 DVY[1] = Beam[L]*Beam[1];
 DVY[2] = Beam[L]*Beam[2];

 DVY[L] = DVY[L] - 1.0;

 double Q, QA, QB;

 QA = PCos[0]*PCos[0] + PCos[1]*PCos[1] + PCos[2]*PCos[2];
 QB = DVY[0]*DVY[0] + DVY[1]*DVY[1] + DVY[2]*DVY[2];
 Q  = 1.0 - R3*R3;

 R1 = R1 * sqrt( Q/QA );
 R2 = R2 * sqrt( Q/QB );

 PCos[0] = R1*PCos[0] + R2*DVY[0] + R3*Beam[0];
 PCos[1] = R1*PCos[1] + R2*DVY[1] + R3*Beam[1];
 PCos[2] = R1*PCos[2] + R2*DVY[2] + R3*Beam[2];

 DIR( PCos[2], PCos[1], PCos[0] );

 DIR = Unitary( DIR );

 return DivType;
}
}
///////////////////////////////////////////////////////////////////////////
