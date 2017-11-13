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
// RANDOM.H
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

// Define Random CLASS ////////////////////////////////////////////////////
#ifndef _INCLUDE_RANDOM_H_
#define _INCLUDE_RANDOM_H_
namespace BCA {
class Random
{
  private:
    unsigned int _i;
    unsigned int _j;
    long int _carry;
    unsigned long int _u[97];
    
  public:
     Random( int Seed1=1802, int Seed2=9373 )
      {
       Init(Seed1,Seed2);
      }
     void   Init( int Seed1, int Seed2 );
     double Generate();
     int    Test();

     int    MakeRotationMatrix( double M[3][3] );
     int    MakeDivergency( Vector &DIR, int DivType, double Divergency );
};


// Generate a random number
inline double Random::Generate()
{
    long int delta = _u[_i]-_u[_j];
    
    if(delta<0) delta += 16777216L;
    _u[_i]=delta;
    if(_i==0) _i=96;
    else      _i--;
    if(_j==0) _j=96;
    else      _j--;
    _carry -= 7654321L;
    if(_carry<0) _carry+= 16777213L;
    delta -= _carry;
    if(delta<0) delta+= 16777216L;
    
    return delta/16777216.0;
}
}
#endif  // _INCLUDE_RANDOM_H_
