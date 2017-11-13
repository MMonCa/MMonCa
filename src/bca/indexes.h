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
// INDEXES.H
//
///////////////////////////////////////////////////////////////////////////
#ifndef _INCLUDE_INDEXES_H_
#define _INCLUDE_INDEXES_H_

#include <math.h>
#include <assert.h>
#include <iostream>

#define NV3D      7

const int N_ALPHA   = 1;
const int N_PHI     = 1;
const int N_THETA   = 1;
const int N_S       = 30;
const int N_ENERGY  = 72;

namespace BCA {
extern double E_GRID[N_ENERGY];
extern double S_GRID[N_S];

// Get the lower index for this variable and his value /////////////////////
inline int GetEnergyIndex( double Value )
    {
       int Index;
       int a = 0;
       int b = N_ENERGY;

       do
       {
        Index = (a+b)/2;
        if( E_GRID[Index] == Value ){ a = Index ;break; };
        if( E_GRID[Index] > Value ) b = Index;
        else            a = Index;

       } while ( b-a!=1 );
       Index = a;
       if( Index< 0       ) Index =       0;
       return( Index ); // Lower index
    }

inline int GetSIndex(double Value )
    {
       int Index;

       int a = 0;
       int b = N_S;
       do
       {
        Index = (a+b)/2;
        if( S_GRID[Index] == Value ){ a = Index ;break; };
        if( S_GRID[Index] > Value ) b = Index;
        else            a = Index;
       } while ( b-a!=1 );
       Index = a;
       if( Index< 0       ) Index =       0;
       return( Index ); // Lower index
    }

// Get the value for this variable and his index //////////////////////////
inline void GetEnergyValue( int Index, double &Value, double &NextValue)
    {
       assert( Index>=0&&Index<=N_ENERGY-1 );

       Value     = E_GRID[ Index     ];
       NextValue = E_GRID[ Index + 1 ];
    }

 inline void GetSValue( int Index, double &Value, double &NextValue)
    {
       assert( Index>=0&&Index<=N_S-1 );

       Value     = S_GRID[ Index     ];
       NextValue = S_GRID[ Index + 1 ];

    }

// Get the lower index for this variable and his value ////////////////////
 inline int GetAlphaIndex( double Value )
    {
       int Index;

       while( Value > 2*M_PI ) Value -= 2*M_PI;
       while( Value < 0.0    ) Value += 2*M_PI;

       Value = Value*N_ALPHA/M_PI/2;
       Index = (int) Value;

       return( Index ); // Lower index
    }

inline int GetPhiIndex( double Value )
    {
       int Index;

       while( Value > 2*M_PI ) Value -= 2*M_PI;
       while( Value < 0.0    ) Value += 2*M_PI;

       Value = Value*N_PHI/M_PI/2;
       Index = (int) Value;

       return( Index ); // Lower index
    }

inline int GetThetaIndex( double Value )
    {
       int Index;
       Value = Value*N_THETA/M_PI;
       Index = (int) Value;

       if( Index< 0       ) Index =         0;
       if( Index> N_THETA  -2 ) Index = N_THETA    -2;

       return( Index ); // Lower index
    }

// Get the value for this variable and his index //////////////////////////
inline void GetAlphaValue( int Index, double &Value, double &NextValue )
    {
       assert( Index>=0&&Index<=N_ALPHA-1 );

       Value     = Index*2*M_PI/N_ALPHA;

       if( Index==N_ALPHA-1 )
        NextValue = 0.0;
       else
        NextValue = Value + 2*M_PI/N_ALPHA;
    }

 inline void GetPhiValue( int Index, double &Value, double &NextValue )
    {
       assert( Index>=0&&Index<=N_PHI-1 );

       Value     = Index*2*M_PI/N_PHI;
       if( Index==N_PHI )
        NextValue = 0.0;
       else
        NextValue = Value + 2*M_PI/N_PHI;
    }

inline void GetThetaValue( int Index, double &Value, double &NextValue )
    {
       assert( Index>=0&&Index<=N_THETA-1 );

       Value     = Index*M_PI/N_THETA;
       NextValue = Value + M_PI/N_THETA;
    }
}
#endif  // _INCLUDE_INDEXES_H_
