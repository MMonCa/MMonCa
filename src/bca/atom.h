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
// ATOM.H
//
///////////////////////////////////////////////////////////////////////////
#ifndef _INCLUDE_ATOM_H_
#define _INCLUDE_ATOM_H_

#include "defs.h"
#include "vector.h"
#include "random.h"
#include "thermal.h"
#include <assert.h>
#include <string.h>
#include <iostream>

namespace BCA {

// Enumerate the status type
extern Vector NullVector;
extern Vector XVector;
extern ThermClass Thermal;

enum AtomType
    {
        LATTICE_ATOM    = 0,
        PROJECTILE,
        VACANCY,
        INTERSTITIAL,
        SECONDARY_PROJECTILE,
        SPUTTERED,
        BACK_SCATTERED,
        IMPLANTED_ION,
        OUTSIDE_ION,
        ION_NOT_FOLLOWED,
        NOT_FOLLOWED,
        OUTSIDE,
        PHONON_EXCITED
    };

// Define the Atom class //////////////////////////////////////////////////

class Atom
{
  public:
    unsigned char Index;    // Active Index
    unsigned char LS;       // Lattice site type
    Vector R;               // Position
    unsigned char Index1;   // First atom index
    unsigned char Index2;   // Second atom index for the same lattice site
    double   X;               // Probability of the first atom index

    // Constructor
    Atom( int IC=0, int L=0, Vector RC=NullVector, int IC2=0, double xx=1.0 );
    // Compare two atoms
    int operator==( Atom &B );
    // Makes a thermal displacement
    void ThermalDisplacement( int AtomIndex, Random &RND )
    {
      double XD, YD, ZD;
      double rTAmp=Thermal.TAmp[AtomIndex];

      XD = rTAmp * Thermal.Table[ (int) (THERMAL_NTBL*RND.Generate()) ];
//      if( RND.Generate() < 0.5 ) XD = -XD;
      if( rand() < RAND_MAX >> 1 ) XD = -XD; // A litle bit faster

      YD = rTAmp * Thermal.Table[ (int) (THERMAL_NTBL*RND.Generate()) ];
//      if( RND.Generate() < 0.5 ) YD = -YD;
      if( rand() < RAND_MAX >> 1 ) YD = -YD;
      
      ZD = rTAmp * Thermal.Table[ (int) (THERMAL_NTBL*RND.Generate()) ];
//      if( RND.Generate() < 0.5 ) ZD = -ZD;
      if( rand() < RAND_MAX >> 1 ) ZD = -ZD;
      
      R.X+=XD;
      R.Y+=YD;
      R.Z+=YD;
    };

    // Distance function in any direction (SIGNED VALUE)
    friend double ByDistance( Atom &A, Vector &V );
    // Deep function in any direction (SIGNED VALUE)
    friend double ByDeep( Atom &A, Vector &V );
    // Print method for the data (For XMOL)
    int PrintXYZ(std::ostream & Out );
    // Output streams overloaded : prototype
    friend std::ostream &operator<<(std::ostream &Out, Atom A);
};

 // Distance function in any direction (SIGNED VALUE)
    inline double ByDistance( Atom &A, Vector &V )
    {
      Vector R;
      R.X = A.R.X-V.X;
      R.Y = A.R.Y-V.Y;
      R.Z = A.R.Z-V.Z;
      return sqrt( R.X*R.X + R.Y*R.Y + R.Z*R.Z );
    };
    // Deep function in any direction (SIGNED VALUE)
    inline double ByDeep( Atom &A, Vector &V )
    {
      V = Unitary( V );
      return( V * A.R );
    };

// Define Projectile class ////////////////////////////////////////////////

class Projectile : public Atom
{
 public:
    Vector   Abs;       // Absolute coordinates vector
    Vector   InitPos;   // Initial absolute coordinates of projectile
    Vector   Dir;       // Direction (unitary vector)
    double     Energy;    // Energy in eV
    double     Ghi;       // Forward interaction distance
    double     S;         // Impact parameter in the last collision
    AtomType Status;    // Status of the atom
    double     Einel;     // Inelastic losses
    short int Zone;     // Splitting zone passed (=0) for rare event
    short int Zone2;    // Splitting superf. zone
    double     Weigth;    // Statistical weigth for rare event algorithm
   unsigned long Iterations; // Number of iterations done
   unsigned long Fly;   // Number of no collisions done
    double     LocalPath; // Partial distance done by projectile

    double     InitialTime;  // Initial time related to cascade start
    double     FinalTime;    // Final time. Time elapsed until stop the projectile
    double     mean_r;	   // Mean radius of the trajectory. To solve dinamically the 
			   //   cross section of the trajectory
    double     max_r;
    
    // Constructor
    Projectile( Atom Default   = Atom(0),
            Vector   Reference = NullVector,
            double     G         = 0.1,
            Vector   DirC      = XVector,
            double     EnergyC   = 0.0,
            AtomType StatusC   = LATTICE_ATOM,
            double     EInel     = 0.0
          ):Atom(Default)
    {
     S   = 0;
     Ghi     = G;           // Initial Ghi
     Abs     = Reference;   // At first time the absolute position
     InitPos = Reference;   //  is equal to original coordinates
     Energy  = EnergyC;
     Status  = StatusC;
     Einel   = EInel;
     Zone    = 0;
     Zone2   = 0;
     Weigth  = 1.0;
     if(DirC != NullVector) Dir = Unitary(DirC);
     
     Iterations = 0;
     Fly        = 0;
     LocalPath  = 0.0;
    };

    // For sorting methods
    friend double ByEnergy( Projectile &P, Vector &V ) { return( P.Energy );};
    // For sorting methods
    friend double ByEnergyAndType( Projectile &P, Vector &V )
    {
     if( P.Status == PROJECTILE ) return( P.Energy + 1e10 );
     else                 return( P.Energy );
    };
    // Deep function in any direction (SIGNED VALUE)
    double Deep( Vector Direction );
    // Print method for the data (For XMOL)
    int PrintXYZ(std::ostream & Out, int Option );
    // Output streams overloaded : prototype
    friend std::ostream &operator<<(std::ostream &Out, Projectile A);
};
}
#endif  // _INCLUDE_EDT1D_H_
