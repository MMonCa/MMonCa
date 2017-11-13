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
// CRYSTAL.H
//
///////////////////////////////////////////////////////////////////////////
#ifndef _INCLUDE_CRYSTAL_H_
#define _INCLUDE_CRYSTAL_H_

#include "defs.h"
#include "list.h"
#include "atom.h"
#include "thermal.h"
#include "edt1d.h"
#include "edt3d.h"

#include <assert.h>
#include <iostream>
#include <math.h>
#include <fstream>
#include <time.h>
namespace BCA {
extern List<AtomDefinition> Atoms;
extern Vector NullVector;
extern Orientacion OR;
extern double RI;

enum AxeType { AXE_X = 1, AXE_Y, AXE_Z };

enum CenteredType       // The same sorting as MW
     { PRIMITIVE = 1, BODY_CENTERED, A_END_CENTERED, B_END_CENTERED,
    C_END_CENTERED, FACE_CENTERED };

// Class Xtal definition //////////////////////////////////////////////////

class XtalDefinition
{
 public:
   Vector   Origin;
   int      Type, Type2;
   double     X;
   CenteredType Center;
   double     BindingEnergy, BindingEnergy2;

   void Set( Vector O, int T = 1 , CenteredType C = FACE_CENTERED, 
        double  BE = 15.0 /*eV*/ )
    {
      Origin = O;
      Type   = Type2 = T;
      Center = C;
      BindingEnergy = BindingEnergy2 = BE;
      X      = 1.0;
    };

   friend std::ostream &operator<<(std::ostream &Out, XtalDefinition A);
   friend std::istream &operator>>(std::istream &In,  XtalDefinition &V );
};

///////////////////////////////////////////////////////////////////////////
//
//  Define the Crystal CLASS, a child class from List<Atom>
//
class Crystal : public List<Atom>
{
  private:
    int ncellsx, ncellsy, ncellsz;
    int NeighIndex;

  // Public data
  public:
    Vector     LatticeParameter;
    double       XTalSize;
    double       InteractionRadius;
    double       GhiLimit;
    double       SimultaneousDistance;// (A)
    double       IR2;
    double       MaxIP;            // (A^2)

    double       CellVolume;       // (A)^3
    double       MeanAtomicRadius; // Its square (A)^2
    double       Dens;             // atoms/cm^3
    Atom Near;                   // Nearest atom to the projectle
                                 //  after Nearest procedure
                                 //  or a SearchAtom procedure
    DensityClass DENS;
    
    // Array of neighbors for each different lattice site
    List<Atom> Neighbors[NEIGHMAX];

  // Constructor
  Crystal( double XS = 0.5, double IR = 2.7155, Vector LP = NullVector,
       double GL = 0.1, double SD = 0.5, double IP = 0.5 )
  {
    LatticeParameter  = LP;        // Measured in Angstroms
    if( LatticeParameter == NullVector ) LatticeParameter(5.431,5.431,5.431);
    XTalSize             = XS;     //      in lattice units
    InteractionRadius    = IR;     //      in (A)
    GhiLimit             = GL;     //      in lattice units
    SimultaneousDistance = SD;     //      in A
    IR2                  = InteractionRadius*InteractionRadius;
    MaxIP                = IR2;
    NeighIndex           = -1;     // The first used index will be 0.
  };

  // Makes a Bravais' lattice. Displaced or not orthogonal cubic lattice.
  int MakeIt1(
          Vector       Origin,  // Origin
          CenteredType Type,    // Bravais' cell type
          int      AtomType1,   // Atom to fill-in.
          int      AtomType2,
          double X
        );
  int MakeIt1( XtalDefinition &X )
       { return MakeIt1( X.Origin, X.Center, X.Type, X.Type2, X.X ); };

  // Makes the transformations to obtain the desired lattice.
  int MakeIt2(                // Lattice parameters for each axe
          double a,             // X axe (A)
          double b,             // Y axe (A)
          double c,             // Z axe (A)
                              // Angles (Degrees) between axes
          double Alpha,         // b^c angle
          double Beta,          // a^c angle
          double Gamma          // a^b angle
         );

  void CreateEDT( char *EDTFile, EDT1D Dens1D[N_ATOMS+1] );
  int Update( double XS = 0.5, double IR = 2.7155, Vector LP = NullVector,
             double GL = 0.1, double SD = 0.5, double IP = 0.5 );
  int MakeNeighbors( Vector Origin );
  int MakeNeighbors( XtalDefinition &X ) { return MakeNeighbors( X.Origin ); };
  int Rotate( double Angle, AxeType Axe, Vector Reference );
  int Rotate( double M[3][3] );
  int UnRotate(  );
  int ShowToXMol(const std::string &FileName = "XMOLFile", Vector VP = NullVector );
  friend std::ostream &operator<<(std::ostream &Out, Crystal X );

  Atom Nearest( double (*fMatch)(Atom &D, Vector &V), Vector R )
   {
    // overload of primitive function to store the atom
    Near = List<Atom>::Nearest( fMatch, R );
    return Near;
   };

  // Search targets from the crystal    // Return the number of targets
  int SearchTargets( List<Projectile> &Targets,
             Projectile &CurrentProjectile,
             double   &IP2Min,
             Random &RND,
             double   FrontSurface,
             double   BackSurface,
             List<Atom> &Lista2 
           );
};

///////////////////////////////////////////////////////////////////////////

// Print a List of ATOMS in XYZ format to represent with XMOL /////////////
int ShowToXMol( List<Atom> LIST ,char* FileName );

// Print a List of PROJECTILES in XYZ format to represent with XMOL ///////
int ShowToXMol( List<Projectile> LIST ,char* FileName, int Option = 0 );
}
#endif  // _INCLUDE_CRYSTAL_H_
