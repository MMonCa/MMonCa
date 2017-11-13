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
// ATOM.CPP
//
///////////////////////////////////////////////////////////////////////////
#include "atom.h"
#include "atomdef.h"

using std::ostream;
using std::istream;
using std::endl;

namespace BCA {

extern Orientacion OR;

char AtomTypeStr[][20] =
    {
        "LATTICE ATOM",
        "PROJECTILE",
        "VACANCY",
        "INTERSTITIAL",
        "SECONDARY_PROJ",
        "SPUTTERED",
        "BACK SCATTERED",
        "IMPLANTED ION",
        "OUTSIDE ION",
        "ION NOT FOLLOWED",
        "NOT FOLLOWED",
        "OUTSIDE ATOM",
        "PHONON EXCITED"
    };

// Define Atom definition class ///////////////////////////////////////////
AtomDefinition::AtomDefinition(const std::string &TC, double WC, int ZC, double TD )
    {
      assert( WC>0.0 && ZC>0 && ZC<103 );

      W  = WC;
      Z  = ZC;
      SetName( TC );
      TAmp   = 0.0;
      TDebye = TD;
    }

// Compare two atoms
int AtomDefinition::operator==( AtomDefinition &B )
{
      if(Type == B.Type && W == B.W && Z == B.Z)
    	  return 1;
      else
    	  return 0;
}

// Change the chemical symbol
void AtomDefinition::SetName(const std::string &TC)
{
      assert( TC[0]>='A' && TC[0]<='Z' );
      Type = TC;
}

// Get the chemical symbol
const char * AtomDefinition::GetName()
{
      static std::string TC;
      TC = Type;
      return TC.c_str();
}

// Output stream overloaded
ostream &operator<<( ostream &Out, AtomDefinition A)
 {
  Out << "[" << A.Type << "]\tW =";
  Out.width(10);
  Out << A.W << " amu, Z =";
  Out.width(4);
  Out << (A.Z&255) << ", TDebye = ";
  Out.width(5);
  Out << A.TDebye << " K, TAmp = ";
  Out.width(10);
  Out << A.TAmp << " A" << endl;
  return Out;
 }

// Input stream overloaded
istream &operator>>( istream &In, AtomDefinition &A )
 {
 std::string TC;

  In >> TC;
  A.SetName( TC );
  In >> A.Z;
  In >> A.W;
  In >> A.TDebye;
  return In;
 }

// Define the Atom class //////////////////////////////////////////////////

// Constructor
Atom::Atom( int IC, int L, Vector RC, int IC2, double xx )
{
      Index  = (unsigned char) IC;
      LS     = (unsigned char) L;
      R      = RC;
      Index1 = Index;
      Index2 = (unsigned char) IC2;
      X      = xx;
}

// Compare two atoms
int Atom::operator==( Atom &B )
{
      if( ( R == B.R ) && ( Index == B.Index ) ) return 1;
      else                       return 0;
}

// Print method for the data (For XMOL)
int Atom::PrintXYZ( ostream & Out )
{
    //  Out  << "Si"; // << (int) Index;
      Out << " ";
      Out.width(10);
      Out.precision(6);
      Out.setf(std::ios::showpoint);
      Out  << OR.ProyectaX(R);
      Out << " ";
      Out.width(10);
      Out.precision(6);
      Out.setf(std::ios::showpoint);
      Out  << OR.ProyectaY(R);
      Out << " ";
      Out.width(10);
      Out.precision(6);
      Out.setf(std::ios::showpoint);
      Out  << OR.ProyectaZ(R) << endl;
      return Index;
}

// Output stream overloaded
ostream &operator<<( ostream &Out, Atom A)
 {
  Out << "[" << (int)A.Index << "]\tR   " << A.R << endl;
  return Out;
 }

// Define Projectile class ////////////////////////////////////////////////

// Deep function in any direction (SIGNED VALUE)
double Projectile::Deep( Vector Direction )
{
     assert( Direction != NullVector );

     Direction = Unitary( Direction );
     return( Direction * Abs );
}

// Print method for the data (For XMOL)
int Projectile::PrintXYZ( ostream & Out, int Option )
{
     if ( Status == INTERSTITIAL )
     {
       Out << "Si" << Index << "\t";        // Print interstitial
       Out.width(WIDTH);
       Out  << OR.ProyectaX(Abs) << "\t"
        << OR.ProyectaY(Abs) << "\t"
        << OR.ProyectaZ(Abs) << endl;
       if (Option!=0)
       {
       Out << "N" << Index << "\t";         // Print vacancy
       Out.width(WIDTH);
       Out  << OR.ProyectaX(InitPos) << "\t"
        << OR.ProyectaY(InitPos) << "\t"
        << OR.ProyectaZ(InitPos) << endl;
       };
     }
     else if ( Status == PROJECTILE )
     {
       Out << "Si" << Index << "\t";        // Print interstitial
       Out.width(WIDTH);
       Out  << OR.ProyectaX(Abs) << "\t"
        << OR.ProyectaY(Abs) << "\t"
        << OR.ProyectaZ(Abs) << endl;
       if (Option!=0)
       {
       Out << "N" << Index << "\t";         // Print vacancy
       Out.width(WIDTH);
       Out  << OR.ProyectaX(InitPos) << "\t"
        << OR.ProyectaY(InitPos) << "\t"
        << OR.ProyectaZ(InitPos) << endl;
       };
     }
     else
     {
       Out << "C" << Index << "\t";         // Print ion
       Out.width(WIDTH);
       Out  << OR.ProyectaX(Abs) << "\t"
        << OR.ProyectaY(Abs) << "\t"
        << OR.ProyectaZ(Abs) << endl;
     };
     return (Index);
}

// Output stream overloaded
ostream &operator<<( ostream &Out, Projectile A)
{
  Out << (Atom) A;
  Out << "Energy: " << A.Energy  << " (eV)\n";
  Out << "Dir "     << A.Dir     << endl;
  Out << "Abs "     << A.Abs     << endl;
  Out << "Org "     << A.InitPos << endl;
  Out << " Ghi = "  << A.Ghi     << endl;
  Out << AtomTypeStr[A.Status] << endl;  
  Out << endl << endl;
  return Out;
}

}
///////////////////////////////////////////////////////////////////////////
